/**
 * ============================================================================
 * CorrelationConfig.h - correlations between variant and nominal fits
 * ============================================================================
 * ===========================================
 * @author:    Thiziri AMEZZA
 * @email:     thiziri.amezza@cern.ch
 * @date:      10 March 2026 
 * ===========================================
 *
 * The variant/nominal double ratios are built from fits run on the same
 * events, so the fitted mu and sigma are strongly correlated between the two.
 * The ratio macros support three treatments, selected with
 * EASYBJETS_CORR_METHOD:
 *   rho1      (default) assume rho = 1 for every shared-event pair
 *   file      read rho values from a text file
 *   boot      same as file, but the file is expected to come from the
 *             bootstrap macro (only the warning message differs)
 *
 * File format: correlations.txt in the working directory, or the path in EASYBJETS_CORRELATIONS. 
 * One line per (domain, channel, WP, pT bin):
 *
 *   mass DB  0p46 pt200to300 rho_mu_dd=0.98 rho_mu_mm=0.95 rho_sig_dd=0.92 rho_sig_mm=0.88
 *   resp DB  0p46 pt200to300 rho_mu_dd=0.97 rho_mu_mm=0.94
 *   DB 0p46 pt200to300 rho_mu_dd=0.98      (no domain = mass, older format)
 *
 * Lines starting with # are ignored. Keys not given on a line, and bins not in the file at all, 
 * fall back to rho = 1 for shared pairs and 0 for cross pairs, so the "full correlation" result 
 * equals the rho = 1 one wherever the file has uncovered bin.
 * ==========================================================================
 * @version:  1.0
 * ==========================================================================
 */

#ifndef EASYBJETS_CORRELATION_CONFIG_H
#define EASYBJETS_CORRELATION_CONFIG_H

#include <cctype>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>

namespace CorrelationConfig {

// ----------------------------------------------------------------------------
// Method dispatch
// ----------------------------------------------------------------------------
enum class CorrMethod { Rho1, File, Bootstrap };

inline CorrMethod CurrentCorrMethod() {
    static CorrMethod cached = [](){
        const char* env = std::getenv("EASYBJETS_CORR_METHOD");
        if (!env || !*env) return CorrMethod::Rho1;
        std::string s(env);
        for (auto& c : s) c = (char)std::tolower((unsigned char)c);
        if (s == "rho1" || s == "1" || s == "default")  return CorrMethod::Rho1;
        if (s == "file" || s == "txt")                   return CorrMethod::File;
        if (s == "bootstrap" || s == "boot")             return CorrMethod::Bootstrap;
        std::cerr << "[Correlations] WARN: unknown EASYBJETS_CORR_METHOD='"
                  << env << "', using default (rho1)\n";
        return CorrMethod::Rho1;
    }();
    return cached;
}

inline std::string CorrMethodTag() {
    switch (CurrentCorrMethod()) {
        case CorrMethod::Rho1:      return "rho1";
        case CorrMethod::File:      return "file";
        case CorrMethod::Bootstrap: return "boot";
    }
    return "unknown";
}

// ----------------------------------------------------------------------------
// Entries of the 4x4 correlation matrix between (data_var, data_nom, mc_var,
// mc_nom), once for mu and once for sigma. Defaults are the rho = 1 model.
// Applies to either mass-peak OR response fits — same structure, same fields.
// The file parser uses a domain prefix ("mass_" / "resp_") to keep them
// separate; see Correlations() below.
// ----------------------------------------------------------------------------
struct CorrelationMatrix {
    // Shared events, different calibration: the large ones
    double rho_mu_dd   = 1.0;   // data mu shared
    double rho_mu_mm   = 1.0;   // MC mu shared
    double rho_sig_dd  = 1.0;   // data sigma shared
    double rho_sig_mm  = 1.0;   // MC sigma shared

    // data vs MC: different samples (typically ~0; different events)
    double rho_mu_dv_mv   = 0.0;
    double rho_mu_dv_mn   = 0.0;
    double rho_mu_dn_mv   = 0.0;
    double rho_mu_dn_mn   = 0.0;
    double rho_sig_dv_mv  = 0.0;
    double rho_sig_dv_mn  = 0.0;
    double rho_sig_dn_mv  = 0.0;
    double rho_sig_dn_mn  = 0.0;
};

// Observable domain: Mass peak fits vs response R_DB/R_MJB fits have
// independent correlation matrices (different fits on different histograms,
// can have different Rho even though the underlying events are the same).
enum class CorrDomain { Mass, Response };

inline std::string DomainTag(CorrDomain d) {
    return (d == CorrDomain::Response) ? "resp" : "mass";
}

// Key = "<domain>_<channel>_<wp>_<ptbin>", e.g.
//   mass_DB_0p46_pt200to300
//   resp_DB_0p46_pt200to300
inline std::string CorrKey(CorrDomain d, const std::string& channel,
                            const std::string& wp, const std::string& ptbin) {
    return DomainTag(d) + "_" + channel + "_" + wp + "_" + ptbin;
}

inline std::string CorrKey(const std::string& channel,
                            const std::string& wp, const std::string& ptbin) {
    return CorrKey(CorrDomain::Mass, channel, wp, ptbin);
}

// ----------------------------------------------------------------------------
// File loader. Parses correlations.txt and caches per-process.
// Empty map in rho1 mode or when the file is missing.
// ----------------------------------------------------------------------------
inline const std::map<std::string, CorrelationMatrix>& Correlations() {
    static std::map<std::string, CorrelationMatrix> cache;
    static bool loaded = false;
    if (loaded) return cache;
    loaded = true;

    CorrMethod method = CurrentCorrMethod();
    if (method == CorrMethod::Rho1) {
        return cache;
    }

    const char* envPath = std::getenv("EASYBJETS_CORRELATIONS");
    std::string path = envPath ? std::string(envPath) : std::string("correlations.txt");

    std::ifstream in(path);
    if (!in.is_open()) {
        if (method == CorrMethod::Bootstrap) {
            std::cerr << "[Correlations] WARN: EASYBJETS_CORR_METHOD=bootstrap but '"
                      << path << "' not found — did you run "
                      << "Bootstrap_Correlations_*.cxx? Falling back to ρ=1.\n";
        } else {
            std::cout << "[Correlations] no file at '" << path
                      << "' — all ρ default to 1.0 (FullCorr = CorrAware)\n";
        }
        return cache;
    }

    std::string line;
    int nLoaded = 0, nMass = 0, nResp = 0;
    while (std::getline(in, line)) {
        size_t p = 0;
        while (p < line.size() && std::isspace((unsigned char)line[p])) ++p;
        if (p >= line.size() || line[p] == '#') continue;

        std::istringstream iss(line.substr(p));
        std::string tok1;
        if (!(iss >> tok1)) continue;

        CorrDomain domain = CorrDomain::Mass;  
        std::string channel;
        if (tok1 == "mass") { domain = CorrDomain::Mass;     if (!(iss >> channel)) continue; }
        else if (tok1 == "resp") { domain = CorrDomain::Response; if (!(iss >> channel)) continue; }
        else { channel = tok1; }  

        std::string wp, ptbin;
        if (!(iss >> wp >> ptbin)) continue;

        CorrelationMatrix C;
        std::string tok;
        while (iss >> tok) {
            auto eq = tok.find('=');
            if (eq == std::string::npos) continue;
            std::string k = tok.substr(0, eq);
            double v;
            try { v = std::stod(tok.substr(eq + 1)); }
            catch (...) { continue; }
            if      (k == "rho_mu_dd")       C.rho_mu_dd       = v;
            else if (k == "rho_mu_mm")       C.rho_mu_mm       = v;
            else if (k == "rho_sig_dd")      C.rho_sig_dd      = v;
            else if (k == "rho_sig_mm")      C.rho_sig_mm      = v;
            else if (k == "rho_mu_dv_mv")    C.rho_mu_dv_mv    = v;
            else if (k == "rho_mu_dv_mn")    C.rho_mu_dv_mn    = v;
            else if (k == "rho_mu_dn_mv")    C.rho_mu_dn_mv    = v;
            else if (k == "rho_mu_dn_mn")    C.rho_mu_dn_mn    = v;
            else if (k == "rho_sig_dv_mv")   C.rho_sig_dv_mv   = v;
            else if (k == "rho_sig_dv_mn")   C.rho_sig_dv_mn   = v;
            else if (k == "rho_sig_dn_mv")   C.rho_sig_dn_mv   = v;
            else if (k == "rho_sig_dn_mn")   C.rho_sig_dn_mn   = v;
        }
        cache[CorrKey(domain, channel, wp, ptbin)] = C;
        ++nLoaded;
        if (domain == CorrDomain::Mass) ++nMass; else ++nResp;
    }

    std::cout << "[Correlations] loaded " << nLoaded
              << " entries from " << path
              << " (mass: " << nMass << ", resp: " << nResp << ")\n";
    return cache;
}

// Look up correlations for a given (domain, channel, wp, ptbin). 
// Missing keys return identity (rho1 model) for anything not in the file.
inline CorrelationMatrix CorrFor(CorrDomain domain,
                                  const std::string& channel,
                                  const std::string& wp,
                                  const std::string& ptbin) {
    const auto& m = Correlations();
    auto it = m.find(CorrKey(domain, channel, wp, ptbin));
    return (it != m.end()) ? it->second : CorrelationMatrix{};
}

inline CorrelationMatrix CorrFor(const std::string& channel,
                                  const std::string& wp,
                                  const std::string& ptbin) {
    return CorrFor(CorrDomain::Mass, channel, wp, ptbin);
}

}  // namespace CorrelationConfig

#endif  // EASYBJETS_CORRELATION_CONFIG_H
