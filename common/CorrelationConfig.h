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
 * File: written by Bootstrap_Correlations, looked for in this order:
 *   1. the path in EASYBJETS_CORRELATIONS
 *   2. ../Run3Plots/Bootstrap/correlations.txt (../Run2Plots for run2), i.e.
 *      next to the plots of the channel, one file per run era
 *   3. correlations.txt in the working directory (older location)
 * One line per (domain, channel, WP, bin):
 *
 *   mass DB 0p46 bJER10v01:pt200to300            rho_mu_dd=0.98 rho_mu_mm=0.95 ...
 *   resp DB 0p46 bJER10v01:inclusive             rho_mu_dd=0.97 ... rho_sig_tt=0.99 ...
 *   resp DB 0p46 bJER10v01:200to300GeV           (pT bin, no mass cut)
 *   resp DB 0p46 bJER10v01:200to300GeV_50-80GeV  (pT bin and mass bin)
 *   DB 0p46 pt200to300 rho_mu_dd=0.98            (older format: no domain = mass,
 *                                                 no correction prefix)
 * The bin label starts with the correction ("bJER10v01:") because each
 * correction has its own correlation with the nominal. CorrFor() with a
 * correction name looks for that label first and then for the label without
 * the prefix (files written by the older bootstrap).
 *
 * Lines starting with # are ignored. Keys not given on a line, and bins not in
 * the file, fall back to the default matrix: rho = 1 between correction and
 * nominal of the same quantity, 0 between data and MC and between MC reco and
 * truth widths.
 * ==========================================================================
 * @version:  1.0
 * ==========================================================================
 */

#ifndef EASYBJETS_CORRELATION_CONFIG_H
#define EASYBJETS_CORRELATION_CONFIG_H

#include "RunEra.h"  // RunConfig::IsRun3, for the default file location

#include <algorithm>
#include <cctype>
#include <cmath>
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

    // Response only: MC truth width sigma_truth (pT_truth / pT_ref).
    // Correction and nominal of the truth width: the corrections change the
    // reco jet, not the truth jet, so this is close to 1.
    double rho_sig_tt     = 1.0;
    // MC reco width vs MC truth width (same MC events). Default 0.
    double rho_sig_mv_tv  = 0.0;
    double rho_sig_mv_tn  = 0.0;
    double rho_sig_mn_tv  = 0.0;
    double rho_sig_mn_tn  = 0.0;
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
// Default location of correlations.txt: in the plot folder of the channel,
// one file per run era (the macros run from <channel>/macros).
inline std::string DefaultCorrelationsPath() {
    return std::string(RunConfig::IsRun3() ? "../Run3Plots" : "../Run2Plots") + "/Bootstrap/correlations.txt";
}

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
    std::string path = envPath ? std::string(envPath) : DefaultCorrelationsPath();
    if (!envPath && !std::ifstream(path).good() && std::ifstream("correlations.txt").good()) {
        std::cout << "[Correlations] " << path << " not found, using ./correlations.txt"
                  << " (older location)\n";
        path = "correlations.txt";
    }

    std::ifstream in(path);
    if (!in.is_open()) {
        if (method == CorrMethod::Bootstrap) {
            std::cerr << "[Correlations] WARN: EASYBJETS_CORR_METHOD=bootstrap but '"
                      << path << "' not found — run Bootstrap_Correlations.cxx first"
                      << " (same EASYBJETS_RUN). Using the default rho.\n";
        } else {
            std::cout << "[Correlations] no file at '" << path
                      << "' — default rho used (1 between correction and nominal, 0 otherwise)\n";
        }
        return cache;
    }

    std::string line;
    int nLoaded = 0, nMass = 0, nResp = 0, nOld = 0;
    while (std::getline(in, line)) {
        // Settings of the bootstrap run that wrote the file, shown once
        if (line.rfind("# Settings:", 0) == 0) std::cout << "[Correlations] " << line.substr(2) << "\n";
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
            else if (k == "rho_sig_tt")      C.rho_sig_tt      = v;
            else if (k == "rho_sig_mv_tv")   C.rho_sig_mv_tv   = v;
            else if (k == "rho_sig_mv_tn")   C.rho_sig_mv_tn   = v;
            else if (k == "rho_sig_mn_tv")   C.rho_sig_mn_tv   = v;
            else if (k == "rho_sig_mn_tn")   C.rho_sig_mn_tn   = v;
        }
        cache[CorrKey(domain, channel, wp, ptbin)] = C;
        ++nLoaded;
        if (ptbin.find(':') == std::string::npos) ++nOld;
        if (domain == CorrDomain::Mass) ++nMass; else ++nResp;
    }

    std::cout << "[Correlations] loaded " << nLoaded
              << " entries from " << path
              << " (mass: " << nMass << ", resp: " << nResp << ")\n";
    if (nOld > 0)
        std::cerr << "[Correlations] WARN: " << nOld << " lines have the old format (no \"<correction>:\""
                  << " in the bin label). Files from the older bootstrap gave rho ~ 0 (independent"
                  << " weights); re-run Bootstrap_Correlations to replace them.\n";
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

// Bin label of the response fits: the pT-bin label, plus "_<mass bin>" when
// the fit is split in mass bins ("MassInclusive" or empty = no mass split).
// Used by Bootstrap_Correlations when writing and by the ratio macros when
// reading, so both sides build exactly the same label.
inline std::string RespBinLabel(const std::string& ptLabel, const std::string& massDir = "") {
    if (massDir.empty() || massDir == "MassInclusive") return ptLabel;
    return ptLabel + "_" + massDir;
}

// Correlations of one correction (e.g. "bJER10v01", the variant path_tag):
// looks for "<correction>:<label>" first, then for "<label>" (older files),
// and returns the default matrix if neither is in the file.
inline CorrelationMatrix CorrFor(CorrDomain domain, const std::string& channel,
                                  const std::string& wp, const std::string& label,
                                  const std::string& correction) {
    const auto& m = Correlations();
    if (!correction.empty()) {
        auto it = m.find(CorrKey(domain, channel, wp, correction + ":" + label));
        if (it != m.end()) return it->second;
    }
    auto it = m.find(CorrKey(domain, channel, wp, label));
    return (it != m.end()) ? it->second : CorrelationMatrix{};
}

// Full-covariance relative error of a resolution double ratio, c^T Sigma c over
// the six fitted widths
//   x = (s_d,corr, s_d,nom, s_m,corr, s_m,nom, s_t,corr, s_t,nom)
// (data, combined MC and truth widths, correction and nominal), with
//   u = s_d^2 - s_t^2 = sigma_insitu^2,   v = s_m^2 - s_t^2 = (sigma^MC_insitu)^2.
// The c are the derivatives of ln(DR):
//   DR(sigma_insitu): d ln/ds_d = s_d/u,  d ln/ds_t = -s_t/u           (no s_m)
//   DR(R_JER):        d ln/ds_d = s_d/u,  d ln/ds_m = -s_m/v,
//                     d ln/ds_t = s_t (1/v - 1/u)
// with opposite signs for the nominal. The rho come from C (data-truth pairs
// are 0: different samples). Returns 0 if a difference of squares is <= 0.
inline double ResolutionDRRelErrFull(bool isRJER,
                                     double sd_c, double sd_c_err, double sd_n, double sd_n_err,
                                     double sm_c, double sm_c_err, double sm_n, double sm_n_err,
                                     double st_c, double st_c_err, double st_n, double st_n_err,
                                     const CorrelationMatrix& C) {
    const double u_c = sd_c * sd_c - st_c * st_c, u_n = sd_n * sd_n - st_n * st_n;
    const double v_c = sm_c * sm_c - st_c * st_c, v_n = sm_n * sm_n - st_n * st_n;
    if (u_c <= 0 || u_n <= 0) return 0.0;
    if (isRJER && (v_c <= 0 || v_n <= 0)) return 0.0;

    double c[6];
    c[0] = sd_c / u_c;
    c[1] = -sd_n / u_n;
    if (isRJER) {
        c[2] = -sm_c / v_c;
        c[3] = sm_n / v_n;
        c[4] = st_c * (1.0 / v_c - 1.0 / u_c);
        c[5] = -st_n * (1.0 / v_n - 1.0 / u_n);
    } else {
        c[2] = 0.0;
        c[3] = 0.0;
        c[4] = -st_c / u_c;
        c[5] = st_n / u_n;
    }
    const double e[6] = {sd_c_err, sd_n_err, sm_c_err, sm_n_err, st_c_err, st_n_err};

    double rho[6][6] = {};
    for (int i = 0; i < 6; ++i) rho[i][i] = 1.0;
    auto set = [&](int i, int j, double r) { rho[i][j] = rho[j][i] = r; };
    set(0, 1, C.rho_sig_dd);
    set(2, 3, C.rho_sig_mm);
    set(4, 5, C.rho_sig_tt);
    set(0, 2, C.rho_sig_dv_mv);
    set(0, 3, C.rho_sig_dv_mn);
    set(1, 2, C.rho_sig_dn_mv);
    set(1, 3, C.rho_sig_dn_mn);
    set(2, 4, C.rho_sig_mv_tv);
    set(2, 5, C.rho_sig_mv_tn);
    set(3, 4, C.rho_sig_mn_tv);
    set(3, 5, C.rho_sig_mn_tn);

    double rel2 = 0.0;
    for (int i = 0; i < 6; ++i)
        for (int j = 0; j < 6; ++j) rel2 += c[i] * c[j] * rho[i][j] * e[i] * e[j];
    return std::sqrt(std::max(0.0, rel2));
}

}  // namespace CorrelationConfig

#endif  // EASYBJETS_CORRELATION_CONFIG_H
