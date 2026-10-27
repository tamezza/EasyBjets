/**
 * ============================================================================
 *  BootstrapErrors.h - double-ratio uncertainties measured by the bootstrap
 * ============================================================================
 * ===========================================
 * @author:    Thiziri AMEZZA
 * @email:     thiziri.amezza@cern.ch
 * @date:      26 October 2026
 * ===========================================
 *
  * Gives the summaries the "bootstrap" error of a double ratio, as written by
 * Bootstrap_Analyze in dr_errors.csv: the double ratio is computed in each of
 * the bootstrap replicas, and its error is half the width of the range that
 * holds the central 68% of them (16th to 84th percentile; the same as one
 * sigma for a Gaussian, but not pulled by a few replicas whose fit failed).
 *
 * The other three treatments propagate the fit errors with a formula. This
 * one does not use the fit errors: it measures how much the final double
 * ratio moves when the data fluctuate, so it includes all the correlations
 * and also the fluctuation of the subtracted background, which the Mass_Fit
 * errors treat as exactly known.
 *
 * Files read (both channels, the era of EASYBJETS_RUN), relative to the
 * <channel>/macros folder the summary runs from:
 * 
 *   ../../ZbbGamma/Run3Plots/Bootstrap/analysis/dr_errors.csv
 *   ../../ZbbJets/Run3Plots/Bootstrap/analysis/dr_errors.csv
 * 
 * (Run2Plots for run2). They are written by Bootstrap_Correlations.cxx then
 * Bootstrap_Analyze.cxx in each channel.
 *
 * A point with no bootstrap value (e.g. MuonCorr, which the bootstrap cannot
 * measure, or a channel not run) keeps its full-correlation error; the
 * summary prints how many points did.
 *
 * ==========================================================================
 * @version:  1.0
 * ==========================================================================
 */

#ifndef EASYBJETS_BOOTSTRAP_ERRORS_H
#define EASYBJETS_BOOTSTRAP_ERRORS_H

#include "RunEra.h"      
#include "CommonUtils.h" 

#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace BootstrapErrors {

// Correction name without the channel prefix: DB_bJER10v01 -> bJER10v01.
inline std::string CorrTag(const std::string& variant) {
    for (const std::string p : {"DB_", "MJB_"})
        if (variant.rfind(p, 0) == 0) return variant.substr(p.size());
    return variant;
}

// Bin label of the mass fits: pt200to450, pt450to750, ...
inline std::string MassBinLabel(double ptLo, double ptHi) {
    char b[64];
    std::snprintf(b, sizeof(b), "pt%dto%d", (int)std::lround(ptLo), (int)std::lround(ptHi));
    return b;
}

// Bin label of the response fits: "inclusive", "200to300GeV" (all masses) or
// "200to300GeV_50-80GeV" (pT x mass bin).
inline std::string RespBinLabel(const std::string& channel, int ptBinIndex, const std::string& massBin) {
    if (massBin == "inclusive") return "inclusive";
    const std::string pt = (channel == "MJB") ? BinConfig::PtBinFileLabel_3(ptBinIndex)
                                              : BinConfig::PtBinFileLabel_DB(ptBinIndex);
    if (massBin.empty() || massBin == "MassInclusive") return pt;
    return pt + "_" + massBin;
}

// Relative bootstrap errors of both channels, keyed by
// "<channel>|<dr>|<correction>|<wp>|<bin>", read once per session.
inline const std::map<std::string, double>& Table() {
    static std::map<std::string, double> t;
    static bool loaded = false;
    if (loaded) return t;
    loaded = true;
    const std::string era = RunConfig::IsRun3() ? "Run3Plots" : "Run2Plots";
    for (const std::string dir : {"ZbbGamma", "ZbbJets"}) {
        const std::string path = "../../" + dir + "/" + era + "/Bootstrap/analysis/dr_errors.csv";
        std::ifstream in(path);
        if (!in.is_open()) {
            std::cout << "[BootstrapErrors] no " << path << " (run Bootstrap_Correlations and"
                      << " Bootstrap_Analyze in " << dir << ")\n";
            continue;
        }
        std::string line;
        std::vector<std::string> head;
        int n = 0;
        while (std::getline(in, line)) {
            std::vector<std::string> tok;
            std::stringstream ss(line);
            std::string x;
            while (std::getline(ss, x, ',')) tok.push_back(x);
            if (head.empty()) {
                head = tok;
                continue;
            }
            auto col = [&](const std::string& name) -> std::string {
                for (size_t i = 0; i < head.size() && i < tok.size(); ++i)
                    if (head[i] == name) return tok[i];
                return "";
            };
            const std::string rel = col("rel_bootstrap68");
            if (rel.empty()) continue;
            const std::string key = col("channel") + "|" + col("dr") + "|" + CorrTag(col("correction")) + "|"
                                    + col("wp") + "|" + col("bin");
            t[key] = std::atof(rel.c_str());
            ++n;
        }
        std::cout << "[BootstrapErrors] " << n << " double-ratio errors read from " << path << "\n";
    }
    return t;
}

// Number of points that found / did not find a bootstrap value.
inline std::pair<int, int>& Counts() {
    static std::pair<int, int> c{0, 0};
    return c;
}

// Relative bootstrap error of one double ratio; NaN if not measured.
//   dr: bJMDP, bJMDR, bJES, sigma_insitu or bJER
inline double RelErr(const std::string& channel, const std::string& dr, const std::string& variant,
                     const std::string& wp, const std::string& bin) {
    const auto& t = Table();
    const auto it = t.find(channel + "|" + dr + "|" + CorrTag(variant) + "|" + wp + "|" + bin);
    if (it == t.end() || !std::isfinite(it->second) || it->second <= 0) {
        ++Counts().second;
        return std::numeric_limits<double>::quiet_NaN();
    }
    ++Counts().first;
    return it->second;
}

// One line at the end of a summary.
inline void PrintCounts(const std::string& ctx) {
    if (Counts().first + Counts().second == 0) return;
    std::cout << "[" << ctx << "] bootstrap errors: " << Counts().first << " points from dr_errors.csv, "
              << Counts().second << " without (full-correlation error kept)\n";
}

}  // namespace BootstrapErrors

#endif  // EASYBJETS_BOOTSTRAP_ERRORS_H
