/**
 * ============================================================================
 *  Mass_Ratios_PeakWidth.h - variant/nominal double ratios of the mass fit
 * ============================================================================
 * ===========================================
 * @author:    Thiziri AMEZZA
 * @email:     thiziri.amezza@cern.ch
 * @date:      15 April 2026
 * ===========================================
 *
 * Data/MC ratios of the Z->bb mass peak and the double ratios that compare each
 * variant (correction) to the nominal one. For every variant, WP and pT bin:
 *   b-JMDP = mu_data / mu_MC          DR(JMDP) = JMDP_variant / JMDP_nominal
 *   b-JMDR = sigma_data / sigma_MC    DR(JMDR) = JMDR_variant / JMDR_nominal
 * b-JMDR is the mass analogue of the JER ratio of the response analysis:
 * dimensionless, > 1 when data is wider than MC.
 *
 * Each double ratio is given with three errors:
 *   uncorrelated   the four fit errors added in quadrature (conservative)
 *   corr-aware     variant and nominal are fitted on the same events, so their
 *                  fluctuations are taken as fully correlated (rho = 1)
 *   full-corr      rho read from correlations.txt; equal to corr-aware unless
 *                  EASYBJETS_CORR_METHOD=file (or bootstrap)
 *
 * How to run (from ZbbGamma/macros or ZbbJets/macros, after Mass_Fit.cxx):
 *   export EASYBJETS_RUN=run3          same as for Mass_Fit
 *   root 'Mass_Ratios_PeakWidth.cxx("inclusive")'   one pT bin
 *   root 'Mass_Ratios_PeakWidth.cxx("ptbins")'      pT bins of the channel
 *   root 'Mass_Ratios_PeakWidth.cxx("ALL")'         both
 * Use the same EASYBJETS_RUN, EASYBJETS_BINNING, EASYBJETS_FIT_MODEL and pT mode
 * as Mass_Fit, otherwise the fit_summary.txt files are not found.
 * Then run Mass_Summary_Combined.cxx.
 *
 * Input:  <PlotDir>/<variant>/Mass_Fit_<channel>/<wp>/<ptbin>_<binning>/fit_summary.txt
 * Output: <PlotDir>/Mass/Mass_Ratios_Width/Combined/WP_<wp>_<ptmode>_<binning>/
 *   Combined_Mass_Ratio_Summary.csv     input of Mass_Summary_Combined
 *   Combined_Mass_Ratio_Summary.txt     same numbers, readable
 *   Combined_{Mass_Scale,JMR}_Ratios*.pdf, Combined_2D_*_Map*.pdf,
 *   Combined_Mass_{Scale,Width}_Absolute.pdf, Summary/Summary_*.pdf
 *
 * Options (environment variables):
 *   EASYBJETS_RUN, EASYBJETS_BINNING, EASYBJETS_FIT_MODEL   must match Mass_Fit
 *   EASYBJETS_CORR_METHOD        rho1 (default) | file | bootstrap
 *   EASYBJETS_CORRELATIONS       path of the rho file (default ./correlations.txt)
 *   EASYBJETS_SUMMARY_VARIANTS   e.g. bJER10v01,MuonCorr: variants of the subset plot
 *
 * ==========================================================================
 * @version:  1.0
 * ==========================================================================
 */

#ifndef MASS_RATIOS_WIDTH_H
#define MASS_RATIOS_WIDTH_H

#include "VariantConfig.h"
#include "CommonUtils.h"
#include "BinningConfig.h"
#include "CorrelationConfig.h"
#include "Math/QuantFuncMathCore.h"
#include <algorithm>
#include <vector>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <map>
#include <string>
#include <cctype>
#include <functional>
#include <regex>
#include "TH2F.h"
#include "TPaletteAxis.h"
#include "TStyle.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TGaxis.h"

// Set once by MassRatios::Run so the helpers need no channel argument.
inline std::string& CurrentChannel() {
    static std::string ch = "DB";
    return ch;
}

inline std::string TopologyLabelFor(const std::string& channel) {
    return (channel == "MJB") ? std::string("Z(#rightarrow b#bar{b}) + jets")
                              : std::string("Z(#rightarrow b#bar{b}) + #gamma");
}

// ============================================================================
// PT BIN CONFIGURATION - uses PtBinConfig from VariantConfig.h
//   Mass_Ratios_MJB("inclusive")  -> single-bin mode
//   Mass_Ratios_MJB("ptbins")     -> multi-bin mode
// ============================================================================
constexpr int kMaxPtBins = 8; // PtBinDef internal cap


// ----------------------------------------------------------------------------
// Per-cell fit results read from Mass Analysis fit_summary.txt files
// ----------------------------------------------------------------------------
struct MassFitResults {
    double mean = 0;
    double mean_err = 0;
    double sigma = 0;
    double sigma_err = 0;
    double chi2_ndf = 0;
};


// Everything computed for one variant, per pT bin.
//   mean_ratio / sigma_ratio   DR(JMDP) / DR(JMDR) with the three errors
//   jms_* / jmr_*              absolute JMDP / JMDR for the variant and the
//                              nominal (the nominal is repeated in every entry)
//   mu_* / sig_*               raw fit components, written to the CSV
struct RatioResults {
    std::string label;
    int color = kBlue;
    int marker_style = 20;
    // Double-ratio quantities (variant / nominal)
    double mean_ratio[kMaxPtBins] = {};
    double mean_ratio_err[kMaxPtBins] = {};
    double mean_ratio_err_corr[kMaxPtBins] = {};
    double mean_ratio_err_full[kMaxPtBins] = {};
    double sigma_ratio[kMaxPtBins] = {};
    double sigma_ratio_err[kMaxPtBins] = {};
    double sigma_ratio_err_corr[kMaxPtBins] = {};
    double sigma_ratio_err_full[kMaxPtBins] = {};
    // Absolute quantities - variant
    double jms_var[kMaxPtBins] = {};
    double jms_var_err[kMaxPtBins] = {};
    double jmr_var[kMaxPtBins] = {};
    double jmr_var_err[kMaxPtBins] = {};
    // Absolute quantities - nominal
    double jms_nom[kMaxPtBins] = {};
    double jms_nom_err[kMaxPtBins] = {};
    double jmr_nom[kMaxPtBins] = {};
    double jmr_nom_err[kMaxPtBins] = {};

    // Raw fit components (peak mu, width sigma) for Data and MC-used by the
    // combined peak/resolution Data-vs-MC plots (the b-JMDP / b-JMDR components).
    double mu_data_var[kMaxPtBins] = {}, mu_data_var_err[kMaxPtBins] = {};
    double mu_mc_var[kMaxPtBins] = {}, mu_mc_var_err[kMaxPtBins] = {};
    double mu_data_nom[kMaxPtBins] = {}, mu_data_nom_err[kMaxPtBins] = {};
    double mu_mc_nom[kMaxPtBins] = {}, mu_mc_nom_err[kMaxPtBins] = {};
    double sig_data_var[kMaxPtBins] = {}, sig_data_var_err[kMaxPtBins] = {};
    double sig_mc_var[kMaxPtBins] = {}, sig_mc_var_err[kMaxPtBins] = {};
    double sig_data_nom[kMaxPtBins] = {}, sig_data_nom_err[kMaxPtBins] = {};
    double sig_mc_nom[kMaxPtBins] = {}, sig_mc_nom_err[kMaxPtBins] = {};
};

// ============================================================================
// b-JMDR "width ratio" quantity, per-variant  (WIDTH-RATIO VERSION)
// ============================================================================
// b-JMDR within one variant: ratio of the data and MC peak widths,
//
//   JMDR = sigma_data / sigma_MC
//
// The data and MC widths come from separate fits to separate histograms, so
// their errors are independent and add in quadrature:
//
//   (dJMDR / JMDR)^2 = (dsigma_data / sigma_data)^2 + (dsigma_MC / sigma_MC)^2
//
// ============================================================================
inline double CalculateResolution(double data_sigma, double mc_sigma) {
    if (mc_sigma == 0) return 0.0;
    return data_sigma / mc_sigma;
}

inline double CalculateResolutionError(double data_sigma, double data_sigma_err, double mc_sigma,
                                       double mc_sigma_err) {
    if (data_sigma == 0 || mc_sigma == 0) return 0.0;
    double ratio = data_sigma / mc_sigma;
    double rel_d = data_sigma_err / data_sigma;
    double rel_m = mc_sigma_err / mc_sigma;
    return std::abs(ratio) * std::sqrt(rel_d * rel_d + rel_m * rel_m);
}


/* ----------------------------------------------------------------------------
 * b-JMDP double ratio errors
 *
 * DR = (mu_d_var * mu_m_nom) / (mu_m_var * mu_d_nom), and d(ln DR) has the
 * coefficients c = (1/mu_d_var, -1/mu_d_nom, -1/mu_m_var, 1/mu_m_nom).
 *
 * Corr-aware: variant and nominal are fitted on the same events, so one
 * common fluctuation moves mu_d_var and mu_d_nom together (rho = 1), with an
 * error taken as the geometric mean of the two fit errors:
 *   (dDR/DR)^2 = (1/mu_d_var - 1/mu_d_nom)^2 * dmu_d_var * dmu_d_nom + (MC)
 * The contribution vanishes when the variant does not move the peak.
 *
 * Full-corr: c^T Sigma c on (mu_d_var, mu_d_nom, mu_m_var, mu_m_nom), with
 * the rho's from correlations.txt. At rho_dd = rho_mm = 1 and no cross terms
 * it gives (dmu_d_var/mu_d_var - dmu_d_nom/mu_d_nom)^2 + (MC), which equals
 * the corr-aware result only when the variant and nominal fit errors are
 * the same.
 * -------------------------------------------------------------------------- */
inline double JMSDoubleRatioErrCorrAware(double DR, double mu_d_var, double mu_d_var_err, double mu_d_nom,
                                         double mu_d_nom_err, double mu_m_var, double mu_m_var_err,
                                         double mu_m_nom, double mu_m_nom_err) {
    if (DR == 0 || mu_d_var == 0 || mu_d_nom == 0 || mu_m_var == 0 || mu_m_nom == 0) return 0.0;

    double dmu_d_common = std::sqrt(std::abs(mu_d_var_err * mu_d_nom_err));
    double dmu_m_common = std::sqrt(std::abs(mu_m_var_err * mu_m_nom_err));

    double bracket_d = (1.0 / mu_d_var) - (1.0 / mu_d_nom);
    double bracket_m = (1.0 / mu_m_nom) - (1.0 / mu_m_var);

    double rel2 = bracket_d * bracket_d * dmu_d_common * dmu_d_common
                  + bracket_m * bracket_m * dmu_m_common * dmu_m_common;

    return std::abs(DR) * std::sqrt(rel2);
}


// b-JMDR corr-aware double ratio error: same formula as the b-JMDP one above,
// with the widths sigma in place of the peaks mu,
//
//   DR = (sigma_d_var * sigma_m_nom) / (sigma_m_var * sigma_d_nom)
//
// The rho = 1 assumption is weaker here than for the peak: the bJER
// calibrations are designed to change the width, so the variant and nominal
// width fits (and their errors) differ more than the peak fits do.
inline double JMRDoubleRatioErrCorrAware(double DR, double sigma_d_var, double sigma_d_var_err,
                                         double sigma_d_nom, double sigma_d_nom_err, double sigma_m_var,
                                         double sigma_m_var_err, double sigma_m_nom, double sigma_m_nom_err) {
    if (DR == 0 || sigma_d_var == 0 || sigma_d_nom == 0 || sigma_m_var == 0 || sigma_m_nom == 0) return 0.0;

    double dsig_d_common = std::sqrt(std::abs(sigma_d_var_err * sigma_d_nom_err));
    double dsig_m_common = std::sqrt(std::abs(sigma_m_var_err * sigma_m_nom_err));

    double bracket_d = (1.0 / sigma_d_var) - (1.0 / sigma_d_nom);
    double bracket_m = (1.0 / sigma_m_nom) - (1.0 / sigma_m_var);

    double rel2 = bracket_d * bracket_d * dsig_d_common * dsig_d_common
                  + bracket_m * bracket_m * dsig_m_common * dsig_m_common;

    return std::abs(DR) * std::sqrt(rel2);
}


// Full-correlation errors: the general version of the corr-aware errors
// above, with the rho's read from a file instead of fixed to 1 (format,
// file location and defaults in CorrelationConfig.h). rho = 1 gives the
// smallest error, so the corr-aware result is a lower limit; rho < 1 when
// the calibration change adds noise of its own or moves events in or out
// of the selection.
using CorrelationConfig::CorrMethod;
using CorrelationConfig::CurrentCorrMethod;
using CorrelationConfig::CorrMethodTag;
using CorrelationConfig::CorrelationMatrix;
using CorrelationConfig::CorrKey;
using CorrelationConfig::Correlations;
using CorrelationConfig::CorrFor;
using CorrelationConfig::CorrDomain;

// ============================================================================
// Full-corr b-JMDP error. The four fitted peaks
//
//   x = (mu_d_var, mu_d_nom, mu_m_var, mu_m_nom)
//
// have covariance Sigma_ij = rho_ij * dx_i * dx_j, and the log-derivatives of
//
//   DR = (mu_d_var * mu_m_nom) / (mu_m_var * mu_d_nom)
//
// are c = (1/mu_d_var, -1/mu_d_nom, -1/mu_m_var, 1/mu_m_nom), so
//   (dDR/DR)^2 = c^T Sigma c
//
// rho for each pair (keys of correlations.txt, see CorrelationConfig.h):
//   rho(1,2) = rho_mu_dd      data variant - data nominal   (same events, ~1)
//   rho(3,4) = rho_mu_mm      MC variant   - MC nominal     (same events, ~1)
//   rho(1,3) = rho_mu_dv_mv   rho(1,4) = rho_mu_dv_mn       data - MC,
//   rho(2,3) = rho_mu_dn_mv   rho(2,4) = rho_mu_dn_mn       ~0 (different samples)
// Defaults are 1 for the same-event pairs and 0 for the others. With those,
// the result is close to the corr-aware one, and equal only when the
// variant and nominal fit errors are the same.
// ============================================================================

inline double JMSDoubleRatioErrFullCorr(double DR, double mu_d_var, double mu_d_var_err, double mu_d_nom,
                                        double mu_d_nom_err, double mu_m_var, double mu_m_var_err,
                                        double mu_m_nom, double mu_m_nom_err, const CorrelationMatrix& C) {
    if (DR == 0 || mu_d_var == 0 || mu_d_nom == 0 || mu_m_var == 0 || mu_m_nom == 0) return 0.0;

    // 4-vector of log-derivatives
    double c[4] = {
        1.0 / mu_d_var,  // d(ln DR)/dmu_d^var
        -1.0 / mu_d_nom, // d(ln DR)/dmu_d^nom
        -1.0 / mu_m_var, // d(ln DR)/dmu_m^var
        1.0 / mu_m_nom   // d(ln DR)/dmu_m^nom
    };
    // 4-vector of fit errors
    double d[4] = {mu_d_var_err, mu_d_nom_err, mu_m_var_err, mu_m_nom_err};

    // 4x4 correlation matrix (diagonal = 1, symmetric)
    double rho[4][4] = {{1.0, C.rho_mu_dd, C.rho_mu_dv_mv, C.rho_mu_dv_mn},
                        {C.rho_mu_dd, 1.0, C.rho_mu_dn_mv, C.rho_mu_dn_mn},
                        {C.rho_mu_dv_mv, C.rho_mu_dn_mv, 1.0, C.rho_mu_mm},
                        {C.rho_mu_dv_mn, C.rho_mu_dn_mn, C.rho_mu_mm, 1.0}};

    // c^T SUM(c) = SUM_ij c_i c_j Rho_ij d_i d_j
    double rel2 = 0.0;
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            rel2 += c[i] * c[j] * rho[i][j] * d[i] * d[j];
        }
    }
    rel2 = std::max(0.0, rel2);
    return std::abs(DR) * std::sqrt(rel2);
}

// Full-corr b-JMDR error: same as the b-JMDP one above, with the widths sigma
// in place of the peaks mu and the rho_sig_* keys of correlations.txt
// (rho_sig_dd, rho_sig_mm, rho_sig_dv_mv, ...) in the same matrix positions.
inline double JMRDoubleRatioErrFullCorr(double DR, double sigma_d_var, double sigma_d_var_err,
                                        double sigma_d_nom, double sigma_d_nom_err, double sigma_m_var,
                                        double sigma_m_var_err, double sigma_m_nom, double sigma_m_nom_err,
                                        const CorrelationMatrix& C) {
    if (DR == 0 || sigma_d_var == 0 || sigma_d_nom == 0 || sigma_m_var == 0 || sigma_m_nom == 0) return 0.0;

    // 4-vector of log-derivatives
    double c[4] = {
        1.0 / sigma_d_var,  // d(ln DR)/d(sigma_d_var)
        -1.0 / sigma_d_nom, // d(ln DR)/d(sigma_d_nom)
        -1.0 / sigma_m_var, // d(ln DR)/d(sigma_m_var)
        1.0 / sigma_m_nom   // d(ln DR)/d(sigma_m_nom)
    };
    // 4-vector of fit errors
    double d[4] = {sigma_d_var_err, sigma_d_nom_err, sigma_m_var_err, sigma_m_nom_err};

    // 4x4 correlation matrix
    double rho[4][4] = {{1.0, C.rho_sig_dd, C.rho_sig_dv_mv, C.rho_sig_dv_mn},
                        {C.rho_sig_dd, 1.0, C.rho_sig_dn_mv, C.rho_sig_dn_mn},
                        {C.rho_sig_dv_mv, C.rho_sig_dn_mv, 1.0, C.rho_sig_mm},
                        {C.rho_sig_dv_mn, C.rho_sig_dn_mn, C.rho_sig_mm, 1.0}};

    double rel2 = 0.0;
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            rel2 += c[i] * c[j] * rho[i][j] * d[i] * d[j];
        }
    }
    rel2 = std::max(0.0, rel2);
    return std::abs(DR) * std::sqrt(rel2);
}


/* ----------------------------------------------------------------------------
 * fit_summary.txt parser
 *
 * Reads the "Gaussian Fits" section of the file written by Mass_Fit.h, up to
 * "Mass Shift". The mu value is on the same line as the "MC:" / "Data:"
 * marker, the sigma value is on the next line.
 *
 * The section names and markers must match what Mass_Fit.h writes (the same
 * format is read by Mass_Ratios_signedJ.h). Any field that cannot be read
 * stays 0 and a warning is printed; zeros then show up as empty ratio plots.
 * -------------------------------------------------------------------------- */
inline void ReadFitResultsEnhanced(const std::string& summaryPath, MassFitResults& data_results,
                                   MassFitResults& mc_results) {
    std::ifstream infile(summaryPath);
    if (!infile.is_open()) {
        std::cerr << "[ReadFit] WARN: Could not open " << summaryPath << std::endl;
        return;
    }

    static const std::regex re_pm(
        R"(([-+]?\d+\.?\d*(?:[eE][-+]?\d+)?)\s*(?:\xC2\xB1|\xb1|\+/-|\+\\-)\s*([-+]?\d+\.?\d*(?:[eE][-+]?\d+)?))");

    auto parse_pair = [](const std::string& line, double& val, double& err) -> bool {
        std::smatch m;
        if (std::regex_search(line, m, re_pm)) {
            val = std::stod(m[1].str());
            err = std::stod(m[2].str());
            return true;
        }
        return false;
    };

    auto looks_like_mu_line = [](const std::string& line) {
        return line.find("\xCE\xBC") != std::string::npos // UTF-8 mu
               || line.find("mu ") != std::string::npos   // ASCII
               || line.find("Mu ") != std::string::npos || line.find("MU ") != std::string::npos;
    };
    auto looks_like_sigma_line = [](const std::string& line) {
        return line.find("\xCF\x83") != std::string::npos // UTF-8 sigma
               || line.find("sigma ") != std::string::npos || line.find("Sigma ") != std::string::npos
               || line.find("SIGMA ") != std::string::npos;
    };

    std::string line;
    bool in_gaussian_section = false;
    bool reading_mc = false;
    bool reading_data = false;
    int parsed_mc_mu = 0, parsed_mc_sg = 0, parsed_dt_mu = 0, parsed_dt_sg = 0;

    while (std::getline(infile, line)) {
        if (line.find("Gaussian Fits") != std::string::npos) {
            in_gaussian_section = true;
            continue;
        }
        if (!in_gaussian_section) continue;
        if (line.find("Mass Shift") != std::string::npos) break;

        // Section markers: "MC:" and "Data:". The mu line for that section
        // is on the same line as the "MC:" / "Data:" marker (writer convention).
        bool has_MC = (line.find("MC:") != std::string::npos);
        bool has_Data = (line.find("Data:") != std::string::npos);

        if (has_MC) {
            reading_mc = true;
            reading_data = false;
            if (looks_like_mu_line(line) && parse_pair(line, mc_results.mean, mc_results.mean_err))
                ++parsed_mc_mu;
            continue;
        }
        if (has_Data) {
            reading_data = true;
            reading_mc = false;
            if (looks_like_mu_line(line) && parse_pair(line, data_results.mean, data_results.mean_err))
                ++parsed_dt_mu;
            continue;
        }

        // Sigma lines are on their OWN line (writer convention) immediately following
        // the MC: or Data: line, identified by section context.
        if (looks_like_sigma_line(line)) {
            if (reading_mc && parse_pair(line, mc_results.sigma, mc_results.sigma_err)) {
                ++parsed_mc_sg;
                reading_mc = false;
            } else if (reading_data && parse_pair(line, data_results.sigma, data_results.sigma_err)) {
                ++parsed_dt_sg;
                reading_data = false;
            }
        }
    }
    infile.close();

    // Diagnostic - one line per file
    std::cout << "[ReadFit] " << summaryPath << "\n"
              << "          MC:   mu=" << mc_results.mean << "+/-" << mc_results.mean_err
              << "  sg=" << mc_results.sigma << "+/-" << mc_results.sigma_err
              << "  (parsed mu=" << parsed_mc_mu << " sg=" << parsed_mc_sg << ")\n"
              << "          Data: mu=" << data_results.mean << "+/-" << data_results.mean_err
              << "  sg=" << data_results.sigma << "+/-" << data_results.sigma_err
              << "  (parsed mu=" << parsed_dt_mu << " sg=" << parsed_dt_sg << ")\n";

    if (parsed_mc_mu == 0 || parsed_mc_sg == 0 || parsed_dt_mu == 0 || parsed_dt_sg == 0) {
        std::cerr << "[ReadFit] WARN: Some fields could not be parsed in " << summaryPath << "\n"
                  << "          Check fit_summary.txt format. Expected lines like:\n"
                  << "            Gaussian Fits:\n"
                  << "              MC: mu = 90.123 +/- 0.045 GeV\n"
                  << "                  sigma = 5.678 +/- 0.012 GeV\n"
                  << "              Data: mu = 91.234 +/- 0.067 GeV\n"
                  << "                    sigma = 6.789 +/- 0.023 GeV\n";
    }
}


// ============================================================================
// CalculateRatios - Double ratios of one variant against the nominal
// ============================================================================
inline RatioResults CalculateRatios(const std::string& dataset_before, const std::string& dataset_after,
                                    const std::string& wp, const std::string& label, int color,
                                    int marker_style, const PtBinConfig::PtBinDef& ptCfg) {
    RatioResults results;
    results.label = label;
    results.color = color;
    results.marker_style = marker_style;

    DatasetConfig config_before = DatasetConfig::FromVariant(GetVariant(dataset_before));
    DatasetConfig config_after = DatasetConfig::FromVariant(GetVariant(dataset_after));
    std::string fileWP = PathBuilder::WPtoFileWP(wp);

    MassFitResults data_before[kMaxPtBins], data_after[kMaxPtBins];
    MassFitResults mc_before[kMaxPtBins], mc_after[kMaxPtBins];

    // Read fit summaries. The pT-bin directory is suffixed with the current
    // binning mode (see common/BinningConfig.h) so this reads from the same
    // folder that Mass_Analysis wrote to under the same EASYBJETS_BINNING.
    // Default is "fixed"; set EASYBJETS_BINNING=fd (or scott/sturges/peakaware...)
    const std::string binTag = BinningConfig::ModeTag();
    for (int iPt = 0; iPt < ptCfg.nBins; ++iPt) {
        std::string base_b = (RunConfig::PlotDir() + "/") + config_before.output_subdir + "/Mass_Fit_"
                             + CurrentChannel() + "/" + fileWP + "/" + ptCfg.fileNames[iPt] + "_" + binTag
                             + "/";
        std::string base_a = (RunConfig::PlotDir() + "/") + config_after.output_subdir + "/Mass_Fit_"
                             + CurrentChannel() + "/" + fileWP + "/" + ptCfg.fileNames[iPt] + "_" + binTag
                             + "/";
        ReadFitResultsEnhanced(base_b + "fit_summary.txt", data_before[iPt], mc_before[iPt]);
        ReadFitResultsEnhanced(base_a + "fit_summary.txt", data_after[iPt], mc_after[iPt]);
    }

    // Per-bin DR computation
    for (int b = 0; b < ptCfg.nBins; ++b) {

        // ---------- JMS double ratio ----------
        double mu_d_var = data_after[b].mean, emu_d_var = data_after[b].mean_err;
        double mu_d_nom = data_before[b].mean, emu_d_nom = data_before[b].mean_err;
        double mu_m_var = mc_after[b].mean, emu_m_var = mc_after[b].mean_err;
        double mu_m_nom = mc_before[b].mean, emu_m_nom = mc_before[b].mean_err;

        double jms_var = 0, jms_nom = 0;
        if (mu_m_var != 0) jms_var = mu_d_var / mu_m_var;
        if (mu_m_nom != 0) jms_nom = mu_d_nom / mu_m_nom;

        if (jms_nom != 0 && jms_var != 0) {
            double DR = jms_var / jms_nom;
            results.mean_ratio[b] = DR;

            // Uncorrelated (naive)
            double rel_d_v = (mu_d_var != 0) ? emu_d_var / mu_d_var : 0;
            double rel_d_n = (mu_d_nom != 0) ? emu_d_nom / mu_d_nom : 0;
            double rel_m_v = (mu_m_var != 0) ? emu_m_var / mu_m_var : 0;
            double rel_m_n = (mu_m_nom != 0) ? emu_m_nom / mu_m_nom : 0;
            results.mean_ratio_err[b] =
                std::abs(DR)
                * std::sqrt(rel_d_v * rel_d_v + rel_d_n * rel_d_n + rel_m_v * rel_m_v + rel_m_n * rel_m_n);

            // Correlation-aware (Rho = 1 assumed everywhere)
            results.mean_ratio_err_corr[b] = JMSDoubleRatioErrCorrAware(
                DR, mu_d_var, emu_d_var, mu_d_nom, emu_d_nom, mu_m_var, emu_m_var, mu_m_nom, emu_m_nom);

            // Full correlation-matrix (Rho from correlations.txt)
            // rho measured by Bootstrap_Correlations for this correction and pT bin
            // ("<correction>:<ptbin>"); default matrix if not in the file.
            CorrelationMatrix C = CorrelationConfig::CorrFor(CorrelationConfig::CorrDomain::Mass, CurrentChannel(),
                                                             fileWP, ptCfg.fileNames[b],
                                                             GetVariant(dataset_after).path_tag);
            results.mean_ratio_err_full[b] = JMSDoubleRatioErrFullCorr(
                DR, mu_d_var, emu_d_var, mu_d_nom, emu_d_nom, mu_m_var, emu_m_var, mu_m_nom, emu_m_nom, C);
        }

        // ---------- Absolute b-JMDP per variant  ----------
        // Save the in-variant b-JMDP = mu_d/mu_m for variant and nominal.
        // Errors are naive uncorrelated propagation (mu_d and mu_m come from
        // independent fits - correct WITHIN one variant).
        // Saved OUTSIDE the DR guard above so we capture the absolute values
        {
            double rel_d_v = (mu_d_var != 0) ? emu_d_var / mu_d_var : 0;
            double rel_m_v = (mu_m_var != 0) ? emu_m_var / mu_m_var : 0;
            results.jms_var[b] = jms_var;
            results.jms_var_err[b] = std::abs(jms_var) * std::sqrt(rel_d_v * rel_d_v + rel_m_v * rel_m_v);

            double rel_d_n = (mu_d_nom != 0) ? emu_d_nom / mu_d_nom : 0;
            double rel_m_n = (mu_m_nom != 0) ? emu_m_nom / mu_m_nom : 0;
            results.jms_nom[b] = jms_nom;
            results.jms_nom_err[b] = std::abs(jms_nom) * std::sqrt(rel_d_n * rel_d_n + rel_m_n * rel_m_n);
        }

        // ---------- JMR double ratio ----------
        double sd_v = data_after[b].sigma, esd_v = data_after[b].sigma_err;
        double sd_n = data_before[b].sigma, esd_n = data_before[b].sigma_err;
        double sm_v = mc_after[b].sigma, esm_v = mc_after[b].sigma_err;
        double sm_n = mc_before[b].sigma, esm_n = mc_before[b].sigma_err;

        // Capture the raw peak (mu) and width (sigma) components for the CSV.
        results.mu_data_var[b] = mu_d_var;
        results.mu_data_var_err[b] = emu_d_var;
        results.mu_mc_var[b] = mu_m_var;
        results.mu_mc_var_err[b] = emu_m_var;
        results.mu_data_nom[b] = mu_d_nom;
        results.mu_data_nom_err[b] = emu_d_nom;
        results.mu_mc_nom[b] = mu_m_nom;
        results.mu_mc_nom_err[b] = emu_m_nom;
        results.sig_data_var[b] = sd_v;
        results.sig_data_var_err[b] = esd_v;
        results.sig_mc_var[b] = sm_v;
        results.sig_mc_var_err[b] = esm_v;
        results.sig_data_nom[b] = sd_n;
        results.sig_data_nom_err[b] = esd_n;
        results.sig_mc_nom[b] = sm_n;
        results.sig_mc_nom_err[b] = esm_n;

        double jmr_var = CalculateResolution(sd_v, sm_v);
        double jmr_nom = CalculateResolution(sd_n, sm_n);
        double ejmr_var = CalculateResolutionError(sd_v, esd_v, sm_v, esm_v);
        double ejmr_nom = CalculateResolutionError(sd_n, esd_n, sm_n, esm_n);

        // ---------- Absolute b-JMDR per variant  ----------
        // Save JMR = sigma_d / sigma_m for variant and nominal separately.
        // Errors reuse CalculateResolutionError (in-variant uncorrelated
        // propagation of sigma_d and sigma_m).
        // Saved outside the DR guard below - we want absolute points even
        // in cells where the DR is degenerate.
        results.jmr_var[b] = jmr_var;
        results.jmr_var_err[b] = ejmr_var;
        results.jmr_nom[b] = jmr_nom;
        results.jmr_nom_err[b] = ejmr_nom;

        // JMR is now the simple width ratio sigma_d / sigma_m (dimensionless,
        // always positive when sigmas are positive). The DR is well-defined
        // whenever both jmr_var and jmr_nom are nonzero.
        if (jmr_nom > 1e-10 && jmr_var > 1e-10) {
            double DR = jmr_var / jmr_nom;
            results.sigma_ratio[b] = DR;

            // Uncorrelated (naive)
            double rel_v = ejmr_var / jmr_var;
            double rel_n = ejmr_nom / jmr_nom;
            results.sigma_ratio_err[b] = std::abs(DR) * std::sqrt(rel_v * rel_v + rel_n * rel_n);

            // Correlation-aware (rho = 1 assumed everywhere)
            results.sigma_ratio_err_corr[b] =
                JMRDoubleRatioErrCorrAware(DR, sd_v, esd_v, sd_n, esd_n, sm_v, esm_v, sm_n, esm_n);

            // Full correlation-matrix (rho from correlations.txt)
            // rho measured by Bootstrap_Correlations for this correction and pT bin
            // ("<correction>:<ptbin>"); default matrix if not in the file.
            CorrelationMatrix C = CorrelationConfig::CorrFor(CorrelationConfig::CorrDomain::Mass, CurrentChannel(),
                                                             fileWP, ptCfg.fileNames[b],
                                                             GetVariant(dataset_after).path_tag);
            results.sigma_ratio_err_full[b] =
                JMRDoubleRatioErrFullCorr(DR, sd_v, esd_v, sd_n, esd_n, sm_v, esm_v, sm_n, esm_n, C);
        } else if (jmr_nom > 1e-10) {
            results.sigma_ratio[b] = 0;
            results.sigma_ratio_err[b] = ejmr_var / jmr_nom;
            results.sigma_ratio_err_corr[b] = results.sigma_ratio_err[b];
            results.sigma_ratio_err_full[b] = results.sigma_ratio_err[b];
        }
    }

    return results;
}


// ============================================================================
// Generic vs-pT summary plot - one TGraphErrors per variant
// ============================================================================
inline void DrawCombinedRatioPlot_vs_pT(const std::vector<RatioResults>& all_results,
                                        const PtBinConfig::PtBinDef& ptCfg,
                                        std::function<double(const RatioResults&, int)> getValue,
                                        std::function<double(const RatioResults&, int)> getError,
                                        const std::string& canvasName, const std::string& yTitle,
                                        const std::string& outputPath, const std::string& wp,
                                        double y_lo_force = std::numeric_limits<double>::quiet_NaN(),
                                        double y_hi_force = std::numeric_limits<double>::quiet_NaN()) {
    TCanvas* c = new TCanvas(canvasName.c_str(), canvasName.c_str(), 800, 600);
    c->SetLeftMargin(0.15);
    c->SetBottomMargin(0.15);
    c->SetTopMargin(0.10);
    c->SetRightMargin(0.05);

    // y range
    double ymin = 1e6, ymax = -1e6;
    int n_valid = 0;
    for (const auto& r : all_results) {
        for (int b = 0; b < ptCfg.nBins; ++b) {
            double v = getValue(r, b);
            double e = std::abs(getError(r, b));
            if (!std::isfinite(v) || !std::isfinite(e)) continue;
            ymin = std::min(ymin, v - e);
            ymax = std::max(ymax, v + e);
            ++n_valid;
        }
    }
    if (n_valid == 0) {
        std::cerr << "[Mass_Ratios] WARN: " << canvasName
                  << " has zero valid points across all variants and pT bins.\n"
                  << "             Check the [ReadFit] diagnostic above to see which "
                  << "fit_summary.txt files failed to parse.\n";
    }
    if (ymin > ymax) {
        ymin = 0.95;
        ymax = 1.05;
    }
    double range = ymax - ymin;
    double y_lo = ymin - 0.30 * range;
    double y_hi = ymax + 0.30 * range;
    if (y_lo > 0.95) y_lo = 0.95;
    if (y_hi < 1.05) y_hi = 1.05;

    // Override (JMS/JMR-specific forced ranges)
    if (!std::isnan(y_lo_force)) y_lo = y_lo_force;
    if (!std::isnan(y_hi_force)) y_hi = y_hi_force;

    TH1F* frame = new TH1F((canvasName + "_frame").c_str(), "", ptCfg.nBins, ptCfg.display_bins);
    frame->SetMinimum(y_lo);
    frame->SetMaximum(y_hi);
    frame->GetXaxis()->SetTitle("p_{T} [GeV]");
    frame->GetYaxis()->SetTitle(yTitle.c_str());
    frame->GetXaxis()->SetLabelSize(0.035);
    frame->GetYaxis()->SetLabelSize(0.035);
    frame->GetXaxis()->SetTitleSize(0.04);
    frame->GetYaxis()->SetTitleSize(0.035);
    frame->GetYaxis()->SetTitleOffset(1.7);
    frame->SetStats(0);
    frame->Draw("AXIS");

    TLine* ref = new TLine(ptCfg.display_bins[0], 1.0, ptCfg.display_bins[ptCfg.nBins], 1.0);
    ref->SetLineStyle(2);
    ref->SetLineColor(kBlack);
    ref->SetLineWidth(2);
    ref->Draw("same");

    gStyle->SetLegendBorderSize(1);
    TLegend* leg = new TLegend(0.65, 0.7, 0.9, 0.85);
    leg->SetBorderSize(1);
    leg->SetLineColor(kBlack);
    leg->SetLineWidth(1);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.025);

    double ptCenters[kMaxPtBins];
    double ptHalfWidths[kMaxPtBins];
    for (int i = 0; i < ptCfg.nBins; ++i) {
        ptCenters[i] = (ptCfg.display_bins[i] + ptCfg.display_bins[i + 1]) / 2.0;
        ptHalfWidths[i] = (ptCfg.display_bins[i + 1] - ptCfg.display_bins[i]) / 2.0;
    }

    // Per-variant horizontal offset for marker visibility.
    const int N = (int)all_results.size();
    auto offsetFor = [&](int idx, int ipt) -> double {
        if (N <= 1) return 0.0;
        double frac = (idx - (N - 1) / 2.0) / (double)N;
        return frac * 0.10 * ptHalfWidths[ipt] * 2.0;
    };

    int idx = 0;
    for (const auto& r : all_results) {
        double xs[kMaxPtBins], ys[kMaxPtBins], exs[kMaxPtBins], eys[kMaxPtBins];
        for (int i = 0; i < ptCfg.nBins; ++i) {
            xs[i] = ptCenters[i] + offsetFor(idx, i);
            ys[i] = getValue(r, i);
            exs[i] = 0;
            eys[i] = std::abs(getError(r, i));
        }
        TGraphErrors* g = new TGraphErrors(ptCfg.nBins, xs, ys, exs, eys);
        g->SetMarkerStyle(r.marker_style);
        g->SetMarkerSize(1.3);
        g->SetMarkerColor(r.color);
        g->SetLineColor(r.color);
        g->SetLineWidth(2);
        g->Draw("P same");
        leg->AddEntry(g, r.label.c_str(), "lp");
        ++idx;
    }

    leg->Draw();
    AtlasUtils::DrawAtlasLabel(0.18, 0.85);

    TLatex wp_label;
    wp_label.SetNDC();
    wp_label.SetTextSize(0.035);
    // Topology under the ATLAS Internal / lumi block (left side).
    wp_label.DrawLatex(0.18, 0.75, TopologyLabelFor(CurrentChannel()).c_str());

    c->SaveAs(outputPath.c_str());

    delete frame;
    delete c;
    delete ref;
    delete leg;
}


// ============================================================================
// Side-by-side error-treatment comparison (uncorr vs corr-aware)
// ============================================================================
// Same central values with the uncorrelated (open, dashed) and the
// correlation-aware (filled) error bars side by side.
// ============================================================================
inline void DrawErrorComparisonPlot_vs_pT(const std::vector<RatioResults>& all_results,
                                          const PtBinConfig::PtBinDef& ptCfg,
                                          std::function<double(const RatioResults&, int)> getValue,
                                          std::function<double(const RatioResults&, int)> getErrUncorr,
                                          std::function<double(const RatioResults&, int)> getErrCorr,
                                          const std::string& canvasName, const std::string& yTitle,
                                          const std::string& outputPath, const std::string& wp) {
    TCanvas* c = new TCanvas(canvasName.c_str(), canvasName.c_str(), 950, 700);
    c->SetLeftMargin(0.15);
    c->SetBottomMargin(0.15);
    c->SetTopMargin(0.10);
    c->SetRightMargin(0.05);

    // y range from the (larger) uncorrelated bars
    double ymin = 1e6, ymax = -1e6;
    for (const auto& r : all_results) {
        for (int b = 0; b < ptCfg.nBins; ++b) {
            double v = getValue(r, b);
            double e = std::abs(getErrUncorr(r, b));
            if (!std::isfinite(v) || !std::isfinite(e)) continue;
            ymin = std::min(ymin, v - e);
            ymax = std::max(ymax, v + e);
        }
    }
    if (ymin > ymax) {
        ymin = 0.95;
        ymax = 1.05;
    }
    double range = ymax - ymin;
    double y_lo = ymin - 0.30 * range;
    double y_hi = ymax + 0.30 * range;
    if (y_lo > 0.95) y_lo = 0.95;
    if (y_hi < 1.05) y_hi = 1.05;

    TH1F* frame = new TH1F((canvasName + "_frame").c_str(), "", ptCfg.nBins, ptCfg.display_bins);
    frame->SetMinimum(y_lo);
    frame->SetMaximum(y_hi);
    frame->GetXaxis()->SetTitle("p_{T} [GeV]");
    frame->GetYaxis()->SetTitle(yTitle.c_str());
    frame->GetXaxis()->SetLabelSize(0.035);
    frame->GetYaxis()->SetLabelSize(0.035);
    frame->GetXaxis()->SetTitleSize(0.04);
    frame->GetYaxis()->SetTitleSize(0.035);
    frame->GetYaxis()->SetTitleOffset(1.7);
    frame->SetStats(0);
    frame->Draw("AXIS");

    TLine* ref = new TLine(ptCfg.display_bins[0], 1.0, ptCfg.display_bins[ptCfg.nBins], 1.0);
    ref->SetLineStyle(2);
    ref->SetLineColor(kBlack);
    ref->SetLineWidth(2);
    ref->Draw("same");

    TLegend* leg = new TLegend(0.55, 0.55, 0.93, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.025);

    double ptCenters[kMaxPtBins];
    double ptHalfWidths[kMaxPtBins];
    for (int i = 0; i < ptCfg.nBins; ++i) {
        ptCenters[i] = (ptCfg.display_bins[i] + ptCfg.display_bins[i + 1]) / 2.0;
        ptHalfWidths[i] = (ptCfg.display_bins[i + 1] - ptCfg.display_bins[i]) / 2.0;
    }

    // offset the variants so they don't overlap
    const int N = (int)all_results.size();
    auto baseOffset = [&](int idx, int ipt) -> double {
        if (N <= 1) return 0.0;
        double frac = (idx - (N - 1) / 2.0) / (double)N;
        return frac * 0.15 * ptHalfWidths[ipt] * 2.0;
    };
    auto intraPairOffset = [&](int ipt) {
        // +/-4% of half-width spaces uncorr/corr around the variant's base
        return 0.04 * ptHalfWidths[ipt];
    };

    int idx = 0;
    for (const auto& r : all_results) {
        double xs_u[kMaxPtBins], xs_c[kMaxPtBins], ys[kMaxPtBins];
        double exs[kMaxPtBins], ey_u[kMaxPtBins], ey_c[kMaxPtBins];
        for (int i = 0; i < ptCfg.nBins; ++i) {
            double base_x = ptCenters[i] + baseOffset(idx, i);
            double dx = intraPairOffset(i);
            xs_u[i] = base_x - dx;
            xs_c[i] = base_x + dx;
            ys[i] = getValue(r, i);
            exs[i] = 0;
            ey_u[i] = std::abs(getErrUncorr(r, i));
            ey_c[i] = std::abs(getErrCorr(r, i));
        }
        TGraphErrors* g_u = new TGraphErrors(ptCfg.nBins, xs_u, ys, exs, ey_u);
        g_u->SetMarkerStyle(r.marker_style + 4);
        g_u->SetMarkerSize(1.0);
        g_u->SetMarkerColor(r.color);
        g_u->SetLineColor(r.color);
        g_u->SetLineWidth(1);
        g_u->SetLineStyle(2);
        g_u->Draw("P same");

        TGraphErrors* g_c = new TGraphErrors(ptCfg.nBins, xs_c, ys, exs, ey_c);
        g_c->SetMarkerStyle(r.marker_style);
        g_c->SetMarkerSize(1.3);
        g_c->SetMarkerColor(r.color);
        g_c->SetLineColor(r.color);
        g_c->SetLineWidth(2);
        g_c->Draw("P same");

        leg->AddEntry(g_u, (r.label + " - uncorr.").c_str(), "lp");
        leg->AddEntry(g_c, (r.label + " - corr-aware").c_str(), "lp");
        ++idx;
    }

    leg->Draw();
    AtlasUtils::DrawAtlasLabel(0.18, 0.85);

    TLatex wp_label;
    wp_label.SetNDC();
    wp_label.SetTextSize(0.035);
    wp_label.DrawLatex(0.18, 0.75, TopologyLabelFor(CurrentChannel()).c_str());

    c->SaveAs(outputPath.c_str());

    delete frame;
    delete c;
    delete ref;
    delete leg;
}


// ============================================================================
// 2-D heatmap (nVariants x nPtBins)
// ============================================================================
//
inline void Draw2DRatioMap(const std::vector<RatioResults>& all_results, const PtBinConfig::PtBinDef& ptCfg,
                           std::function<double(const RatioResults&, int)> getValue,
                           std::function<double(const RatioResults&, int)> getError,
                           const std::string& canvasName, const std::string& titleText,
                           const std::string& outputPath, const std::string& wp) {
    int nv = (int)all_results.size();
    if (nv == 0) return;

    TCanvas* c = new TCanvas(canvasName.c_str(), canvasName.c_str(), 1000, 80 + 90 * nv);
    c->SetTopMargin(0.12);
    c->SetRightMargin(0.12);
    c->SetLeftMargin(0.18);
    c->SetBottomMargin(0.15);

    TH2F* h = new TH2F((canvasName + "_h").c_str(), ";p_{T} Bin;Configuration", ptCfg.nBins, -0.5,
                       ptCfg.nBins - 0.5, nv, -0.5, nv - 0.5);

    for (int b = 0; b < ptCfg.nBins; ++b) {
        for (int v = 0; v < nv; ++v) {
            h->SetBinContent(b + 1, v + 1, getValue(all_results[v], b));
            h->SetBinError(b + 1, v + 1, std::abs(getError(all_results[v], b)));
        }
    }

    h->SetStats(0);
    h->SetMarkerSize(1.5);

    for (int b = 0; b < ptCfg.nBins; ++b) {
        TString lbl = Form("%d-%d GeV", (int)ptCfg.display_bins[b], (int)ptCfg.display_bins[b + 1]);
        h->GetXaxis()->SetBinLabel(b + 1, lbl.Data());
    }
    for (int v = 0; v < nv; ++v) h->GetYaxis()->SetBinLabel(v + 1, all_results[v].label.c_str());

    h->GetXaxis()->SetLabelSize(0.04);
    h->GetYaxis()->SetLabelSize(0.04);

    // Center palette on 1.0 with symmetric extent
    double minVal = 1e6, maxVal = -1e6;
    for (int i = 1; i <= h->GetNbinsX(); ++i) {
        for (int j = 1; j <= h->GetNbinsY(); ++j) {
            double content = h->GetBinContent(i, j);
            if (content != 0) {
                minVal = std::min(minVal, content);
                maxVal = std::max(maxVal, content);
            }
        }
    }
    double extent = std::max(std::abs(1.0 - minVal), std::abs(maxVal - 1.0));
    if (extent < 0.01) extent = 0.01;
    h->SetMinimum(1.0 - extent * 1.1);
    h->SetMaximum(1.0 + extent * 1.1);

    h->Draw("COLZ");

    // Overlay per-cell "value +/- error" labels with TLatex.
    TLatex cellTxt;
    cellTxt.SetTextAlign(22);
    cellTxt.SetTextSize(0.045);
    cellTxt.SetTextFont(42);
    cellTxt.SetTextColor(kBlack);
    for (int b = 0; b < ptCfg.nBins; ++b) {
        for (int v = 0; v < nv; ++v) {
            double val = h->GetBinContent(b + 1, v + 1);
            double err = h->GetBinError(b + 1, v + 1);
            if (!std::isfinite(val) || val == 0) continue;
            TString s = Form("%.3f #pm %.3f", val, err);
            cellTxt.DrawLatex((double)b, (double)v, s.Data());
        }
    }

    TLatex t;
    t.SetNDC();
    t.SetTextSize(0.040);
    t.DrawLatex(0.18, 0.94, Form("%s - WP: %s%%", titleText.c_str(), wp.c_str()));

    AtlasUtils::DrawAtlasLabel(0.18, 0.85);
    c->SaveAs(outputPath.c_str());

    delete h;
    delete c;
}


// ============================================================================
// Absolute per-variant plot - b-JMDP or b-JMDR per variant vs pT
// ============================================================================
inline void DrawAbsolutePlot_vs_pT(const std::vector<RatioResults>& all_results,
                                   const PtBinConfig::PtBinDef& ptCfg,
                                   std::function<double(const RatioResults&, int)> getValueVar,
                                   std::function<double(const RatioResults&, int)> getErrorVar,
                                   std::function<double(const RatioResults&, int)> getValueNom,
                                   std::function<double(const RatioResults&, int)> getErrorNom,
                                   const std::string& canvasName, const std::string& yTitle,
                                   const std::string& outputPath, const std::string& wp) {
    if (all_results.empty()) return;

    auto rawVariantLabel = [](const std::string& full) {
        const std::string suffix = "/Nominal";
        if (full.size() >= suffix.size()
            && full.compare(full.size() - suffix.size(), suffix.size(), suffix) == 0) {
            return full.substr(0, full.size() - suffix.size());
        }
        return full;
    };

    TCanvas* c = new TCanvas(canvasName.c_str(), canvasName.c_str(), 800, 600);
    c->SetLeftMargin(0.15);
    c->SetBottomMargin(0.15);
    c->SetTopMargin(0.10);
    c->SetRightMargin(0.05);

    // y range from ALL variants AND the nominal
    double ymin = 1e6, ymax = -1e6;
    int n_valid = 0;
    auto includeInRange = [&](double v, double e) {
        if (!std::isfinite(v) || !std::isfinite(e)) return;
        ymin = std::min(ymin, v - e);
        ymax = std::max(ymax, v + e);
        ++n_valid;
    };
    for (const auto& r : all_results) {
        for (int b = 0; b < ptCfg.nBins; ++b) includeInRange(getValueVar(r, b), std::abs(getErrorVar(r, b)));
    }
    for (int b = 0; b < ptCfg.nBins; ++b) {
        includeInRange(getValueNom(all_results[0], b), std::abs(getErrorNom(all_results[0], b)));
    }
    if (n_valid == 0) {
        std::cerr << "[Mass_Ratios] WARN: " << canvasName
                  << " has zero valid points (variants + nominal all empty).\n";
        delete c;
        return;
    }
    if (ymin > ymax) {
        ymin = -1;
        ymax = 1;
    }
    double range = ymax - ymin;
    double y_lo = ymin - 0.30 * range;
    double y_hi = ymax + 0.30 * range;

    TH1F* frame = new TH1F((canvasName + "_frame").c_str(), "", ptCfg.nBins, ptCfg.display_bins);
    frame->SetMinimum(y_lo);
    frame->SetMaximum(y_hi);
    frame->GetXaxis()->SetTitle("p_{T} [GeV]");
    frame->GetYaxis()->SetTitle(yTitle.c_str());
    frame->GetXaxis()->SetLabelSize(0.035);
    frame->GetYaxis()->SetLabelSize(0.035);
    frame->GetXaxis()->SetTitleSize(0.04);
    frame->GetYaxis()->SetTitleSize(0.035);
    frame->GetYaxis()->SetTitleOffset(1.7);
    frame->SetStats(0);
    frame->Draw("AXIS");

    gStyle->SetLegendBorderSize(1);
    TLegend* leg = new TLegend(0.62, 0.65, 0.92, 0.88);
    leg->SetBorderSize(1);
    leg->SetLineColor(kBlack);
    leg->SetLineWidth(1);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.025);

    double ptCenters[kMaxPtBins];
    double ptHalfWidths[kMaxPtBins];
    for (int i = 0; i < ptCfg.nBins; ++i) {
        ptCenters[i] = (ptCfg.display_bins[i] + ptCfg.display_bins[i + 1]) / 2.0;
        ptHalfWidths[i] = (ptCfg.display_bins[i + 1] - ptCfg.display_bins[i]) / 2.0;
    }

    const int N = (int)all_results.size();
    const int N_total = N + 1;
    auto offsetFor = [&](int idx, int ipt) -> double {
        if (N_total <= 1) return 0.0;
        double frac = (idx - (N_total - 1) / 2.0) / (double)N_total;
        return frac * 0.10 * ptHalfWidths[ipt] * 2.0;
    };

    // Variants
    int idx = 0;
    for (const auto& r : all_results) {
        double xs[kMaxPtBins], ys[kMaxPtBins], exs[kMaxPtBins], eys[kMaxPtBins];
        for (int i = 0; i < ptCfg.nBins; ++i) {
            xs[i] = ptCenters[i] + offsetFor(idx, i);
            ys[i] = getValueVar(r, i);
            exs[i] = 0;
            eys[i] = std::abs(getErrorVar(r, i));
        }
        TGraphErrors* g = new TGraphErrors(ptCfg.nBins, xs, ys, exs, eys);
        g->SetMarkerStyle(r.marker_style);
        g->SetMarkerSize(1.3);
        g->SetMarkerColor(r.color);
        g->SetLineColor(r.color);
        g->SetLineWidth(2);
        g->Draw("P same");
        leg->AddEntry(g, rawVariantLabel(r.label).c_str(), "lp");
        ++idx;
    }

    {
        double xs[kMaxPtBins], ys[kMaxPtBins], exs[kMaxPtBins], eys[kMaxPtBins];
        for (int i = 0; i < ptCfg.nBins; ++i) {
            xs[i] = ptCenters[i] + offsetFor(N, i);
            ys[i] = getValueNom(all_results[0], i);
            exs[i] = 0;
            eys[i] = std::abs(getErrorNom(all_results[0], i));
        }
        TGraphErrors* g_nom = new TGraphErrors(ptCfg.nBins, xs, ys, exs, eys);
        g_nom->SetMarkerStyle(21);
        g_nom->SetMarkerSize(1.4);
        g_nom->SetMarkerColor(kBlack);
        g_nom->SetLineColor(kBlack);
        g_nom->SetLineWidth(2);
        g_nom->Draw("P same");
        leg->AddEntry(g_nom, "Nominal", "lp");
    }

    leg->Draw();
    AtlasUtils::DrawAtlasLabel(0.18, 0.85);

    TLatex wp_label;
    wp_label.SetNDC();
    wp_label.SetTextSize(0.035);
    wp_label.DrawLatex(0.18, 0.75, TopologyLabelFor(CurrentChannel()).c_str());

    c->SaveAs(outputPath.c_str());

    delete frame;
    delete c;
    delete leg;
}

// Summary plot vs pT, on one canvas:
//   left axis   data/MC single ratio (filled markers, nominal in black)
//   right axis  variant/nominal double ratio (open markers)
//
// The six get* arguments return, for one variant and one pT bin, the value
// and error to draw. Passing the mu getters gives the b-JMDP plot, passing
// the sigma getters gives the b-JMDR plot, with the same code.
//
// selected: positions in all_results of the variants to draw, e.g. {k} for
// one variant or the EASYBJETS_SUMMARY_VARIANTS subset. Empty = all variants.
inline void DrawDualAxisSummary(const std::vector<RatioResults>& all_results,
                                const PtBinConfig::PtBinDef& ptCfg,
                                std::function<double(const RatioResults&, int)> getSingleVar,
                                std::function<double(const RatioResults&, int)> getSingleVarErr,
                                std::function<double(const RatioResults&, int)> getSingleNom,
                                std::function<double(const RatioResults&, int)> getSingleNomErr,
                                std::function<double(const RatioResults&, int)> getDR,
                                std::function<double(const RatioResults&, int)> getDRErr,
                                const std::string& canvasName, const std::string& leftTitle,
                                const std::string& rightTitle, const std::string& outputPath,
                                const std::string& wp, const std::vector<size_t>& selected = {}) {
    if (all_results.empty()) return;

    // Resolve which entries to draw.
    std::vector<size_t> idxs;
    if (selected.empty()) {
        for (size_t i = 0; i < all_results.size(); ++i) idxs.push_back(i);
    } else {
        idxs = selected;
    }
    if (idxs.empty()) return;

    auto rawVariantLabel = [](const std::string& full) {
        const std::string suffix = "/Nominal";
        if (full.size() >= suffix.size()
            && full.compare(full.size() - suffix.size(), suffix.size(), suffix) == 0)
            return full.substr(0, full.size() - suffix.size());
        return full;
    };

    TCanvas* c = new TCanvas(canvasName.c_str(), canvasName.c_str(), 900, 650);
    c->SetLeftMargin(0.13);
    c->SetBottomMargin(0.14);
    c->SetTopMargin(0.08);
    c->SetRightMargin(0.14); // room for the right axis

    // left axis: single ratios, with headroom for the legend
    double sMin = 1e9, sMax = -1e9;
    auto addS = [&](double v, double e) {
        if (!std::isfinite(v) || !std::isfinite(e)) return;
        sMin = std::min(sMin, v - e);
        sMax = std::max(sMax, v + e);
    };
    for (int b = 0; b < ptCfg.nBins; ++b) {
        addS(getSingleNom(all_results[idxs[0]], b), std::abs(getSingleNomErr(all_results[idxs[0]], b)));
        for (size_t k : idxs)
            addS(getSingleVar(all_results[k], b), std::abs(getSingleVarErr(all_results[k], b)));
    }
    if (sMin > sMax) {
        sMin = 0.9;
        sMax = 1.1;
    }
    double sRange = sMax - sMin;
    double sLo = sMin - 0.25 * sRange;
    double sHi = sMax + 0.40 * sRange; // headroom for the legend

    // right axis: double ratios
    double dMin = 1e9, dMax = -1e9;
    auto addD = [&](double v, double e) {
        if (!std::isfinite(v) || !std::isfinite(e)) return;
        dMin = std::min(dMin, v - e);
        dMax = std::max(dMax, v + e);
    };
    for (int b = 0; b < ptCfg.nBins; ++b)
        for (size_t k : idxs) addD(getDR(all_results[k], b), std::abs(getDRErr(all_results[k], b)));
    if (dMin > dMax) {
        dMin = 0.95;
        dMax = 1.05;
    }
    double dRange = dMax - dMin;
    if (dRange < 1e-6) dRange = 0.1;
    double dLo = dMin - 0.25 * dRange;
    double dHi = dMax + 0.25 * dRange;

    // DR values are drawn in left-axis coordinates, the right TGaxis relabels them
    auto mapDR = [&](double dr) {
        return sLo + (dr - dLo) / (dHi - dLo) * (sHi - sLo);
    };

    // ---- Frame ----
    TH1F* frame = new TH1F((canvasName + "_frame").c_str(), "", ptCfg.nBins, ptCfg.display_bins);
    frame->SetMinimum(sLo);
    frame->SetMaximum(sHi);
    frame->GetXaxis()->SetTitle("p_{T} [GeV]");
    frame->GetYaxis()->SetTitle(leftTitle.c_str());
    frame->GetXaxis()->SetLabelSize(0.035);
    frame->GetYaxis()->SetLabelSize(0.035);
    frame->GetXaxis()->SetTitleSize(0.04);
    frame->GetYaxis()->SetTitleSize(0.038);
    frame->GetYaxis()->SetTitleOffset(1.45);
    frame->SetStats(0);
    frame->Draw("AXIS");

    // Dashed reference line at single-ratio = 1
    if (sLo < 1.0 && sHi > 1.0) {
        TLine* ref = new TLine(ptCfg.display_bins[0], 1.0, ptCfg.display_bins[ptCfg.nBins], 1.0);
        ref->SetLineStyle(3);
        ref->SetLineColor(kGray + 1);
        ref->Draw("same");
    }

    // ---- Geometry: bin centers + per-variant x-offsets ----
    double ptCenters[kMaxPtBins], ptHalf[kMaxPtBins];
    for (int i = 0; i < ptCfg.nBins; ++i) {
        ptCenters[i] = (ptCfg.display_bins[i] + ptCfg.display_bins[i + 1]) / 2.0;
        ptHalf[i] = (ptCfg.display_bins[i + 1] - ptCfg.display_bins[i]) / 2.0;
    }
    const int N = (int)idxs.size();
    auto offset = [&](int slot, int ipt) {
        if (N <= 1) return -0.06 * ptHalf[ipt] * 2.0; // small left shift for the pair
        double frac = (slot - (N - 1) / 2.0) / (double)N;
        return frac * 0.10 * ptHalf[ipt] * 2.0;
    };

    // two columns: nominal on its own row, then [single | DR] per variant
    TLegend* leg = new TLegend(0.50, 0.68, 0.86, 0.90);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.024);
    leg->SetNColumns(2);
    leg->SetColumnSeparation(0.05);

    // ---- Nominal single ratio (black filled squares) ----
    {
        double xs[kMaxPtBins], ys[kMaxPtBins], exs[kMaxPtBins], eys[kMaxPtBins];
        for (int i = 0; i < ptCfg.nBins; ++i) {
            xs[i] = ptCenters[i];
            ys[i] = getSingleNom(all_results[idxs[0]], i);
            exs[i] = 0;
            eys[i] = std::abs(getSingleNomErr(all_results[idxs[0]], i));
        }
        TGraphErrors* g = new TGraphErrors(ptCfg.nBins, xs, ys, exs, eys);
        g->SetMarkerStyle(21);
        g->SetMarkerSize(1.3);
        g->SetMarkerColor(kBlack);
        g->SetLineColor(kBlack);
        g->SetLineWidth(2);
        g->Draw("P same");
        leg->AddEntry(g, "Nominal data/MC", "lp");
        leg->AddEntry((TObject*)nullptr, "", "");
    }

    // Per-variant: single ratio (filled marker, left axis) + DR (open, right axis)
    int slot = 0;
    for (size_t k : idxs) {
        const RatioResults& r = all_results[k];
        std::string lbl = rawVariantLabel(r.label);

        // single ratio
        double xs[kMaxPtBins], ys[kMaxPtBins], exs[kMaxPtBins], eys[kMaxPtBins];
        for (int i = 0; i < ptCfg.nBins; ++i) {
            xs[i] = ptCenters[i] + offset(slot, i) - 0.02 * ptHalf[i];
            ys[i] = getSingleVar(r, i);
            exs[i] = 0;
            eys[i] = std::abs(getSingleVarErr(r, i));
        }
        TGraphErrors* gS = new TGraphErrors(ptCfg.nBins, xs, ys, exs, eys);
        gS->SetMarkerStyle(r.marker_style);
        gS->SetMarkerSize(1.2);
        gS->SetMarkerColor(r.color);
        gS->SetLineColor(r.color);
        gS->SetLineWidth(2);
        gS->Draw("P same");
        leg->AddEntry(gS, (lbl + " data/MC").c_str(), "lp");

        // double ratio
        double xsd[kMaxPtBins], ysd[kMaxPtBins], eysd[kMaxPtBins];
        for (int i = 0; i < ptCfg.nBins; ++i) {
            xsd[i] = ptCenters[i] + offset(slot, i) + 0.02 * ptHalf[i];
            ysd[i] = mapDR(getDR(r, i));
            eysd[i] = std::abs(getDRErr(r, i)) / (dHi - dLo) * (sHi - sLo);
        }
        TGraphErrors* gD = new TGraphErrors(ptCfg.nBins, xsd, ysd, exs, eysd);
        gD->SetMarkerStyle(r.marker_style + 4);
        gD->SetMarkerSize(1.2);
        gD->SetMarkerColor(r.color);
        gD->SetLineColor(r.color);
        gD->SetLineWidth(1);
        gD->SetLineStyle(2);
        gD->Draw("P same");
        leg->AddEntry(gD, (lbl + " DR").c_str(), "lp");

        ++slot;
    }

    // ---- Right axis for DR ----
    double xRight = ptCfg.display_bins[ptCfg.nBins];
    TGaxis* axR = new TGaxis(xRight, sLo, xRight, sHi, dLo, dHi, 510, "+L");
    axR->SetTitle(rightTitle.c_str());
    axR->SetTitleColor(kBlack);
    axR->SetLabelColor(kBlack);
    axR->SetLabelSize(0.035);
    axR->SetTitleSize(0.038);
    axR->SetTitleOffset(1.4);
    axR->Draw();

    leg->Draw();
    AtlasUtils::DrawAtlasLabel(0.16, 0.88);

    TLatex topo;
    topo.SetNDC();
    topo.SetTextSize(0.032);
    topo.DrawLatex(0.16, 0.70, TopologyLabelFor(CurrentChannel()).c_str());

    c->SaveAs(outputPath.c_str());

    delete frame;
    delete c;
    delete leg;
    delete axR;
}

/* ----------------------------------------------------------------------------
 * Text outputs
 * -------------------------------------------------------------------------- */

// CSV for Mass_Summary_Combined.h, one row per variant and pT bin. Columns
// (0-based): 0-6 channel, WP, variant and pT bin; 7-18 single and double
// ratios of the peak (b-JMDP) and width (b-JMDR) with uncorrelated errors;
// 19-34 the data and MC mu and sigma behind them; 35-38 the corr-aware and
// full-corr double ratio errors.
// Mass_Summary_Combined reads the columns by position: add new columns only
// at the end, and update the reader at the same time.
inline void WriteMassRatiosCSV(const std::vector<RatioResults>& all_results,
                               const PtBinConfig::PtBinDef& ptCfg, const std::string& csvPath,
                               const std::string& channel, const std::string& wp_file) {
    std::ofstream out(csvPath);
    if (!out.is_open()) {
        std::cerr << "[Mass_Ratios_Width] WARN: cannot open CSV " << csvPath << "\n";
        return;
    }
    out << "channel,wp,variant,pt_bin_index,pt_lo,pt_hi,pt_center,"
        << "jmdp_nom,jmdp_nom_err,jmdp_var,jmdp_var_err,"
        << "jmdr_nom,jmdr_nom_err,jmdr_var,jmdr_var_err,"
        << "jmdp_dr,jmdp_dr_err,jmdr_dr,jmdr_dr_err,"
        << "mu_data_nom,mu_data_nom_err,mu_mc_nom,mu_mc_nom_err,"
        << "mu_data_var,mu_data_var_err,mu_mc_var,mu_mc_var_err,"
        << "sig_data_nom,sig_data_nom_err,sig_mc_nom,sig_mc_nom_err,"
        << "sig_data_var,sig_data_var_err,sig_mc_var,sig_mc_var_err,"
        << "jmdp_dr_err_corr,jmdp_dr_err_full,jmdr_dr_err_corr,jmdr_dr_err_full\n";

    out << std::setprecision(8);
    for (const auto& r : all_results) {
        // r.label is a display string like "bJER10v01/Nominal"; keep only the
        // correction suffix so the variant column reads "DB_bJER10v01".
        std::string suffix = r.label;
        size_t cut = suffix.find_first_of(" /");
        if (cut != std::string::npos) suffix = suffix.substr(0, cut);
        const std::string variantKey = channel + "_" + suffix;
        for (int b = 0; b < ptCfg.nBins; ++b) {
            double lo = ptCfg.display_bins[b];
            double hi = ptCfg.display_bins[b + 1];
            double mid = 0.5 * (lo + hi);
            out << channel << "," << wp_file << "," << variantKey << "," << b << "," << lo << "," << hi << ","
                << mid << "," << r.jms_nom[b] << "," << r.jms_nom_err[b] << "," << r.jms_var[b] << ","
                << r.jms_var_err[b] << "," << r.jmr_nom[b] << "," << r.jmr_nom_err[b] << "," << r.jmr_var[b]
                << "," << r.jmr_var_err[b] << "," << r.mean_ratio[b] << "," << r.mean_ratio_err[b] << ","
                << r.sigma_ratio[b] << "," << r.sigma_ratio_err[b] << "," << r.mu_data_nom[b] << ","
                << r.mu_data_nom_err[b] << "," << r.mu_mc_nom[b] << "," << r.mu_mc_nom_err[b] << ","
                << r.mu_data_var[b] << "," << r.mu_data_var_err[b] << "," << r.mu_mc_var[b] << ","
                << r.mu_mc_var_err[b] << "," << r.sig_data_nom[b] << "," << r.sig_data_nom_err[b] << ","
                << r.sig_mc_nom[b] << "," << r.sig_mc_nom_err[b] << "," << r.sig_data_var[b] << ","
                << r.sig_data_var_err[b] << "," << r.sig_mc_var[b] << "," << r.sig_mc_var_err[b] << ","
                << r.mean_ratio_err_corr[b] << "," << r.mean_ratio_err_full[b] << ","
                << r.sigma_ratio_err_corr[b] << "," << r.sigma_ratio_err_full[b] << "\n";
        }
    }
    out.close();
    std::cout << "[Mass_Ratios_Width] CSV saved: " << csvPath << "\n";
}

inline void WriteCombinedSummary(const std::vector<RatioResults>& all_results,
                                 const PtBinConfig::PtBinDef& ptCfg, const std::string& outputPath,
                                 const std::string& wp) {
    std::ofstream out(outputPath);
    if (!out.is_open()) {
        std::cerr << "Error: cannot write " << outputPath << std::endl;
        return;
    }
    // Remove "/Nominal" suffix from variant labels
    auto rawVariantLabel = [](const std::string& full) {
        const std::string suffix = "/Nominal";
        if (full.size() >= suffix.size()
            && full.compare(full.size() - suffix.size(), suffix.size(), suffix) == 0) {
            return full.substr(0, full.size() - suffix.size());
        }
        return full;
    };

    out << "============================================================================\n";
    out << "  Combined mass-scale and mass-resolution double ratios - WP: " << wp << "%\n";
    out << "============================================================================\n\n";
    out << "Definitions (per pT bin):\n";
    out << "  JMS = mu_data / mu_MC                   (within one variant)\n";
    out << "  JMR = sigma_data / sigma_MC             (within one variant, width ratio)\n";
    out << "  DR(X) = X_variant / X_nominal\n\n";
    out << "Error treatments:\n";
    out << "  uncorr   : naive - all 4 inputs (mu/sigma for data and MC,\n";
    out << "             variant and nominal) treated as independent.\n";
    out << "             Conservative; double-counts fluctuations that\n";
    out << "             actually move coherently between variants.\n";
    out << "  corr-aw  : data and MC fit fluctuations on the SAME events\n";
    out << "             between variants are treated as fully correlated\n";
    out << "             and cancel; only the var/nom fitted-value\n";
    out << "             differences enter the residual error.\n\n";

    out << std::fixed << std::setprecision(4);

    // ---- JMS double ratios ----
    out << "----------------------------------------------------------------------------\n";
    out << "  DR(JMS)  =  (mu_d/mu_m)_var / (mu_d/mu_m)_nom\n";
    out << "----------------------------------------------------------------------------\n";
    for (int b = 0; b < ptCfg.nBins; ++b) {
        out << "\n  pT bin: " << ptCfg.labels[b] << "\n";
        out << "  " << std::left << std::setw(28) << "Variant" << std::right << std::setw(28)
            << "DR(JMS) - uncorr" << std::setw(28) << "DR(JMS) - corr-aware" << std::setw(28)
            << "DR(JMS) - full-corr" << "\n";
        out << "  " << std::string(110, '-') << "\n";
        for (const auto& r : all_results) {
            out << "  " << std::left << std::setw(28) << r.label << std::right << std::setw(28)
                << Form("%.4f +/- %.4f", r.mean_ratio[b], r.mean_ratio_err[b]) << std::setw(28)
                << Form("%.4f +/- %.4f", r.mean_ratio[b], r.mean_ratio_err_corr[b]) << std::setw(28)
                << Form("%.4f +/- %.4f", r.mean_ratio[b], r.mean_ratio_err_full[b]) << "\n";
        }
    }
    out << "\n";

    // ---- JMR double ratios ----
    out << "----------------------------------------------------------------------------\n";
    out << "  DR(JMR)  =  (sigma_d/sigma_m)_var / (sigma_d/sigma_m)_nom\n";
    out << "----------------------------------------------------------------------------\n";
    for (int b = 0; b < ptCfg.nBins; ++b) {
        out << "\n  pT bin: " << ptCfg.labels[b] << "\n";
        out << "  " << std::left << std::setw(28) << "Variant" << std::right << std::setw(28)
            << "DR(JMR) - uncorr" << std::setw(28) << "DR(JMR) - corr-aware" << std::setw(28)
            << "DR(JMR) - full-corr" << "\n";
        out << "  " << std::string(110, '-') << "\n";
        for (const auto& r : all_results) {
            out << "  " << std::left << std::setw(28) << r.label << std::right << std::setw(28)
                << Form("%.4f +/- %.4f", r.sigma_ratio[b], r.sigma_ratio_err[b]) << std::setw(28)
                << Form("%.4f +/- %.4f", r.sigma_ratio[b], r.sigma_ratio_err_corr[b]) << std::setw(28)
                << Form("%.4f +/- %.4f", r.sigma_ratio[b], r.sigma_ratio_err_full[b]) << "\n";
        }
    }

    // ========================================================================
    // ABSOLUTE QUANTITIES
    // ========================================================================
    // Per-variant absolute b-JMDP and b-JMDR. The nominal entry is duplicated
    // across every variant's RatioResults (so any single entry carries the
    // nominal reference); we pull it from all_results[0] and print it ONCE
    // per pT bin instead of repeating it per variant. Errors here are the
    // in-variant naive uncorrelated propagation - there's no notion of
    // "corr-aware" or "full-corr" for an absolute single-variant value.
    // ========================================================================

    out << "\n\n";
    out << "============================================================================\n";
    out << "  ABSOLUTE QUANTITIES - per variant + nominal reference\n";
    out << "============================================================================\n\n";
    out << "Definitions (per pT bin):\n";
    out << "  b-JMDP = mu_data / mu_MC                                  (dimensionless)\n";
    out << "  b-JMDR = (sigma_data - sigma_MC) * sqrt|sigma_data^2 - sigma_MC^2|       \n\n";
    out << "Each variant carries its own b-JMDP / b-JMDR. The 'Nominal' row is the\n";
    out << "reference shared by all variants. Errors are naive uncorrelated propagation\n";
    out << "of the in-variant fit errors (mu_d and mu_m are from independent fits to\n";
    out << "data and MC histograms - uncorrelated within one variant).\n\n";

    // ---- Absolute b-JMDP ----
    out << "----------------------------------------------------------------------------\n";
    out << "  b-JMDP absolute  =  mu_data / mu_MC\n";
    out << "----------------------------------------------------------------------------\n";
    for (int b = 0; b < ptCfg.nBins; ++b) {
        out << "\n  pT bin: " << ptCfg.labels[b] << "\n";
        out << "  " << std::left << std::setw(28) << "Variant" << std::right << std::setw(28)
            << "b-JMDP +/- error" << "\n";
        out << "  " << std::string(56, '-') << "\n";
        // Nominal row first (from any entry - all duplicate the same nominal)
        out << "  " << std::left << std::setw(28) << "Nominal" << std::right << std::setw(28)
            << Form("%.4f +/- %.4f", all_results[0].jms_nom[b], all_results[0].jms_nom_err[b]) << "\n";
        // Then each variant
        for (const auto& r : all_results) {
            out << "  " << std::left << std::setw(28) << rawVariantLabel(r.label) << std::right
                << std::setw(28) << Form("%.4f +/- %.4f", r.jms_var[b], r.jms_var_err[b]) << "\n";
        }
    }
    out << "\n";

    // ---- Absolute b-JMDR ----
    out << "----------------------------------------------------------------------------\n";
    out << "  b-JMDR absolute  =  sigma_data / sigma_MC          (dimensionless width ratio)\n";
    out << "----------------------------------------------------------------------------\n";
    for (int b = 0; b < ptCfg.nBins; ++b) {
        out << "\n  pT bin: " << ptCfg.labels[b] << "\n";
        out << "  " << std::left << std::setw(28) << "Variant" << std::right << std::setw(28)
            << "b-JMDR +/- error" << "\n";
        out << "  " << std::string(56, '-') << "\n";
        out << "  " << std::left << std::setw(28) << "Nominal" << std::right << std::setw(28)
            << Form("%.4f +/- %.4f", all_results[0].jmr_nom[b], all_results[0].jmr_nom_err[b]) << "\n";
        for (const auto& r : all_results) {
            out << "  " << std::left << std::setw(28) << rawVariantLabel(r.label) << std::right
                << std::setw(28) << Form("%.4f +/- %.4f", r.jmr_var[b], r.jmr_var_err[b]) << "\n";
        }
    }

    out.close();
    std::cout << "Summary table saved: " << outputPath << std::endl;
}


// ============================================================================
// Main per-WP analysis
// ============================================================================
inline void PerformCombinedMassRatioWidthAnalysis(const std::string& wp, const PtBinConfig::PtBinDef& ptCfg) {
    std::cout << "\n***********************************************\n";
    std::cout << "COMBINED MASS RATIO ANALYSIS - WP: " << wp << "%\n";
    std::cout << "***********************************************\n\n";

    std::string fileWP = PathBuilder::WPtoFileWP(wp);
    // Tag the output directory with the pT mode so inclusive and ptbins
    // runs land in distinct folders and never overwrite each other.
    // Also tag with the binning mode (EASYBJETS_BINNING) so runs with
    // different binnings don't overlap.
    std::string ptTag = (ptCfg.nBins == 1) ? "inclusive" : "ptbins";
    std::string binTag = BinningConfig::ModeTag();
    std::string outputPath =
        Form((RunConfig::PlotDir() + "/Mass/Mass_Ratios_Width/Combined/WP_%s_%s_%s/").c_str(), fileWP.c_str(),
             ptTag.c_str(), binTag.c_str());
    gSystem->Exec(("mkdir -p " + outputPath).c_str());

    // every variant of the channel against the nominal, nominal itself skipped
    std::vector<RatioResults> all_results;
    auto allKeys = GetPrimaryVariantKeys(CurrentChannel());
    int colorIdx = 0;
    const int colors[] = {kBlue, kRed, kGreen + 2, kViolet, kOrange + 7, kAzure + 1, kMagenta + 1};
    const int markers[] = {20, 21, 22, 23, 29, 33, 34};
    const int nColors = sizeof(colors) / sizeof(colors[0]);

    for (const auto& key : allKeys) {
        if (key == CurrentChannel()) continue;
        // cosmetic, for the legend / summary plots only.
        std::string short_key = key; // "DB_bJER10v00" -> "bJER10v00"
        std::string prefix = CurrentChannel() + "_";
        if (short_key.rfind(prefix, 0) == 0) {
            short_key = short_key.substr(prefix.size());
        }
        std::string label = short_key + "/Nominal";
        all_results.push_back(CalculateRatios(CurrentChannel(), key, wp, label, colors[colorIdx % nColors],
                                              markers[colorIdx % nColors], ptCfg));
        ++colorIdx;
    }

    if (all_results.empty()) {
        std::cerr << "No non-nominal variants found - nothing to plot.\n";
        return;
    }

    // Accessors for the value-error combinations
    auto valJMS = [](const RatioResults& r, int b) {
        return r.mean_ratio[b];
    };
    auto errJMSu = [](const RatioResults& r, int b) {
        return r.mean_ratio_err[b];
    };
    auto errJMSc = [](const RatioResults& r, int b) {
        return r.mean_ratio_err_corr[b];
    };
    auto errJMSf = [](const RatioResults& r, int b) {
        return r.mean_ratio_err_full[b];
    };
    auto valJMR = [](const RatioResults& r, int b) {
        return r.sigma_ratio[b];
    };
    auto errJMRu = [](const RatioResults& r, int b) {
        return r.sigma_ratio_err[b];
    };
    auto errJMRc = [](const RatioResults& r, int b) {
        return r.sigma_ratio_err_corr[b];
    };
    auto errJMRf = [](const RatioResults& r, int b) {
        return r.sigma_ratio_err_full[b];
    };

    DrawCombinedRatioPlot_vs_pT(
        all_results, ptCfg, valJMS, errJMSu, "c_jms_uncorr",
        "#frac{b-JMDP_{Corrected}^{Data} / b-JMDP_{Corrected}^{MC}}{b-JMDP_{Nominal}^{Data} / b-JMDP_{Nominal}^{MC}}",
        outputPath + "Combined_Mass_Scale_Ratios.pdf", wp);

    // JMS - correlation-aware (Rho = 1 assumed)
    DrawCombinedRatioPlot_vs_pT(
        all_results, ptCfg, valJMS, errJMSc, "c_jms_corr",
        "#frac{b-JMDP_{Corrected}^{Data} / b-JMDP_{Corrected}^{MC}}{b-JMDP_{Nominal}^{Data} / b-JMDP_{Nominal}^{MC}} ",
        outputPath + "Combined_Mass_Scale_Ratios_CorrAware.pdf", wp);

    // JMS - full correlation matrix (Rho read from correlations.txt)
    DrawCombinedRatioPlot_vs_pT(
        all_results, ptCfg, valJMS, errJMSf, "c_jms_full",
        "#frac{b-JMDP_{Corrected}^{Data} / b-JMDP_{Corrected}^{MC}}{b-JMDP_{Nominal}^{Data} / b-JMDP_{Nominal}^{MC}} (full-corr)",
        outputPath + "Combined_Mass_Scale_Ratios_FullCorr.pdf", wp);

    // JMS - error-treatment comparison (uncorr vs corr-aware)
    DrawErrorComparisonPlot_vs_pT(all_results, ptCfg, valJMS, errJMSu, errJMSc, "c_jms_cmp",
                                  "DR(JMS) - error-treatment comparison",
                                  outputPath + "Combined_Mass_Scale_Ratios_ErrorComparison.pdf", wp);

    // ---------- JMR plots ----------

    // JMR - uncorrelated
    DrawCombinedRatioPlot_vs_pT(all_results, ptCfg, valJMR, errJMRu, "c_jmr_uncorr",
                                "b-JMDR_{Corrected}/b-JMDR_{Nominal}", outputPath + "Combined_JMR_Ratios.pdf",
                                wp);

    // JMR - correlation-aware (Rho = 1 assumed)
    DrawCombinedRatioPlot_vs_pT(all_results, ptCfg, valJMR, errJMRc, "c_jmr_corr",
                                "b-JMDR_{Corrected}/b-JMDR_{Nominal} ",
                                outputPath + "Combined_JMR_Ratios_CorrAware.pdf", wp);

    // JMR - full correlation matrix (Rho read from correlations.txt)
    DrawCombinedRatioPlot_vs_pT(all_results, ptCfg, valJMR, errJMRf, "c_jmr_full",
                                "b-JMDR_{Corrected}/b-JMDR_{Nominal} (full-corr)",
                                outputPath + "Combined_JMR_Ratios_FullCorr.pdf", wp);

    // JMR - error-treatment comparison (uncorr vs corr-aware)
    DrawErrorComparisonPlot_vs_pT(all_results, ptCfg, valJMR, errJMRu, errJMRc, "c_jmr_cmp",
                                  "DR(JMR) - error-treatment comparison",
                                  outputPath + "Combined_JMR_Ratios_ErrorComparison.pdf", wp);

    // ---- 2-D heatmaps ----

    // JMS - uncorrelated
    Draw2DRatioMap(all_results, ptCfg, valJMS, errJMSu, "c_2d_jms_uncorr", " b-JMDP Ratios",
                   outputPath + "Combined_2D_Mass_Scale_Map.pdf", wp);

    // JMS - correlation-aware
    Draw2DRatioMap(all_results, ptCfg, valJMS, errJMSc, "c_2d_jms_corr", " b-JMDP Ratios ",
                   outputPath + "Combined_2D_Mass_Scale_Map_CorrAware.pdf", wp);

    // JMS - full correlation matrix
    Draw2DRatioMap(all_results, ptCfg, valJMS, errJMSf, "c_2d_jms_full", " b-JMDP Ratios (full-corr)",
                   outputPath + "Combined_2D_Mass_Scale_Map_FullCorr.pdf", wp);

    // JMR - uncorrelated
    Draw2DRatioMap(all_results, ptCfg, valJMR, errJMRu, "c_2d_jmr_uncorr", " b-JMDR Ratios",
                   outputPath + "Combined_2D_JMR_Map.pdf", wp);

    // JMR - correlation-aware
    Draw2DRatioMap(all_results, ptCfg, valJMR, errJMRc, "c_2d_jmr_corr", " b-JMDR Ratios ",
                   outputPath + "Combined_2D_JMR_Map_CorrAware.pdf", wp);

    // JMR - full correlation matrix
    Draw2DRatioMap(all_results, ptCfg, valJMR, errJMRf, "c_2d_jmr_full", " b-JMDR Ratios (full-corr)",
                   outputPath + "Combined_2D_JMR_Map_FullCorr.pdf", wp);

    // ---- Absolute per-variant plots  ----
    auto valJMSAbs = [](const RatioResults& r, int b) {
        return r.jms_var[b];
    };
    auto errJMSAbs = [](const RatioResults& r, int b) {
        return r.jms_var_err[b];
    };
    auto valJMSAbsNom = [](const RatioResults& r, int b) {
        return r.jms_nom[b];
    };
    auto errJMSAbsNom = [](const RatioResults& r, int b) {
        return r.jms_nom_err[b];
    };
    auto valJMRAbs = [](const RatioResults& r, int b) {
        return r.jmr_var[b];
    };
    auto errJMRAbs = [](const RatioResults& r, int b) {
        return r.jmr_var_err[b];
    };
    auto valJMRAbsNom = [](const RatioResults& r, int b) {
        return r.jmr_nom[b];
    };
    auto errJMRAbsNom = [](const RatioResults& r, int b) {
        return r.jmr_nom_err[b];
    };

    DrawAbsolutePlot_vs_pT(all_results, ptCfg, valJMSAbs, errJMSAbs, valJMSAbsNom, errJMSAbsNom, "c_jms_abs",
                           "b-JMDP = #mu_{data} / #mu_{MC}", outputPath + "Combined_Mass_Scale_Absolute.pdf",
                           wp);

    DrawAbsolutePlot_vs_pT(all_results, ptCfg, valJMRAbs, errJMRAbs, valJMRAbsNom, errJMRAbsNom, "c_jmr_abs",
                           "b-JMDR = #sigma_{data} / #sigma_{MC}",
                           outputPath + "Combined_Mass_Width_Absolute.pdf", wp);

    // Summary plots (DrawDualAxisSummary) for b-JMDP and b-JMDR: single ratio
    // and double ratio with the corr-aware error, one per variant, one with
    // all variants, and one for the EASYBJETS_SUMMARY_VARIANTS subset.
    // Written to Summary/Summary_{Scale,Width}_<variant|AllVariants|Subset_...>.pdf
    //
    std::string summaryDir = outputPath + "Summary/";
    gSystem->Exec(("mkdir -p " + summaryDir).c_str());

    // Accessors: single ratios and double ratios for scale (JMS) and width (JMR)
    auto sJMSvar = [](const RatioResults& r, int b) {
        return r.jms_var[b];
    };
    auto sJMSvarE = [](const RatioResults& r, int b) {
        return r.jms_var_err[b];
    };
    auto sJMSnom = [](const RatioResults& r, int b) {
        return r.jms_nom[b];
    };
    auto sJMSnomE = [](const RatioResults& r, int b) {
        return r.jms_nom_err[b];
    };
    auto drJMS = [](const RatioResults& r, int b) {
        return r.mean_ratio[b];
    };
    auto drJMSE = [](const RatioResults& r, int b) {
        return r.mean_ratio_err_corr[b];
    };

    auto sJMRvar = [](const RatioResults& r, int b) {
        return r.jmr_var[b];
    };
    auto sJMRvarE = [](const RatioResults& r, int b) {
        return r.jmr_var_err[b];
    };
    auto sJMRnom = [](const RatioResults& r, int b) {
        return r.jmr_nom[b];
    };
    auto sJMRnomE = [](const RatioResults& r, int b) {
        return r.jmr_nom_err[b];
    };
    auto drJMR = [](const RatioResults& r, int b) {
        return r.sigma_ratio[b];
    };
    auto drJMRE = [](const RatioResults& r, int b) {
        return r.sigma_ratio_err_corr[b];
    };

    // ---- (1) One plot per variant ----
    for (size_t k = 0; k < all_results.size(); ++k) {
        std::string vlabel = all_results[k].label;
        std::string suffix = "/Nominal";
        if (vlabel.size() >= suffix.size()
            && vlabel.compare(vlabel.size() - suffix.size(), suffix.size(), suffix) == 0)
            vlabel = vlabel.substr(0, vlabel.size() - suffix.size());

        std::vector<size_t> one = {k};
        DrawDualAxisSummary(all_results, ptCfg, sJMSvar, sJMSvarE, sJMSnom, sJMSnomE, drJMS, drJMSE,
                            ("c_sum_jms_" + vlabel), "b-JMDP = #mu_{data} / #mu_{MC}",
                            "DR (b-JMDP_{var}/b-JMDP_{nom})", summaryDir + "Summary_Scale_" + vlabel + ".pdf",
                            wp, one);
        DrawDualAxisSummary(all_results, ptCfg, sJMRvar, sJMRvarE, sJMRnom, sJMRnomE, drJMR, drJMRE,
                            ("c_sum_jmr_" + vlabel), "b-JMDR = #sigma_{data} / #sigma_{MC}",
                            "DR (b-JMDR_{var}/b-JMDR_{nom})", summaryDir + "Summary_Width_" + vlabel + ".pdf",
                            wp, one);
    }

    // ---- (2) One plot for ALL variants ----
    DrawDualAxisSummary(all_results, ptCfg, sJMSvar, sJMSvarE, sJMSnom, sJMSnomE, drJMS, drJMSE,
                        "c_sum_jms_all", "b-JMDP = #mu_{data} / #mu_{MC}", "DR (b-JMDP_{var}/b-JMDP_{nom})",
                        summaryDir + "Summary_Scale_AllVariants.pdf", wp);
    DrawDualAxisSummary(all_results, ptCfg, sJMRvar, sJMRvarE, sJMRnom, sJMRnomE, drJMR, drJMRE,
                        "c_sum_jmr_all", "b-JMDR = #sigma_{data} / #sigma_{MC}",
                        "DR (b-JMDR_{var}/b-JMDR_{nom})", summaryDir + "Summary_Width_AllVariants.pdf", wp);

    // (3) Subset chosen with EASYBJETS_SUMMARY_VARIANTS, a comma-separated
    // list of variant names without channel, e.g.
    //
    //   export EASYBJETS_SUMMARY_VARIANTS=bJER10v01,MuonCorr
    //
    // If no name matches, a warning is printed and no subset plot is made.
    {
        const char* envSub = std::getenv("EASYBJETS_SUMMARY_VARIANTS");
        if (envSub && *envSub) {
            // parse comma list
            std::vector<std::string> want;
            std::string s(envSub), tok;
            std::stringstream ss(s);
            while (std::getline(ss, tok, ',')) {
                // trim whitespace
                size_t a = tok.find_first_not_of(" \t");
                size_t b = tok.find_last_not_of(" \t");
                if (a != std::string::npos) want.push_back(tok.substr(a, b - a + 1));
            }
            // resolve indices
            std::vector<size_t> sel;
            std::string tag;
            for (size_t k = 0; k < all_results.size(); ++k) {
                std::string vlabel = all_results[k].label;
                std::string suffix = "/Nominal";
                if (vlabel.size() >= suffix.size()
                    && vlabel.compare(vlabel.size() - suffix.size(), suffix.size(), suffix) == 0)
                    vlabel = vlabel.substr(0, vlabel.size() - suffix.size());
                for (const auto& w : want) {
                    if (vlabel == w) {
                        sel.push_back(k);
                        tag += (tag.empty() ? "" : "_") + vlabel;
                        break;
                    }
                }
            }
            if (!sel.empty()) {
                DrawDualAxisSummary(all_results, ptCfg, sJMSvar, sJMSvarE, sJMSnom, sJMSnomE, drJMS, drJMSE,
                                    "c_sum_jms_subset", "b-JMDP = #mu_{data} / #mu_{MC}",
                                    "DR (b-JMDP_{var}/b-JMDP_{nom})",
                                    summaryDir + "Summary_Scale_Subset_" + tag + ".pdf", wp, sel);
                DrawDualAxisSummary(all_results, ptCfg, sJMRvar, sJMRvarE, sJMRnom, sJMRnomE, drJMR, drJMRE,
                                    "c_sum_jmr_subset", "b-JMDR = #sigma_{data} / #sigma_{MC}",
                                    "DR (b-JMDR_{var}/b-JMDR_{nom})",
                                    summaryDir + "Summary_Width_Subset_" + tag + ".pdf", wp, sel);
                std::cout << "[Mass_Ratios_Width] Custom-subset summary plots written for: " << tag << "\n";
            } else {
                std::cerr << "[Mass_Ratios_Width] WARN: EASYBJETS_SUMMARY_VARIANTS='" << envSub
                          << "' matched no variants; subset plot skipped.\n";
            }
        }
    }

    // ---- Summary text ----
    WriteCombinedSummary(all_results, ptCfg, outputPath + "Combined_Mass_Ratio_Summary.txt", wp);
    WriteMassRatiosCSV(all_results, ptCfg, outputPath + "Combined_Mass_Ratio_Summary.csv", CurrentChannel(),
                       PathBuilder::WPtoFileWP(wp));

    std::cout << "\nCombined analysis complete. Outputs in: " << outputPath << "\n";
}


// Entry point, called by Mass_Ratios_PeakWidth.cxx in each channel's macros/
// folder (the folder sets DB or MJB):
//
//   DB : cd ZbbGamma/macros            |       MJB: cd ZbbJets/macros
//   root 'Mass_Ratios_PeakWidth.cxx("inclusive")' # one pT bin
//   root 'Mass_Ratios_PeakWidth.cxx("ptbins")'    # DB pT bins
//   root 'Mass_Ratios_PeakWidth.cxx("ALL")'       # both
//
// Run Mass_Fit.cxx first with the same pT mode and EASYBJETS_BINNING.
inline void RunOnePtMode_Width(const std::string& modeLabel, const PtBinConfig::PtBinDef& ptCfg) {
    std::cout << "\n*******************************************\n";
    std::cout << "COMBINED MASS RATIO ANALYSIS - WIDTH RATIO version\n";
    std::cout << "Comparing all calibrations to Nominal\n";
    std::cout << "JMR = sigma_data / sigma_MC   (simple width ratio,\n";
    std::cout << "      dimensionless, NOT the signed-J definition)\n";
    std::cout << "pT mode: " << modeLabel << " (" << ptCfg.nBins << " bins)\n";
    std::cout << "Correlation method: " << CorrMethodTag() << "\n";
    std::cout << "*******************************************\n\n";

    auto workingPoints = AnalysisWPConfig::StandardWPs();
    for (const auto& wp : workingPoints) PerformCombinedMassRatioWidthAnalysis(wp, ptCfg);
}

// Runs the given pT mode for the current channel:
//   "inclusive" (or "1bin")  one pT bin over the full channel range
//   "ptbins"                 the channel's pT bins (PtBinConfig in VariantConfig.h)
//   "ALL"                    inclusive, then ptbins
// The mode is case-insensitive; any other value prints an error and stops.
inline void Mass_Ratios_Width_Run(const std::string& ptMode, const std::string& channel) {
    CurrentChannel() = channel; // establish channel state
    gROOT->SetBatch(kTRUE);
    AtlasUtils::LoadAtlasStyle();

    std::string mode = ptMode;
    for (auto& c : mode) c = (char)std::tolower((unsigned char)c);

    if (mode == "inclusive" || mode == "1bin") {
        RunOnePtMode_Width("inclusive", PtBinConfig::InclusiveBins(channel));
    } else if (mode == "ptbins") {
        RunOnePtMode_Width("ptbins", PtBinConfig::DefaultBins(channel));
    } else if (mode == "all") {
        std::cout << "\n[Mass_Ratios_Width] mode=ALL -> running inclusive, then ptbins\n";
        RunOnePtMode_Width("inclusive", PtBinConfig::InclusiveBins(channel));
        RunOnePtMode_Width("ptbins", PtBinConfig::DefaultBins(channel));
    } else {
        std::cerr << "\n[Mass_Ratios_Width] ERROR: unknown pT mode '" << ptMode << "'\n"
                  << "               Valid modes: \"inclusive\", \"ptbins\", \"ALL\"\n\n";
        return;
    }

    std::cout << "\nAll width-ratio analyses completed.\n";
}

#endif // MASS_RATIOS_WIDTH_H