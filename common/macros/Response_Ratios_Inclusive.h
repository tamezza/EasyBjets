/**
 * ============================================================================
 *  Response_Ratios_Inclusive.h - inclusive R_JES / R_JER double ratios
 * ============================================================================
 * ===========================================
 * @author:    Thiziri AMEZZA
 * @email:     thiziri.amezza@cern.ch
 * @date:      13 April 2026
 * ===========================================
 *
 * Shared code for the DB (ZbbGamma) and MJB (ZbbJets) channels. Each channel
 * has a Response_Ratios_Inclusive.cxx that calls
 * Ratio_Response_Inclusive_Run(variant, "DB") or (variant, "MJB").
 *
 * It reads the inclusive response fits (one fit over the whole pT range, no
 * mass bins) of the nominal and of each correction (bJER10v00, bJER10v01,
 * MuonCorr) and builds, per correction:
 *
 *   R_JES           = <R>_data / <R>_combMC         (R = R_DB or R_MJB)
 *   sigma_insitu    = sqrt(sigma_data^2   - sigma_truth^2)    
 *   sigma^MC_insitu = sqrt(sigma_combMC^2 - sigma_truth^2)
 *   R_JER           = sigma_insitu / sigma^MC_insitu
 *
 * sigma_truth is the width of the truth/reference fit in combined MC
 * (truth_ref). Then the double ratios correction / nominal:
 *
 *   DR(R_JES), DR(sigma_insitu), DR(R_JER)       DR(X) = X_corr / X_nom
 *
 * DR(sigma_insitu) is a data-only cross-check; DR(R_JER) is the real result.
 *
 * Error treatments - each double ratio is given with three errors:
 *
 *   uncorrelated   All fit errors (data and MC, correction and nominal) are
 *                  treated as independent and added in quadrature. This ignores
 *                  that correction and nominal use the same events, so it is the
 *                  largest error: a conservative upper bound (we overestimate the error).
 *
 *   corr-aware     Correction and nominal fluctuate together (same events), so
 *                  much of the error cancels in the double ratio:
 *                    R_JES       rho = 1 between correction and nominal, with one
 *                                common error per sample (geometric mean of the
 *                                two fit errors). Close to 0 if the correction does
 *                                not move the peak.
 *                    resolution  the sigma_data and sigma_MC errors are dropped as
 *                                fully cancelling; only the sigma_truth fit error
 *                                is kept (independent between the two).
 *                  An approximation of the rho = 1 case that can underestimate
 *                  the error.
 *
 *   full           Uses the correlations measured by the bootstrap
 *                  (correlations.txt) for R_JES, sigma_insitu and R_JER.
 *                  If none is found (e.g. MuonCorr), rho = 1 is used.
 * 
 *                  Expected: corr-aware <= full <= uncorrelated.
 * 
 *
 * How to run (from ZbbGamma/macros or ZbbJets/macros, after Response_Fit_Inclusive.cxx):
 * 
 *   export EASYBJETS_RUN=run3          same as for the fit  (run2 or run3)
 *   root 'Response_Ratios_Inclusive.cxx()'                  all corrections
 *   root 'Response_Ratios_Inclusive.cxx("DB_bJER10v00")'    one correction (DB)
 *   root 'Response_Ratios_Inclusive.cxx("MJB_MuonCorr")'    one correction (MJB)
 * 
 * Then run Response_Summary_Combined_Inclusive.cxx.
 *
 * Input: <PlotDir>/<DB|MJB>/FitResponse_<path_tag>/<wp>/{data,combined,truth_ref}_fit_parameters.txt
 *        with PlotDir = ../Run3Plots or ../Run2Plots (plus _<model> for a
 *        non-default EASYBJETS_FIT_MODEL), path_tag = nominal, bJER10v00, ...
 * 
 * Output: <PlotDir>/Response/Ratios/<DB|MJB>_Inclusive_Ratio/<wp>/
 *   Summary_Inclusive_Ratios_<wp>.txt     table of all double ratios and fit inputs
 *   Summary_Inclusive_Ratios_<wp>.csv     one row per correction, 44 columns, read
 *                                         by Response_Summary_Combined_Inclusive.cxx
 * 
 *   R_JES_Inclusive_DoubleRatio_<wp>.pdf
 *   sigma_insitu_Inclusive_DoubleRatio_<wp>{,_CorrAware,_ErrorComparison}.pdf
 *   R_JER_Inclusive_DoubleRatio_<wp>{,_CorrAware,_ErrorComparison}.pdf
 *   R_JER_Inclusive_Absolute_<wp>{,_CorrAware}.pdf   R_JER per correction, not a DR
 *   R_JES_Inclusive_Absolute_<wp>{,_CorrAware}.pdf   R_JES per correction, not a DR
 *   Summary/Summary_{Scale,Resolution}_<label>.pdf   data/MC + DR, one correction
 *   Summary/Summary_{Scale,Resolution}_AllVariants.pdf
 *   Summary/Summary_{Scale,Resolution}_Subset_<labels>.pdf   only with EASYBJETS_SUMMARY_VARIANTS
 *
 * Options (environment variables):
 *   EASYBJETS_RUN               run2 | run3 (default run3); picks Run2Plots or Run3Plots
 *   EASYBJETS_FIT_MODEL         empty = Gaussian (default); bukin, cb, dscb, dgauss
 *                               or all add _<tag> to the plot folder, so it must match the fit
 *   EASYBJETS_CORR_METHOD       rho1 (default) | file | bootstrap; only changes the
 *                               "full" R_JES error
 *   EASYBJETS_CORRELATIONS      path of the rho file (default ./correlations.txt),
 *                               read only when the method is file or bootstrap
 *   EASYBJETS_SUMMARY_VARIANTS  comma-separated correction labels for the extra
 *                               Summary_*_Subset plots (default: none)
 * 
 * EASYBJETS_BINNING and EASYBJETS_MASS_CUT are not read here. They change the
 * fit results, but the fit writes to the same folder whatever their value, so
 * run this macro right after the fit you want.
 *
 * ==========================================================================
 * @version:  1.0
 * ==========================================================================
 */

#ifndef RATIO_RESPONSE_INCLUSIVE_H
#define RATIO_RESPONSE_INCLUSIVE_H

#include "VariantConfig.h"
#include "CommonUtils.h"
#include "CorrelationConfig.h"

#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>
#include <functional>
#include "TGaxis.h"

// Reader for the <type>_fit_parameters.txt files written by
// FitOutputUtils::SaveFitParameters ("Mean: x +/- dx", "Sigma: s +/- ds").
struct FitParams {
    double mean = 0.0;
    double mean_err = 0.0;
    double sigma = 0.0;
    double sigma_err = 0.0;
    bool valid = false;
};

inline bool ReadFitParams(const std::string& filepath, FitParams& params) {
    std::ifstream infile(filepath);
    if (!infile.is_open()) {
        Log::warning("Ratio_Response_Inclusive", "Cannot open: " + filepath);
        return false;
    }

    std::string line;
    while (std::getline(infile, line)) {
        auto parseLine = [](const std::string& l, double& val, double& err) {
            size_t colon = l.find(':');
            if (colon == std::string::npos) return false;
            std::string rest = l.substr(colon + 1);

            size_t sep = rest.find("+/-");
            size_t sep_len = 3;
            if (sep == std::string::npos) {
                sep = rest.find("±");
                sep_len = 2; // UTF-8 "±" is 2 bytes
            }
            if (sep == std::string::npos) {
                val = std::atof(rest.c_str());
                err = 0.0;
                return true;
            }

            std::string val_str = rest.substr(0, sep);
            std::string err_str = rest.substr(sep + sep_len);
            val = std::atof(val_str.c_str());
            err = std::atof(err_str.c_str());
            return true;
        };

        if (line.substr(0, 5) == "Mean:") {
            parseLine(line, params.mean, params.mean_err);
        } else if (line.substr(0, 6) == "Sigma:") {
            parseLine(line, params.sigma, params.sigma_err);
        }
    }
    infile.close();
    params.valid = true;
    return true;
}

// sigma_insitu = sqrt(sigma_data^2 - sigma_truth^2), and its error
// with sigma_data and sigma_truth taken as independent. 
inline double JERFromEq92(double sigma_data, double sigma_truth) {
    double diff = sigma_data * sigma_data - sigma_truth * sigma_truth;
    return (diff > 0) ? std::sqrt(diff) : 0.0;
}

inline double JERErrFromEq92(double sigma_data, double sigma_data_err, double sigma_truth,
                             double sigma_truth_err) {
    double jer = JERFromEq92(sigma_data, sigma_truth);
    if (jer <= 0) return 0.0;
    double t1 = sigma_data * sigma_data_err / jer;
    double t2 = sigma_truth * sigma_truth_err / jer;
    return std::sqrt(t1 * t1 + t2 * t2);
}

// Corr-aware error on DR(sigma_insitu) = sigma_insitu_corr / sigma_insitu_nom.
inline double SigmaInsituDoubleRatioErrCorrAware(double jer_var, double sigma_truth_var,
                                                 double sigma_truth_var_err, double jer_nom,
                                                 double sigma_truth_nom, double sigma_truth_nom_err) {
    if (jer_var <= 0 || jer_nom <= 0) return 0.0;
    double ratio = jer_var / jer_nom;

    double rel_var = (sigma_truth_var * sigma_truth_var_err) / (jer_var * jer_var);
    double rel_nom = (sigma_truth_nom * sigma_truth_nom_err) / (jer_nom * jer_nom);

    return std::abs(ratio) * std::sqrt(rel_var * rel_var + rel_nom * rel_nom);
}

// Corr-aware error on DR(R_JER) = R_JER_corr / R_JER_nom.
// With R = sqrt((s_d^2 - s_t^2) / (s_m^2 - s_t^2)):
//   dR/ds_t = s_t (s_d^2 - s_m^2) / (sigma_insitu * sigma^MC_insitu^3)
// so the relative part from the truth width is
//   dR/R |s_t = |s_t (s_d^2 - s_m^2)| * ds_t / ((s_d^2 - s_t^2)(s_m^2 - s_t^2))
// The s_d and s_m errors are left out (same events), and
//   (dDR/DR)^2 = (dR_corr/R_corr)^2|s_t + (dR_nom/R_nom)^2|s_t
inline double RJERRelErrFromSigmaTruth(double sigma_d, double sigma_m, double sigma_t, double sigma_t_err) {
    double u = sigma_d * sigma_d - sigma_t * sigma_t; // sigma_insitu^2
    double v = sigma_m * sigma_m - sigma_t * sigma_t; // sigma^MC_insitu^2
    if (u <= 0 || v <= 0) return 0.0;
    double num = std::abs(sigma_t * (sigma_d * sigma_d - sigma_m * sigma_m));
    return num * sigma_t_err / (u * v);
}

inline double RJERDoubleRatioErrCorrAware(double R_var, double sigma_d_var, double sigma_m_var,
                                          double sigma_t_var, double sigma_t_var_err, double R_nom,
                                          double sigma_d_nom, double sigma_m_nom, double sigma_t_nom,
                                          double sigma_t_nom_err) {
    if (R_var <= 0 || R_nom <= 0) return 0.0;
    double DR = R_var / R_nom;
    double rel_var = RJERRelErrFromSigmaTruth(sigma_d_var, sigma_m_var, sigma_t_var, sigma_t_var_err);
    double rel_nom = RJERRelErrFromSigmaTruth(sigma_d_nom, sigma_m_nom, sigma_t_nom, sigma_t_nom_err);
    return std::abs(DR) * std::sqrt(rel_var * rel_var + rel_nom * rel_nom);
}

struct DoubleRatio {
    double value = 0.0;
    double err = 0.0;
};

// DR(R_JES) = (<R>_data / <R>_combMC)_corr / (<R>_data / <R>_combMC)_nom with
// the four means taken as independent. This ignores that correction and
// nominal share events, so it is an upper bound on the error.
inline DoubleRatio ComputeJESDoubleRatio(double d_var, double d_var_err, double c_var, double c_var_err,
                                         double d_nom, double d_nom_err, double c_nom, double c_nom_err) {
    DoubleRatio r;
    if (c_var == 0 || c_nom == 0 || d_nom == 0) return r;

    double ratio_var = d_var / c_var;
    double ratio_nom = d_nom / c_nom;
    if (ratio_nom == 0) return r;

    r.value = ratio_var / ratio_nom;

    double rel_d_var = (d_var != 0) ? d_var_err / d_var : 0;
    double rel_c_var = c_var_err / c_var;
    double rel_d_nom = (d_nom != 0) ? d_nom_err / d_nom : 0;
    double rel_c_nom = c_nom_err / c_nom;

    double rel_total = std::sqrt(rel_d_var * rel_d_var + rel_c_var * rel_c_var + rel_d_nom * rel_d_nom
                                 + rel_c_nom * rel_c_nom);
    r.err = std::abs(r.value) * rel_total;
    return r;
}

// Corr-aware DR(R_JES), rho = 1 between correction and nominal.
//   DR = (mu_d_corr * mu_m_nom) / (mu_m_corr * mu_d_nom),  mu = <pT_jet/pT_ref>
// The corrections change the jet pT of the same events, so mu_d_corr and
// mu_d_nom move together; the same holds in MC. Data and MC are independent.
// A common shift eps_d of both data means gives
//   d(ln DR)/d(eps_d) = 1/mu_d_corr - 1/mu_d_nom   (bracket_d)
//   d(ln DR)/d(eps_m) = 1/mu_m_nom  - 1/mu_m_corr  (bracket_m)
//   (dDR/DR)^2 = bracket_d^2 * dmu_d^2 + bracket_m^2 * dmu_m^2
// with dmu = sqrt(err_corr * err_nom), the geometric mean of the two fit
// errors. Same method as JMSDoubleRatioErrCorrAware in Mass_Ratios_PeakWidth.h.
// The error goes to 0 when the correction does not move the peak.
inline DoubleRatio ComputeJESDoubleRatioCorrAware(double d_var, double d_var_err, double c_var,
                                                  double c_var_err, double d_nom, double d_nom_err,
                                                  double c_nom, double c_nom_err) {
    DoubleRatio r;
    if (c_var == 0 || c_nom == 0 || d_nom == 0 || d_var == 0) return r;

    double ratio_var = d_var / c_var;
    double ratio_nom = d_nom / c_nom;
    if (ratio_nom == 0) return r;
    r.value = ratio_var / ratio_nom;

    double dd_common = std::sqrt(std::abs(d_var_err * d_nom_err));
    double dc_common = std::sqrt(std::abs(c_var_err * c_nom_err));

    double bracket_d = (1.0 / d_var) - (1.0 / d_nom);
    double bracket_c = (1.0 / c_nom) - (1.0 / c_var);

    double rel2 =
        bracket_d * bracket_d * dd_common * dd_common + bracket_c * bracket_c * dc_common * dc_common;

    r.err = std::abs(r.value) * std::sqrt(rel2);
    return r;
}

// Full-correlation DR(R_JES): c^T Sigma c with any 4x4 correlation matrix.
// Inputs (mu_d_corr, mu_d_nom, mu_m_corr, mu_m_nom); the rho values come from
// the response ("resp") lines of correlations.txt via CorrelationConfig.
//   c = (1/mu_d_corr, -1/mu_d_nom, -1/mu_m_corr, 1/mu_m_nom)
//   (dDR/DR)^2 = sum_ij c_i c_j rho_ij err_i err_j
// With the default matrix (rho_dd = rho_mm = 1, data-MC = 0) this gives
// (err_d_corr/mu_d_corr - err_d_nom/mu_d_nom)^2 + (same for MC). That equals
// the corr-aware result above only when the correction and nominal fit
// errors are equal.
inline DoubleRatio ComputeJESDoubleRatioFullCorr(double d_var, double d_var_err, double c_var,
                                                 double c_var_err, double d_nom, double d_nom_err,
                                                 double c_nom, double c_nom_err,
                                                 const CorrelationConfig::CorrelationMatrix& C) {
    DoubleRatio r;
    if (c_var == 0 || c_nom == 0 || d_nom == 0 || d_var == 0) return r;

    double ratio_var = d_var / c_var;
    double ratio_nom = d_nom / c_nom;
    if (ratio_nom == 0) return r;
    r.value = ratio_var / ratio_nom;

    // d(ln DR) / d(input)
    double c_vec[4] = {
        1.0 / d_var,  // mu_d corr
        -1.0 / d_nom, // mu_d nom
        -1.0 / c_var, // mu_m corr
        1.0 / c_nom   // mu_m nom
    };
    double d_vec[4] = {d_var_err, d_nom_err, c_var_err, c_nom_err};

    // Correlation matrix from the mu fields of C.
    double rho[4][4] = {{1.0, C.rho_mu_dd, C.rho_mu_dv_mv, C.rho_mu_dv_mn},
                        {C.rho_mu_dd, 1.0, C.rho_mu_dn_mv, C.rho_mu_dn_mn},
                        {C.rho_mu_dv_mv, C.rho_mu_dn_mv, 1.0, C.rho_mu_mm},
                        {C.rho_mu_dv_mn, C.rho_mu_dn_mn, C.rho_mu_mm, 1.0}};

    double rel2 = 0.0;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j) rel2 += c_vec[i] * c_vec[j] * rho[i][j] * d_vec[i] * d_vec[j];

    rel2 = std::max(0.0, rel2); // rounding can make it slightly negative when rho is close to 1
    r.err = std::abs(r.value) * std::sqrt(rel2);
    return r;
}

// Generic double ratio of two quantities with independent errors.
inline DoubleRatio ComputeDoubleRatio(double num, double num_err, double den, double den_err) {
    DoubleRatio r;
    if (den == 0) return r;
    r.value = num / den;
    double rel_num = (num != 0) ? num_err / num : 0;
    double rel_den = den_err / den;
    r.err = std::abs(r.value) * std::sqrt(rel_num * rel_num + rel_den * rel_den);
    return r;
}

// Everything computed for one correction at one WP.
struct VariantRatio {
    std::string key;   // e.g. "DB_bJER10v00"
    std::string label; // legend label from VariantConfig, e.g. "bJR10v00", "#mu-in-jet"

    // Double ratios correction / nominal
    DoubleRatio jes;             // R_JES, uncorrelated
    DoubleRatio jes_corr;        // R_JES, rho = 1
    DoubleRatio jes_full;        // R_JES, rho from CorrelationConfig
    DoubleRatio sig_insitu;      // sigma_insitu, uncorrelated
    DoubleRatio sig_insitu_corr; // sigma_insitu, corr-aware (sigma_truth error only)
    DoubleRatio rjer;            // R_JER, uncorrelated
    DoubleRatio rjer_corr;       // R_JER, corr-aware (sigma_truth error only)
    DoubleRatio sig_insitu_full; // sigma_insitu, rho from CorrelationConfig (six widths)
    DoubleRatio rjer_full;       // R_JER, rho from CorrelationConfig (six widths)

    // R_JER of the correction and of the nominal
    double R_JER_var = 0, R_JER_var_err = 0;
    double R_JER_nom = 0, R_JER_nom_err = 0;

    // Split of the R_JER error for the CorrAware absolute plot:
    //   err_unc^2 + err_corr^2 = err^2
    // err_unc comes from sigma_truth only; it is the part kept in the
    // corr-aware DR(R_JER). err_corr is the sigma_data + sigma_MC part.
    double R_JER_var_err_unc = 0;  // sigma_truth part
    double R_JER_var_err_corr = 0; // sigma_data and sigma_MC parts, in quadrature
    double R_JER_nom_err_unc = 0;
    double R_JER_nom_err_corr = 0;

    // R_JES = <R>_data / <R>_combMC of the correction and of the nominal (not
    // a double ratio), with (dR/R)^2 = (dmu_d/mu_d)^2 + (dmu_m/mu_m)^2.
    // Data and MC share no input, so err_unc = err and err_corr = 0. The
    // fields exist only so the R_JES CorrAware plot can reuse the R_JER code.
    double R_JES_var = 0, R_JES_var_err = 0;
    double R_JES_nom = 0, R_JES_nom_err = 0;
    double R_JES_var_err_unc = 0;
    double R_JES_var_err_corr = 0;
    double R_JES_nom_err_unc = 0;
    double R_JES_nom_err_corr = 0;

    // Fit inputs, for the summary table
    double sigma_data_var = 0, sigma_data_nom = 0;
    double sigma_comb_var = 0, sigma_comb_nom = 0;
    double sigma_truth_var = 0, sigma_truth_nom = 0;
    double mean_data_var = 0, mean_data_nom = 0;
    double mean_comb_var = 0, mean_comb_nom = 0;
    double jer_var_val = 0, jer_nom_val = 0;
    double jer_var_err = 0, jer_nom_err = 0;

    // b-JES and b-JER inputs for the CSV: b_JES = response mean (data, MC),
    // b_JER = sigma_insitu (data) and sigma^MC_insitu (MC).
    double bjes_data_var = 0, bjes_data_var_err = 0, bjes_mc_var = 0, bjes_mc_var_err = 0;
    double bjes_data_nom = 0, bjes_data_nom_err = 0, bjes_mc_nom = 0, bjes_mc_nom_err = 0;
    double bjer_data_var = 0, bjer_data_var_err = 0, bjer_mc_var = 0, bjer_mc_var_err = 0;
    double bjer_data_nom = 0, bjer_data_nom_err = 0, bjer_mc_nom = 0, bjer_mc_nom_err = 0;

    bool valid = false;
};

// Folder written by Response_Fit_Inclusive.cxx. Must match FitRdbDir() in
// Response_Fit_Inclusive.h.
inline std::string FitResponseDir(const VariantConfig& cfg) {
    return cfg.path_tag.empty() ? "FitResponse" : "FitResponse_" + cfg.path_tag;
}

inline std::string FitParamsPath(const VariantConfig& cfg, const std::string& wp_file,
                                 const std::string& type) {
    return (RunConfig::PlotDir() + "/") + cfg.plot_subdir + "/" + FitResponseDir(cfg) + "/" + wp_file + "/"
           + type + "_fit_parameters.txt";
}

// The three fits needed per variant and WP: data, combined MC, truth/ref.
struct VariantFitBundle {
    FitParams data, comb, truth_ref;
    bool valid = false;
};

inline VariantFitBundle ReadVariantBundle(const VariantConfig& cfg, const std::string& wp_file) {
    VariantFitBundle b;
    bool ok_d = ReadFitParams(FitParamsPath(cfg, wp_file, "data"), b.data);
    bool ok_c = ReadFitParams(FitParamsPath(cfg, wp_file, "combined"), b.comb);
    bool ok_t = ReadFitParams(FitParamsPath(cfg, wp_file, "truth_ref"), b.truth_ref);
    b.valid = ok_d && ok_c && ok_t;
    if (!b.valid) {
        Log::warning("Ratio_Response_Inclusive", "Incomplete fit bundle for " + cfg.key + " WP " + wp_file);
    }
    return b;
}

// Plot of an absolute quantity per correction (not a ratio to nominal). The
// nominal is the first point, at x = 0, with a dashed line at its value so
// the shifts can be read directly. Used for both R_JER and R_JES; the caller
// passes the nominal value and error.
inline void DrawRJERAbsoluteSummary(const std::vector<VariantRatio>& results,
                                    std::function<double(const VariantRatio&)> getValue,
                                    std::function<double(const VariantRatio&)> getError, double nomValue,
                                    double nomError, const std::string& canvasTag,
                                    const std::string& plotTitle, const std::string& yTitle,
                                    const std::string& outputPath, const std::string& wpLabel,
                                    Color_t nomColor, Color_t varColor) {
    const int nVar = (int)results.size();
    if (nVar == 0) return;
    const int n = nVar + 1; // nominal at x = 0

    TCanvas* c = new TCanvas(("c_" + canvasTag).c_str(), plotTitle.c_str(), 900, 600);
    c->SetLeftMargin(0.15);
    c->SetBottomMargin(0.20);
    c->SetTopMargin(0.10);
    c->SetRightMargin(0.05);
    c->SetGridy();

    // y range from all points including their errors
    std::vector<double> vals(n), errs(n);
    vals[0] = nomValue;
    errs[0] = std::abs(nomError);
    for (int i = 0; i < nVar; ++i) {
        vals[i + 1] = getValue(results[i]);
        errs[i + 1] = std::abs(getError(results[i]));
    }
    double ymin = 1e9, ymax = -1e9;
    for (int i = 0; i < n; ++i) {
        double lo = vals[i] - errs[i];
        double hi = vals[i] + errs[i];
        if (lo < ymin) ymin = lo;
        if (hi > ymax) ymax = hi;
    }
    double range = (ymax - ymin);
    if (range < 0.02) range = 0.05;
    double y_lo = ymin - 0.35 * range;
    double y_hi = ymax + 0.35 * range;

    TH1F* frame = new TH1F(("frame_" + canvasTag).c_str(), plotTitle.c_str(), n, -0.5, n - 0.5);
    frame->SetMinimum(y_lo);
    frame->SetMaximum(y_hi);
    frame->GetYaxis()->SetTitle(yTitle.c_str());
    frame->GetYaxis()->SetTitleSize(0.045);
    frame->GetYaxis()->SetTitleOffset(1.4);
    frame->GetYaxis()->SetLabelSize(0.040);
    frame->GetXaxis()->SetLabelSize(0.045);
    frame->GetXaxis()->SetTickLength(0);
    frame->SetStats(0);
    frame->GetXaxis()->SetBinLabel(1, "Nominal");
    for (int i = 0; i < nVar; ++i) frame->GetXaxis()->SetBinLabel(i + 2, results[i].label.c_str());
    frame->Draw("AXIS");

    // Dashed line at the nominal value
    TLine* ref = new TLine(-0.5, nomValue, n - 0.5, nomValue);
    ref->SetLineStyle(2);
    ref->SetLineColor(kGray + 2);
    ref->SetLineWidth(2);
    ref->Draw("same");

    // Nominal point, own graph so it gets its own marker and colour
    double x_nom[1] = {0.0};
    double y_nom[1] = {nomValue};
    double ex_nom[1] = {0.0};
    double ey_nom[1] = {std::abs(nomError)};
    TGraphErrors* g_nom = new TGraphErrors(1, x_nom, y_nom, ex_nom, ey_nom);
    g_nom->SetMarkerStyle(21); // filled square
    g_nom->SetMarkerSize(1.8);
    g_nom->SetMarkerColor(nomColor);
    g_nom->SetLineColor(nomColor);
    g_nom->SetLineWidth(2);
    g_nom->Draw("P same");

    // Variant points
    std::vector<double> xs(nVar), ys(nVar), exs(nVar, 0.0), eys(nVar);
    for (int i = 0; i < nVar; ++i) {
        xs[i] = i + 1; // x = 0 is the nominal
        ys[i] = getValue(results[i]);
        eys[i] = std::abs(getError(results[i]));
    }
    TGraphErrors* g_var = new TGraphErrors(nVar, xs.data(), ys.data(), exs.data(), eys.data());
    g_var->SetMarkerStyle(20); // filled circle
    g_var->SetMarkerSize(1.6);
    g_var->SetMarkerColor(varColor);
    g_var->SetLineColor(varColor);
    g_var->SetLineWidth(2);
    g_var->Draw("P same");

    // Value labels above each marker
    TLatex tx;
    tx.SetTextSize(0.028);
    tx.SetTextAlign(22);
    tx.DrawLatex(0, nomValue + std::abs(nomError) + 0.03 * range,
                 Form("%.3f #pm %.3f", nomValue, std::abs(nomError)));
    for (int i = 0; i < nVar; ++i) {
        tx.DrawLatex(i + 1, ys[i] + eys[i] + 0.03 * range, Form("%.3f #pm %.3f", ys[i], eys[i]));
    }

    AtlasUtils::DrawAtlasLabel(0.18, 0.85);
    TLatex wpTx;
    wpTx.SetNDC();
    wpTx.SetTextFont(42);
    wpTx.SetTextSize(0.035);
    wpTx.DrawLatex(0.18, 0.75, Form("WP: %s", wpLabel.c_str()));

    c->SaveAs(outputPath.c_str());
    delete frame;
    delete c;
}

// Same plot as DrawRJERAbsoluteSummary, with the error of each correction
// split in two:
//   err_total^2 = err^2|s_d + err^2|s_m + err^2|s_t
//   err_unc     = |dR/ds_t| * ds_t = R * (dR/R)|s_t   (from RJERRelErrFromSigmaTruth)
//   err_corr^2  = err_total^2 - err_unc^2
// The s_d and s_m parts come from the same events as the nominal, so they
// mostly cancel against it; the s_t part does not.
// Drawn as: outer bar = err_total, inner black tick = err_unc, grey band
// around the nominal line with half-width nomCorrErr (not drawn when 0).
inline void DrawRJERAbsoluteSummaryCorrAware(const std::vector<VariantRatio>& results,
                                             std::function<double(const VariantRatio&)> getValue,
                                             std::function<double(const VariantRatio&)> getError, // err_total
                                             std::function<double(const VariantRatio&)> getUncErr, // err_unc
                                             double nomValue, double nomError, double nomCorrErr,
                                             const std::string& canvasTag, const std::string& plotTitle,
                                             const std::string& yTitle, const std::string& outputPath,
                                             const std::string& wpLabel, Color_t nomColor, Color_t varColor) {
    const int nVar = (int)results.size();
    if (nVar == 0) return;
    const int n = nVar + 1;

    TCanvas* c = new TCanvas(("c_" + canvasTag).c_str(), plotTitle.c_str(), 900, 600);
    c->SetLeftMargin(0.15);
    c->SetBottomMargin(0.20);
    c->SetTopMargin(0.10);
    c->SetRightMargin(0.05);
    c->SetGridy();

    // y range from the total errors and the band
    std::vector<double> vals(n), errs(n);
    vals[0] = nomValue;
    errs[0] = std::abs(nomError);
    for (int i = 0; i < nVar; ++i) {
        vals[i + 1] = getValue(results[i]);
        errs[i + 1] = std::abs(getError(results[i]));
    }
    double ymin = 1e9, ymax = -1e9;
    for (int i = 0; i < n; ++i) {
        double lo = vals[i] - errs[i];
        double hi = vals[i] + errs[i];
        if (lo < ymin) ymin = lo;
        if (hi > ymax) ymax = hi;
    }
    double bandLo = nomValue - std::abs(nomCorrErr);
    double bandHi = nomValue + std::abs(nomCorrErr);
    if (bandLo < ymin) ymin = bandLo;
    if (bandHi > ymax) ymax = bandHi;

    double range = (ymax - ymin);
    if (range < 0.02) range = 0.05;
    double y_lo = ymin - 0.35 * range;
    double y_hi = ymax + 0.35 * range;

    TH1F* frame = new TH1F(("frame_" + canvasTag).c_str(), plotTitle.c_str(), n, -0.5, n - 0.5);
    frame->SetMinimum(y_lo);
    frame->SetMaximum(y_hi);
    frame->GetYaxis()->SetTitle(yTitle.c_str());
    frame->GetYaxis()->SetTitleSize(0.045);
    frame->GetYaxis()->SetTitleOffset(1.4);
    frame->GetYaxis()->SetLabelSize(0.040);
    frame->GetXaxis()->SetLabelSize(0.045);
    frame->GetXaxis()->SetTickLength(0);
    frame->SetStats(0);
    frame->GetXaxis()->SetBinLabel(1, "Nominal");
    for (int i = 0; i < nVar; ++i) frame->GetXaxis()->SetBinLabel(i + 2, results[i].label.c_str());
    frame->Draw("AXIS");

    // Band for the correlated part of the nominal error
    if (nomCorrErr > 0) {
        TGraphErrors* g_band = new TGraphErrors(2);
        g_band->SetPoint(0, -0.5, nomValue);
        g_band->SetPoint(1, n - 0.5, nomValue);
        g_band->SetPointError(0, 0, std::abs(nomCorrErr));
        g_band->SetPointError(1, 0, std::abs(nomCorrErr));
        g_band->SetFillColorAlpha(kGray + 1, 0.30);
        g_band->SetFillStyle(1001);
        g_band->SetLineWidth(0);
        g_band->Draw("3 same");
    }

    TLine* ref = new TLine(-0.5, nomValue, n - 0.5, nomValue);
    ref->SetLineStyle(2);
    ref->SetLineColor(kGray + 2);
    ref->SetLineWidth(2);
    ref->Draw("same");

    // Nominal point, total error
    double x_nom[1] = {0.0};
    double y_nom[1] = {nomValue};
    double ex_nom[1] = {0.0};
    double ey_nom[1] = {std::abs(nomError)};
    TGraphErrors* g_nom = new TGraphErrors(1, x_nom, y_nom, ex_nom, ey_nom);
    g_nom->SetMarkerStyle(21);
    g_nom->SetMarkerSize(1.8);
    g_nom->SetMarkerColor(nomColor);
    g_nom->SetLineColor(nomColor);
    g_nom->SetLineWidth(2);
    g_nom->Draw("P same");

    // Outer bar: total error
    std::vector<double> xs(nVar), ys(nVar), exs(nVar, 0.0), eys_total(nVar), eys_unc(nVar);
    for (int i = 0; i < nVar; ++i) {
        xs[i] = i + 1;
        ys[i] = getValue(results[i]);
        eys_total[i] = std::abs(getError(results[i]));
        eys_unc[i] = std::abs(getUncErr(results[i]));
    }
    TGraphErrors* g_var_total = new TGraphErrors(nVar, xs.data(), ys.data(), exs.data(), eys_total.data());
    g_var_total->SetMarkerStyle(20);
    g_var_total->SetMarkerSize(1.6);
    g_var_total->SetMarkerColor(varColor);
    g_var_total->SetLineColor(varColor);
    g_var_total->SetLineWidth(2);
    g_var_total->Draw("P same");

    // Inner tick: uncorrelated error, thick black with end caps so it stands
    // out from the outer bar
    TGraphErrors* g_var_unc = new TGraphErrors(nVar, xs.data(), ys.data(), exs.data(), eys_unc.data());
    g_var_unc->SetMarkerStyle(20);
    g_var_unc->SetMarkerSize(1.6);
    g_var_unc->SetMarkerColor(varColor);
    g_var_unc->SetLineColor(kBlack);
    g_var_unc->SetLineWidth(4);
    g_var_unc->Draw("[]  same");

    TLatex tx;
    tx.SetTextSize(0.028);
    tx.SetTextAlign(22);
    tx.DrawLatex(0, nomValue + std::abs(nomError) + 0.03 * range,
                 Form("%.3f #pm %.3f", nomValue, std::abs(nomError)));
    for (int i = 0; i < nVar; ++i) {
        tx.DrawLatex(i + 1, ys[i] + eys_total[i] + 0.03 * range, Form("%.3f #pm %.3f", ys[i], eys_total[i]));
    }

    TLegend* leg = new TLegend(0.62, 0.74, 0.94, 0.89);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.028);
    leg->AddEntry(g_var_total, "outer bar: total #delta R_{JER}", "le");
    leg->AddEntry(g_var_unc, "inner tick: uncorr. (#delta#sigma_{t})", "le");
    if (nomCorrErr > 0) {
        TGraphErrors* g_band_legend = new TGraphErrors(1);
        g_band_legend->SetFillColorAlpha(kGray + 1, 0.30);
        g_band_legend->SetFillStyle(1001);
        leg->AddEntry(g_band_legend, "grey band: corr. nominal", "f");
    }
    leg->Draw();

    AtlasUtils::DrawAtlasLabel(0.18, 0.85);
    TLatex wpTxC;
    wpTxC.SetNDC();
    wpTxC.SetTextFont(42);
    wpTxC.SetTextSize(0.035);
    wpTxC.DrawLatex(0.18, 0.75, Form("WP: %s", wpLabel.c_str()));

    c->SaveAs(outputPath.c_str());
    delete frame;
    delete c;
}

// One double ratio per correction, around 1. The accessor picks which
// DoubleRatio of VariantRatio is drawn.
inline void DrawGenericRatioSummary(const std::vector<VariantRatio>& results,
                                    std::function<DoubleRatio(const VariantRatio&)> getDR,
                                    const std::string& canvasTag, // canvas and frame names
                                    const std::string& plotTitle, const std::string& yTitle,
                                    const std::string& outputPath, const std::string& wpLabel,
                                    Color_t markerColor) {
    const int n = (int)results.size();
    if (n == 0) return;

    TCanvas* c = new TCanvas(("c_" + canvasTag).c_str(), plotTitle.c_str(), 900, 600);
    c->SetLeftMargin(0.15);
    c->SetBottomMargin(0.20);
    c->SetTopMargin(0.10);
    c->SetRightMargin(0.05);
    c->SetGridy();

    // y range always includes 0.95-1.05 so the line at 1 is visible
    double ymin = 1e9, ymax = -1e9;
    for (const auto& r : results) {
        DoubleRatio d = getDR(r);
        double lo = d.value - std::abs(d.err);
        double hi = d.value + std::abs(d.err);
        if (lo < ymin) ymin = lo;
        if (hi > ymax) ymax = hi;
    }
    double range = (ymax - ymin);
    if (range < 0.05) range = 0.1;
    double y_lo = ymin - 0.25 * range;
    double y_hi = ymax + 0.25 * range;
    if (y_lo > 0.95) y_lo = 0.95;
    if (y_hi < 1.05) y_hi = 1.05;

    TH1F* frame = new TH1F(("frame_" + canvasTag).c_str(), plotTitle.c_str(), n, -0.5, n - 0.5);
    frame->SetMinimum(y_lo);
    frame->SetMaximum(y_hi);
    frame->GetYaxis()->SetTitle(yTitle.c_str());
    frame->GetYaxis()->SetTitleSize(0.045);
    frame->GetYaxis()->SetTitleOffset(1.4);
    frame->GetYaxis()->SetLabelSize(0.040);
    frame->GetXaxis()->SetLabelSize(0.045);
    frame->GetXaxis()->SetTickLength(0);
    frame->SetStats(0);
    for (int i = 0; i < n; ++i) frame->GetXaxis()->SetBinLabel(i + 1, results[i].label.c_str());
    frame->Draw("AXIS");

    TLine* ref = new TLine(-0.5, 1.0, n - 0.5, 1.0);
    ref->SetLineStyle(2);
    ref->SetLineColor(kBlack);
    ref->SetLineWidth(2);
    ref->Draw("same");

    std::vector<double> xs(n), ys(n), exs(n, 0.0), eys(n);
    for (int i = 0; i < n; ++i) {
        DoubleRatio d = getDR(results[i]);
        xs[i] = i;
        ys[i] = d.value;
        eys[i] = std::abs(d.err);
    }

    TGraphErrors* g = new TGraphErrors(n, xs.data(), ys.data(), exs.data(), eys.data());
    g->SetMarkerStyle(20);
    g->SetMarkerSize(1.6);
    g->SetMarkerColor(markerColor);
    g->SetLineColor(markerColor);
    g->SetLineWidth(2);
    g->Draw("P same");

    TLatex tx;
    tx.SetTextSize(0.028);
    tx.SetTextAlign(22);
    for (int i = 0; i < n; ++i) {
        DoubleRatio d = getDR(results[i]);
        tx.DrawLatex(i, ys[i] + eys[i] + 0.03 * range, Form("%.3f #pm %.3f", d.value, d.err));
    }

    AtlasUtils::DrawAtlasLabel(0.18, 0.85);
    TLatex wpTx;
    wpTx.SetNDC();
    wpTx.SetTextFont(42);
    wpTx.SetTextSize(0.035);
    wpTx.DrawLatex(0.18, 0.75, Form("WP: %s", wpLabel.c_str()));

    c->SaveAs(outputPath.c_str());
    delete frame;
    delete c;
}

// Uncorrelated and corr-aware errors side by side: each correction is drawn
// twice (x - 0.10 and x + 0.10) with the same central value, to show how
// much the correlation treatment changes the error.
inline void DrawResolutionErrorComparison(const std::vector<VariantRatio>& results,
                                          std::function<DoubleRatio(const VariantRatio&)> getUncorr,
                                          std::function<DoubleRatio(const VariantRatio&)> getCorr,
                                          const std::string& canvasTag, const std::string& plotTitle,
                                          const std::string& yTitle, const std::string& outputPath,
                                          const std::string& wpLabel) {
    const int n = (int)results.size();
    if (n == 0) return;

    TCanvas* c = new TCanvas(("c_" + canvasTag).c_str(), plotTitle.c_str(), 950, 600);
    c->SetLeftMargin(0.15);
    c->SetBottomMargin(0.20);
    c->SetTopMargin(0.10);
    c->SetRightMargin(0.05);
    c->SetGridy();

    // y range from the uncorrelated errors, normally the larger ones
    double ymin = 1e9, ymax = -1e9;
    for (const auto& r : results) {
        DoubleRatio d = getUncorr(r);
        double lo = d.value - std::abs(d.err);
        double hi = d.value + std::abs(d.err);
        if (lo < ymin) ymin = lo;
        if (hi > ymax) ymax = hi;
    }
    double range = (ymax - ymin);
    if (range < 0.05) range = 0.1;
    double y_lo = ymin - 0.30 * range;
    double y_hi = ymax + 0.30 * range;
    if (y_lo > 0.95) y_lo = 0.95;
    if (y_hi < 1.05) y_hi = 1.05;

    TH1F* frame = new TH1F(("frame_" + canvasTag).c_str(), plotTitle.c_str(), n, -0.5, n - 0.5);
    frame->SetMinimum(y_lo);
    frame->SetMaximum(y_hi);
    frame->GetYaxis()->SetTitle(yTitle.c_str());
    frame->GetYaxis()->SetTitleSize(0.045);
    frame->GetYaxis()->SetTitleOffset(1.4);
    frame->GetYaxis()->SetLabelSize(0.040);
    frame->GetXaxis()->SetLabelSize(0.045);
    frame->GetXaxis()->SetTickLength(0);
    frame->SetStats(0);
    for (int i = 0; i < n; ++i) frame->GetXaxis()->SetBinLabel(i + 1, results[i].label.c_str());
    frame->Draw("AXIS");

    TLine* ref = new TLine(-0.5, 1.0, n - 0.5, 1.0);
    ref->SetLineStyle(2);
    ref->SetLineColor(kBlack);
    ref->SetLineWidth(2);
    ref->Draw("same");

    // Uncorrelated: open grey square at x - 0.10
    std::vector<double> xs_u(n), ys_u(n), exs_u(n, 0.0), eys_u(n);
    for (int i = 0; i < n; ++i) {
        DoubleRatio d = getUncorr(results[i]);
        xs_u[i] = i - 0.10;
        ys_u[i] = d.value;
        eys_u[i] = std::abs(d.err);
    }
    TGraphErrors* g_uncorr = new TGraphErrors(n, xs_u.data(), ys_u.data(), exs_u.data(), eys_u.data());
    g_uncorr->SetMarkerStyle(25);
    g_uncorr->SetMarkerSize(1.5);
    g_uncorr->SetMarkerColor(kGray + 2);
    g_uncorr->SetLineColor(kGray + 2);
    g_uncorr->SetLineWidth(2);
    g_uncorr->SetLineStyle(2);
    g_uncorr->Draw("P same");

    // Corr-aware: filled circle at x + 0.10
    std::vector<double> xs_c(n), ys_c(n), exs_c(n, 0.0), eys_c(n);
    for (int i = 0; i < n; ++i) {
        DoubleRatio d = getCorr(results[i]);
        xs_c[i] = i + 0.10;
        ys_c[i] = d.value;
        eys_c[i] = std::abs(d.err);
    }
    TGraphErrors* g_corr = new TGraphErrors(n, xs_c.data(), ys_c.data(), exs_c.data(), eys_c.data());
    g_corr->SetMarkerStyle(20);
    g_corr->SetMarkerSize(1.6);
    g_corr->SetMarkerColor(kAzure + 2);
    g_corr->SetLineColor(kAzure + 2);
    g_corr->SetLineWidth(2);
    g_corr->Draw("P same");

    TLegend* leg = new TLegend(0.60, 0.78, 0.93, 0.90);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.030);
    leg->AddEntry(g_uncorr, "Uncorrelated (conservative)", "lp");
    leg->AddEntry(g_corr, "Correlation-aware (truth-only)", "lp");
    leg->Draw();

    AtlasUtils::DrawAtlasLabel(0.18, 0.85);
    TLatex wpTx;
    wpTx.SetNDC();
    wpTx.SetTextFont(42);
    wpTx.SetTextSize(0.035);
    wpTx.DrawLatex(0.18, 0.75, Form("WP: %s", wpLabel.c_str()));

    c->SaveAs(outputPath.c_str());
    delete frame;
    delete c;
}

// Readable table of all double ratios and fit inputs. The nominal
// block is taken from results[0], so it is missing if that correction failed.
inline void WriteSummaryTable(const std::vector<VariantRatio>& results, const std::string& outputPath,
                              const std::string& wpLabel, const VariantConfig& cfg_nom) {
    std::ofstream out(outputPath);
    if (!out.is_open()) {
        Log::error("Ratio_Response_Inclusive", "Cannot write: " + outputPath);
        return;
    }

    out << "========================================================================================\n";
    out << "  Inclusive JES / JER Double Ratios vs Nominal — WP: " << wpLabel << "\n";
    out << "========================================================================================\n\n";
    out << "Definitions:\n";
    out << "  R_JES     = <" + cfg_nom.branch_R + ">_data / <" + cfg_nom.branch_R + ">_combMC\n";
    out << "  sigma_insitu    = sqrt(sigma^2_data - sigma^2_truth_ref)       (data-side)\n";
    out << "  sigma^MC_insitu = sqrt(sigma^2_combMC - sigma^2_truth_ref)\n";
    out << "  R_JER     = sigma_insitu / sigma^MC_insitu                    (data/MC JER)\n";
    out << "  DR(X)     = X_variant / X_nominal\n\n";
    out << "Error treatments for resolution double ratios:\n";
    out << "  uncorr   : all fit errors propagated as independent (conservative).\n";
    out << "  corr-aw  : sigma_data/sigma_MC fluctuations largely shared between variants;\n";
    out << "             only sigma_truth fit-statistical noise is uncorrelated.\n\n";

    out << std::fixed << std::setprecision(5);
    out << std::left << std::setw(18) << "Variant" << std::right << std::setw(22) << "DR(R_JES) unc"
        << std::setw(22) << "DR(R_JES) cor" << std::setw(22) << "DR(R_JES) full" << std::setw(24)
        << "DR(sigma_insitu) unc" << std::setw(24) << "DR(sigma_insitu) cor" << std::setw(22)
        << "DR(R_JER) unc" << std::setw(22) << "DR(R_JER) cor"
        << "\n";
    out << std::string(178, '-') << "\n";

    for (const auto& r : results) {
        if (!r.valid) {
            out << std::left << std::setw(18) << r.label << "   [incomplete fit bundle — skipped]\n";
            continue;
        }
        out << std::left << std::setw(18) << r.label << std::right << std::setw(22)
            << Form("%.4f +/- %.4f", r.jes.value, r.jes.err) << std::setw(22)
            << Form("%.4f +/- %.4f", r.jes_corr.value, r.jes_corr.err) << std::setw(22)
            << Form("%.4f +/- %.4f", r.jes_full.value, r.jes_full.err) << std::setw(24)
            << Form("%.4f +/- %.4f", r.sig_insitu.value, r.sig_insitu.err) << std::setw(24)
            << Form("%.4f +/- %.4f", r.sig_insitu_corr.value, r.sig_insitu_corr.err) << std::setw(22)
            << Form("%.4f +/- %.4f", r.rjer.value, r.rjer.err) << std::setw(22)
            << Form("%.4f +/- %.4f", r.rjer_corr.value, r.rjer_corr.err) << "\n";
    }

    out << "\n--- Raw fit inputs per variant ---\n";
    out << std::left << std::setw(18) << "Variant" << std::right << std::setw(12) << "mu_d" << std::setw(12)
        << "mu_c" << std::setw(12) << "sig_d" << std::setw(12) << "sig_c" << std::setw(12) << "sig_t"
        << std::setw(14) << "sig_insitu" << std::setw(20) << "R_JER (data/MC)"
        << "\n";
    out << std::string(110, '-') << "\n";
    for (const auto& r : results) {
        if (!r.valid) continue;
        out << std::left << std::setw(18) << r.label << std::right << std::setw(12) << r.mean_data_var
            << std::setw(12) << r.mean_comb_var << std::setw(12) << r.sigma_data_var << std::setw(12)
            << r.sigma_comb_var << std::setw(12) << r.sigma_truth_var << std::setw(14) << r.jer_var_val
            << std::setw(20) << Form("%.4f +/- %.4f", r.R_JER_var, r.R_JER_var_err) << "\n";
    }

    out << "\nNominal reference values:\n";
    if (!results.empty() && results[0].valid) {
        out << Form(("  <" + cfg_nom.branch_R + ">_data_nom    = %.5f\n").c_str(), results[0].mean_data_nom);
        out << Form(("  <" + cfg_nom.branch_R + ">_comb_nom    = %.5f\n").c_str(), results[0].mean_comb_nom);
        out << Form("  sigma_data_nom     = %.5f\n", results[0].sigma_data_nom);
        out << Form("  sigma_comb_nom     = %.5f\n", results[0].sigma_comb_nom);
        out << Form("  sigma_truth_nom    = %.5f\n", results[0].sigma_truth_nom);
        out << Form("  sigma_insitu_nom   = %.5f +/- %.5f\n", results[0].jer_nom_val, results[0].jer_nom_err);
        out << Form("  R_JER_nom          = %.5f +/- %.5f\n", results[0].R_JER_nom, results[0].R_JER_nom_err);
    }
    out.close();
    Log::info("Ratio_Response_Inclusive", "Summary table: " + outputPath);
}

// Defined after ProcessWP, which calls it.
inline void DrawDualAxisInclusiveSummary(
    const std::vector<VariantRatio>& results, std::function<double(const VariantRatio&)> getSingleVar,
    std::function<double(const VariantRatio&)> getSingleVarErr,
    std::function<double(const VariantRatio&)> getSingleNom,
    std::function<double(const VariantRatio&)> getSingleNomErr,
    std::function<double(const VariantRatio&)> getDR, std::function<double(const VariantRatio&)> getDRErr,
    const std::string& canvasName, const std::string& leftTitle, const std::string& rightTitle,
    const std::string& outputPath, const std::string& wpLabel, const std::vector<size_t>& selected = {});

// CSV with one row per correction, in the same 44-column layout as the
// pT-binned Summary_Ratios CSV of Response_Ratios_pTbins_Mcuts.h, so the same
// reader (RatioResponseCombined::ReadCSV in Response_Ratios_Combined.h) can be
// used by Response_Summary_Combined_Inclusive.cxx. Here mass_bin = "inclusive",
// pt_bin_index = 0, and pt_lo/pt_hi are the channel's full range (DB 200-450,
// MJB 450-1200 GeV).
inline void WriteSummaryInclusiveCSV(const std::vector<VariantRatio>& results, const std::string& csvPath,
                                     const std::string& channel, const std::string& wp_file) {
    std::ofstream out(csvPath);
    if (!out.is_open()) {
        Log::error("Ratio_Response_Inclusive", "Cannot write: " + csvPath);
        return;
    }
    out << "channel,wp,variant,mass_bin,pt_bin_index,pt_lo,pt_hi,pt_center,"
        << "rjes_dr,rjes_dr_err,"
        << "sig_insitu_dr,sig_insitu_dr_err_unc,sig_insitu_dr_err_corr,"
        << "rjer_dr,rjer_dr_err_unc,rjer_dr_err_corr,"
        << "rjes_nom,rjes_nom_err,rjes_var,rjes_var_err,"
        << "rjer_nom,rjer_nom_err,rjer_var,rjer_var_err,"
        << "bjes_data_nom,bjes_data_nom_err,bjes_mc_nom,bjes_mc_nom_err,"
        << "bjes_data_var,bjes_data_var_err,bjes_mc_var,bjes_mc_var_err,"
        << "bjer_data_nom,bjer_data_nom_err,bjer_mc_nom,bjer_mc_nom_err,"
        << "bjer_data_var,bjer_data_var_err,bjer_mc_var,bjer_mc_var_err,"
        // Columns 40-43 (counting from 0): extra DR errors. 
        << "rjes_dr_err_corr,rjes_dr_err_full,sig_insitu_dr_err_full,rjer_dr_err_full\n";

    const bool isMJB = (channel == "MJB");
    const double lo = isMJB ? 450.0 : 200.0;
    const double hi = isMJB ? 1200.0 : 450.0;
    const double mid = 0.5 * (lo + hi);

    out << std::setprecision(8);
    for (const auto& r : results) {
        if (!r.valid) continue;
        out << channel << "," << wp_file << "," << r.key << ",inclusive,0," << lo << "," << hi << "," << mid
            << "," << r.jes.value << "," << r.jes.err << "," << r.sig_insitu.value << "," << r.sig_insitu.err
            << "," << r.sig_insitu_corr.err << "," << r.rjer.value << "," << r.rjer.err << ","
            << r.rjer_corr.err << "," << r.R_JES_nom << "," << r.R_JES_nom_err << "," << r.R_JES_var << ","
            << r.R_JES_var_err << "," << r.R_JER_nom << "," << r.R_JER_nom_err << "," << r.R_JER_var << ","
            << r.R_JER_var_err << "," << r.bjes_data_nom << "," << r.bjes_data_nom_err << "," << r.bjes_mc_nom
            << "," << r.bjes_mc_nom_err << "," << r.bjes_data_var << "," << r.bjes_data_var_err << ","
            << r.bjes_mc_var << "," << r.bjes_mc_var_err << "," << r.bjer_data_nom << ","
            << r.bjer_data_nom_err << "," << r.bjer_mc_nom << "," << r.bjer_mc_nom_err << ","
            << r.bjer_data_var << "," << r.bjer_data_var_err << "," << r.bjer_mc_var << ","
            << r.bjer_mc_var_err << "," << r.jes_corr.err << "," << r.jes_full.err << ","
            << r.sig_insitu_full.err << "," << r.rjer_full.err << "\n";
    }
    out.close();
    Log::info("Ratio_Response_Inclusive", "CSV summary saved: " + csvPath);
}

// One WP: read the nominal once, then build all ratios for each correction
// and write the table, CSV and plots.
inline void ProcessWP(const std::vector<std::string>& variantKeys, const std::string& wp_raw,
                      const std::string& channel) {
    std::string wp_file = PathBuilder::WPtoFileWP(wp_raw);
    std::string wp_label = PathBuilder::WPtoLabel(wp_raw);

    std::string outDir =
        (RunConfig::PlotDir() + "/Response/Ratios/") + channel + "_Inclusive_Ratio/" + wp_file + "/";
    PathBuilder::MakeDir(outDir);

    Log::info("Ratio_Response_Inclusive", Form("=== WP: %s  (%s)  — variants: %d ===", wp_file.c_str(),
                                               wp_label.c_str(), (int)variantKeys.size()));

    // Nominal
    VariantConfig cfg_nom = GetVariant(channel);
    VariantFitBundle bnom = ReadVariantBundle(cfg_nom, wp_file);
    if (!bnom.valid) {
        Log::error("Ratio_Response_Inclusive", "Nominal fit bundle missing — aborting WP " + wp_file);
        return;
    }

    // sigma_insitu, data only, for the DR(sigma_insitu) cross-check
    double jer_nom = JERFromEq92(bnom.data.sigma, bnom.truth_ref.sigma);
    double jer_nom_err =
        JERErrFromEq92(bnom.data.sigma, bnom.data.sigma_err, bnom.truth_ref.sigma, bnom.truth_ref.sigma_err);

    // R_JER data/MC. Inside one variant sigma_truth enters both sigma_insitu
    // and sigma^MC_insitu, so ComputeJERRatio propagates it as one variable.
    FitUtils::JERRatio rjer_nom =
        FitUtils::ComputeJERRatio(bnom.data.sigma, bnom.data.sigma_err, bnom.comb.sigma, bnom.comb.sigma_err,
                                  bnom.truth_ref.sigma, bnom.truth_ref.sigma_err);

    Log::info("Ratio_Response_Inclusive",
              Form(("Nominal: <" + cfg_nom.branch_R + ">_d=%.4f <" + cfg_nom.branch_R
                    + ">_c=%.4f sig_d=%.4f sig_c=%.4f sig_t=%.4f  "
                      "sig_insitu=%.4f+/-%.4f  R_JER=%.4f+/-%.4f%s")
                       .c_str(),
                   bnom.data.mean, bnom.comb.mean, bnom.data.sigma, bnom.comb.sigma, bnom.truth_ref.sigma,
                   jer_nom, jer_nom_err, rjer_nom.ratio, rjer_nom.ratioErr,
                   rjer_nom.valid ? "" : " [INVALID — negative variance]"));

    // Corrections
    std::vector<VariantRatio> results;
    for (const auto& key : variantKeys) {
        VariantConfig cfg_var = GetVariant(key);
        VariantFitBundle bvar = ReadVariantBundle(cfg_var, wp_file);

        VariantRatio r;
        r.key = key;
        r.label = cfg_var.label;
        if (!bvar.valid) {
            results.push_back(r);
            continue;
        }

        double jer_var = JERFromEq92(bvar.data.sigma, bvar.truth_ref.sigma);
        double jer_var_err = JERErrFromEq92(bvar.data.sigma, bvar.data.sigma_err, bvar.truth_ref.sigma,
                                            bvar.truth_ref.sigma_err);

        FitUtils::JERRatio rjer_var =
            FitUtils::ComputeJERRatio(bvar.data.sigma, bvar.data.sigma_err, bvar.comb.sigma,
                                      bvar.comb.sigma_err, bvar.truth_ref.sigma, bvar.truth_ref.sigma_err);

        // DR(R_JES) with the three error treatments. 
        r.jes = ComputeJESDoubleRatio(bvar.data.mean, bvar.data.mean_err, bvar.comb.mean, bvar.comb.mean_err,
                                      bnom.data.mean, bnom.data.mean_err, bnom.comb.mean, bnom.comb.mean_err);
        r.jes_corr = ComputeJESDoubleRatioCorrAware(bvar.data.mean, bvar.data.mean_err, bvar.comb.mean,
                                                    bvar.comb.mean_err, bnom.data.mean, bnom.data.mean_err,
                                                    bnom.comb.mean, bnom.comb.mean_err);
        {
            // rho measured by Bootstrap_Correlations for this correction
            // ("<correction>:inclusive"); default matrix if not in the file.
            CorrelationConfig::CorrelationMatrix Cresp = CorrelationConfig::CorrFor(
                CorrelationConfig::CorrDomain::Response, channel, wp_file, "inclusive", cfg_var.path_tag);
            r.jes_full = ComputeJESDoubleRatioFullCorr(bvar.data.mean, bvar.data.mean_err, bvar.comb.mean,
                                                       bvar.comb.mean_err, bnom.data.mean, bnom.data.mean_err,
                                                       bnom.comb.mean, bnom.comb.mean_err, Cresp);
        }

        // DR(sigma_insitu): same central value, two errors
        r.sig_insitu = ComputeDoubleRatio(jer_var, jer_var_err, jer_nom, jer_nom_err);
        r.sig_insitu_corr.value = r.sig_insitu.value;
        r.sig_insitu_corr.err =
            SigmaInsituDoubleRatioErrCorrAware(jer_var, bvar.truth_ref.sigma, bvar.truth_ref.sigma_err,
                                               jer_nom, bnom.truth_ref.sigma, bnom.truth_ref.sigma_err);

        // Full errors of the resolution double ratios: c^T Sigma c over the six
        // widths (data, combined MC, truth; correction and nominal), with the
        // bootstrap correlations of this correction (default matrix if none).
        const CorrelationConfig::CorrelationMatrix Cres = CorrelationConfig::CorrFor(
            CorrelationConfig::CorrDomain::Response, channel, wp_file, "inclusive", cfg_var.path_tag);
        r.sig_insitu_full.value = r.sig_insitu.value;
        r.sig_insitu_full.err =
            std::abs(r.sig_insitu.value)
            * CorrelationConfig::ResolutionDRRelErrFull(
                false, bvar.data.sigma, bvar.data.sigma_err, bnom.data.sigma, bnom.data.sigma_err, bvar.comb.sigma,
                bvar.comb.sigma_err, bnom.comb.sigma, bnom.comb.sigma_err, bvar.truth_ref.sigma,
                bvar.truth_ref.sigma_err, bnom.truth_ref.sigma, bnom.truth_ref.sigma_err, Cres);

        // DR(R_JER). The uncorrelated error combines the two R_JER errors from
        // ComputeJERRatio as if correction and nominal were independent.
        if (rjer_var.valid && rjer_nom.valid) {
            r.rjer = ComputeDoubleRatio(rjer_var.ratio, rjer_var.ratioErr, rjer_nom.ratio, rjer_nom.ratioErr);
            r.rjer_corr.value = r.rjer.value;
            r.rjer_corr.err = RJERDoubleRatioErrCorrAware(rjer_var.ratio, bvar.data.sigma, bvar.comb.sigma,
                                                          bvar.truth_ref.sigma, bvar.truth_ref.sigma_err,
                                                          rjer_nom.ratio, bnom.data.sigma, bnom.comb.sigma,
                                                          bnom.truth_ref.sigma, bnom.truth_ref.sigma_err);
            r.rjer_full.value = r.rjer.value;
            r.rjer_full.err =
                std::abs(r.rjer.value)
                * CorrelationConfig::ResolutionDRRelErrFull(
                    true, bvar.data.sigma, bvar.data.sigma_err, bnom.data.sigma, bnom.data.sigma_err,
                    bvar.comb.sigma, bvar.comb.sigma_err, bnom.comb.sigma, bnom.comb.sigma_err,
                    bvar.truth_ref.sigma, bvar.truth_ref.sigma_err, bnom.truth_ref.sigma,
                    bnom.truth_ref.sigma_err, Cres);
        } else {
            Log::warning("Ratio_Response_Inclusive",
                         "R_JER invalid (negative variance) for " + cfg_var.key + " — leaving DR(R_JER) = 0");
        }

        r.sigma_data_var = bvar.data.sigma;
        r.sigma_data_nom = bnom.data.sigma;
        r.sigma_comb_var = bvar.comb.sigma;
        r.sigma_comb_nom = bnom.comb.sigma;
        r.sigma_truth_var = bvar.truth_ref.sigma;
        r.sigma_truth_nom = bnom.truth_ref.sigma;
        r.mean_data_var = bvar.data.mean;
        r.mean_data_nom = bnom.data.mean;
        r.mean_comb_var = bvar.comb.mean;
        r.mean_comb_nom = bnom.comb.mean;

        // b-JES / b-JER inputs for the CSV
        r.bjes_data_var = bvar.data.mean;
        r.bjes_data_var_err = bvar.data.mean_err;
        r.bjes_mc_var = bvar.comb.mean;
        r.bjes_mc_var_err = bvar.comb.mean_err;
        r.bjes_data_nom = bnom.data.mean;
        r.bjes_data_nom_err = bnom.data.mean_err;
        r.bjes_mc_nom = bnom.comb.mean;
        r.bjes_mc_nom_err = bnom.comb.mean_err;
        r.bjer_data_var = rjer_var.sig_insitu;
        r.bjer_data_var_err = rjer_var.sig_insitu_err;
        r.bjer_mc_var = rjer_var.sig_MC_insitu;
        r.bjer_mc_var_err = rjer_var.sig_MC_insitu_err;
        r.bjer_data_nom = rjer_nom.sig_insitu;
        r.bjer_data_nom_err = rjer_nom.sig_insitu_err;
        r.bjer_mc_nom = rjer_nom.sig_MC_insitu;
        r.bjer_mc_nom_err = rjer_nom.sig_MC_insitu_err;

        r.jer_var_val = jer_var;
        r.jer_nom_val = jer_nom;
        r.jer_var_err = jer_var_err;
        r.jer_nom_err = jer_nom_err;
        r.R_JER_var = rjer_var.ratio;
        r.R_JER_var_err = rjer_var.ratioErr;
        r.R_JER_nom = rjer_nom.ratio;
        r.R_JER_nom_err = rjer_nom.ratioErr;

        // Split the R_JER errors into the sigma_truth part and the rest, for the
        // CorrAware absolute plot. Only possible when ComputeJERRatio succeeded.
        if (rjer_var.valid && r.R_JER_var > 0) {
            double rel_unc_var = RJERRelErrFromSigmaTruth(bvar.data.sigma, bvar.comb.sigma,
                                                          bvar.truth_ref.sigma, bvar.truth_ref.sigma_err);
            r.R_JER_var_err_unc = std::abs(r.R_JER_var) * rel_unc_var;
            // err_corr^2 = err_total^2 - err_unc^2
            double diff = r.R_JER_var_err * r.R_JER_var_err - r.R_JER_var_err_unc * r.R_JER_var_err_unc;
            r.R_JER_var_err_corr = (diff > 0) ? std::sqrt(diff) : 0.0;
        }
        if (rjer_nom.valid && r.R_JER_nom > 0) {
            double rel_unc_nom = RJERRelErrFromSigmaTruth(bnom.data.sigma, bnom.comb.sigma,
                                                          bnom.truth_ref.sigma, bnom.truth_ref.sigma_err);
            r.R_JER_nom_err_unc = std::abs(r.R_JER_nom) * rel_unc_nom;
            double diff = r.R_JER_nom_err * r.R_JER_nom_err - r.R_JER_nom_err_unc * r.R_JER_nom_err_unc;
            r.R_JER_nom_err_corr = (diff > 0) ? std::sqrt(diff) : 0.0;
        }

        // Absolute R_JES of correction and nominal. Data and MC are different
        // samples, so the whole error is "uncorrelated" and err_corr = 0.
        FitUtils::JESRatio rjes_var =
            FitUtils::ComputeJESRatio(bvar.data.mean, bvar.data.mean_err, bvar.comb.mean, bvar.comb.mean_err);
        FitUtils::JESRatio rjes_nom =
            FitUtils::ComputeJESRatio(bnom.data.mean, bnom.data.mean_err, bnom.comb.mean, bnom.comb.mean_err);
        if (rjes_var.valid) {
            r.R_JES_var = rjes_var.ratio;
            r.R_JES_var_err = rjes_var.ratioErr;
            r.R_JES_var_err_unc = rjes_var.ratioErr;
            r.R_JES_var_err_corr = 0.0;
        }
        if (rjes_nom.valid) {
            r.R_JES_nom = rjes_nom.ratio;
            r.R_JES_nom_err = rjes_nom.ratioErr;
            r.R_JES_nom_err_unc = rjes_nom.ratioErr;
            r.R_JES_nom_err_corr = 0.0;
        }
        r.valid = true;

        Log::info("Ratio_Response_Inclusive",
                  Form("  %-14s  R_JES=%.4f+/-%.4f  sig_insitu=%.4f+/-%.4f (cor %.4f)  "
                       "R_JER=%.4f+/-%.4f (cor %.4f)",
                       r.label.c_str(), r.jes.value, r.jes.err, r.sig_insitu.value, r.sig_insitu.err,
                       r.sig_insitu_corr.err, r.rjer.value, r.rjer.err, r.rjer_corr.err));

        results.push_back(r);
    }

    // Outputs
    WriteSummaryTable(results, outDir + "Summary_Inclusive_Ratios_" + wp_file + ".txt", wp_label, cfg_nom);
    WriteSummaryInclusiveCSV(results, outDir + "Summary_Inclusive_Ratios_" + wp_file + ".csv", channel,
                             wp_file);

    // DR(R_JES), uncorrelated errors only. The corr-aware R_JES error is shown
    // in the Summary/ plots below and both are in the table and CSV.
    DrawGenericRatioSummary(
        results, [](const VariantRatio& v) { return v.jes; }, "RJES_summary_" + wp_file,
        "Inclusive R_{JES} Corrected/Nominal", "R_{JES}^{Corrected} / R_{JES}^{Nominal}",
        outDir + "R_JES_Inclusive_DoubleRatio_" + wp_file + ".pdf", wp_label, kRed + 1);

    // DR(sigma_insitu): uncorrelated, corr-aware, and both side by side
    DrawGenericRatioSummary(
        results, [](const VariantRatio& v) { return v.sig_insitu; }, "sig_insitu_uncorr_" + wp_file,
        "Inclusive #sigma_{insitu} Corrected/Nominal (uncorrelated errors)",
        "#sigma_{insitu}^{Corrected} / #sigma_{insitu}^{Nominal}",
        outDir + "sigma_insitu_Inclusive_DoubleRatio_" + wp_file + ".pdf", wp_label, kGreen + 3);

    DrawGenericRatioSummary(
        results, [](const VariantRatio& v) { return v.sig_insitu_corr; }, "sig_insitu_corr_" + wp_file,
        "Inclusive #sigma_{insitu} Corrected/Nominal (correlation-aware errors)",
        "#sigma_{insitu}^{Corrected} / #sigma_{insitu}^{Nominal}",
        outDir + "sigma_insitu_Inclusive_DoubleRatio_" + wp_file + "_CorrAware.pdf", wp_label, kAzure + 2);

    DrawResolutionErrorComparison(
        results, [](const VariantRatio& v) { return v.sig_insitu; },
        [](const VariantRatio& v) { return v.sig_insitu_corr; }, "sig_insitu_cmp_" + wp_file,
        "#sigma_{insitu} Corrected/Nominal — error-treatment comparison",
        "#sigma_{insitu}^{Corrected} / #sigma_{insitu}^{Nominal}",
        outDir + "sigma_insitu_Inclusive_DoubleRatio_" + wp_file + "_ErrorComparison.pdf", wp_label);

    // DR(R_JER): uncorrelated, corr-aware, and both side by side
    DrawGenericRatioSummary(
        results, [](const VariantRatio& v) { return v.rjer; }, "RJER_uncorr_" + wp_file,
        "Inclusive R_{JER} Corrected/Nominal (uncorrelated errors)",
        "R_{JER}^{Corrected} / R_{JER}^{Nominal}", outDir + "R_JER_Inclusive_DoubleRatio_" + wp_file + ".pdf",
        wp_label, kGreen + 3);

    DrawGenericRatioSummary(
        results, [](const VariantRatio& v) { return v.rjer_corr; }, "RJER_corr_" + wp_file,
        "Inclusive R_{JER} Corrected/Nominal (correlation-aware errors)",
        "R_{JER}^{Corrected} / R_{JER}^{Nominal}",
        outDir + "R_JER_Inclusive_DoubleRatio_" + wp_file + "_CorrAware.pdf", wp_label, kAzure + 2);

    DrawResolutionErrorComparison(
        results, [](const VariantRatio& v) { return v.rjer; },
        [](const VariantRatio& v) { return v.rjer_corr; }, "RJER_cmp_" + wp_file,
        "R_{JER} Corrected/Nominal — error-treatment comparison", "R_{JER}^{Corrected} / R_{JER}^{Nominal}",
        outDir + "R_JER_Inclusive_DoubleRatio_" + wp_file + "_ErrorComparison.pdf", wp_label);

    // Absolute R_JER and R_JES per correction, next to the nominal. The
    // nominal values are copied in every VariantRatio; results[0] is used.
    if (!results.empty()) {
        DrawRJERAbsoluteSummary(
            results, [](const VariantRatio& v) { return v.R_JER_var; },
            [](const VariantRatio& v) { return v.R_JER_var_err; }, results[0].R_JER_nom,
            results[0].R_JER_nom_err, "RJER_abs_" + wp_file,
            "Inclusive R_{JER} = #sigma_{insitu} / #sigma^{MC}_{insitu} per variant", "R_{JER}",
            outDir + "R_JER_Inclusive_Absolute_" + wp_file + ".pdf", wp_label, kBlack, kRed + 1);

        DrawRJERAbsoluteSummaryCorrAware(
            results, [](const VariantRatio& v) { return v.R_JER_var; },
            [](const VariantRatio& v) { return v.R_JER_var_err; },
            [](const VariantRatio& v) { return v.R_JER_var_err_unc; }, results[0].R_JER_nom,
            results[0].R_JER_nom_err, results[0].R_JER_nom_err_corr, "RJER_abs_corr_" + wp_file,
            "Inclusive R_{JER} (correlation-aware errors)", "R_{JER}",
            outDir + "R_JER_Inclusive_Absolute_" + wp_file + "_CorrAware.pdf", wp_label, kBlack, kRed + 1);

        DrawRJERAbsoluteSummary(
            results, [](const VariantRatio& v) { return v.R_JES_var; },
            [](const VariantRatio& v) { return v.R_JES_var_err; }, results[0].R_JES_nom,
            results[0].R_JES_nom_err, "RJES_abs_" + wp_file,
            "Inclusive R_{JES} = <R>_{data} / <R>_{comb. MC} per variant", "R_{JES}",
            outDir + "R_JES_Inclusive_Absolute_" + wp_file + ".pdf", wp_label, kBlack, kRed + 1);

        // For R_JES err_corr = 0, so this CorrAware plot has no band and the
        // inner tick equals the outer bar. It is kept to match the R_JER set.
        DrawRJERAbsoluteSummaryCorrAware(
            results, [](const VariantRatio& v) { return v.R_JES_var; },
            [](const VariantRatio& v) { return v.R_JES_var_err; },
            [](const VariantRatio& v) { return v.R_JES_var_err_unc; }, results[0].R_JES_nom,
            results[0].R_JES_nom_err, results[0].R_JES_nom_err_corr, "RJES_abs_corr_" + wp_file,
            "Inclusive R_{JES} (correlation-aware errors)", "R_{JES}",
            outDir + "R_JES_Inclusive_Absolute_" + wp_file + "_CorrAware.pdf", wp_label, kBlack, kRed + 1);
    }

    // Summary/ plots: data/MC ratio on the left axis and the corr-aware DR on
    // the right axis, for R_JES (Scale) and R_JER (Resolution). Made per
    // correction, for all corrections, and for the EASYBJETS_SUMMARY_VARIANTS
    // subset if set.
    std::string summaryDir = outDir + "Summary/";
    PathBuilder::MakeDir(summaryDir);

    auto sJESvar = [](const VariantRatio& v) {
        return v.R_JES_var;
    };
    auto sJESvarE = [](const VariantRatio& v) {
        return v.R_JES_var_err;
    };
    auto sJESnom = [](const VariantRatio& v) {
        return v.R_JES_nom;
    };
    auto sJESnomE = [](const VariantRatio& v) {
        return v.R_JES_nom_err;
    };
    auto drJES = [](const VariantRatio& v) {
        return v.jes_corr.value;
    };
    auto drJESE = [](const VariantRatio& v) {
        return v.jes_corr.err;
    };

    auto sJERvar = [](const VariantRatio& v) {
        return v.R_JER_var;
    };
    auto sJERvarE = [](const VariantRatio& v) {
        return v.R_JER_var_err;
    };
    auto sJERnom = [](const VariantRatio& v) {
        return v.R_JER_nom;
    };
    auto sJERnomE = [](const VariantRatio& v) {
        return v.R_JER_nom_err;
    };
    auto drJER = [](const VariantRatio& v) {
        return v.rjer_corr.value;
    };
    auto drJERE = [](const VariantRatio& v) {
        return v.rjer_corr.err;
    };

    // One plot per correction. The label goes into the file name.
    for (size_t k = 0; k < results.size(); ++k) {
        if (!results[k].valid) continue;
        std::vector<size_t> one = {k};
        const std::string& vl = results[k].label;
        DrawDualAxisInclusiveSummary(results, sJESvar, sJESvarE, sJESnom, sJESnomE, drJES, drJESE,
                                     "c_sum_jes_" + vl, "R_{JES} = <R>_{data} / <R>_{MC}",
                                     "DR (R_{JES}^{var}/R_{JES}^{nom})",
                                     summaryDir + "Summary_Scale_" + vl + ".pdf", wp_label, one);
        DrawDualAxisInclusiveSummary(
            results, sJERvar, sJERvarE, sJERnom, sJERnomE, drJER, drJERE, "c_sum_jer_" + vl,
            "R_{JER} = #sigma_{insitu}^{data} / #sigma_{insitu}^{MC}", "DR (R_{JER}^{var}/R_{JER}^{nom})",
            summaryDir + "Summary_Resolution_" + vl + ".pdf", wp_label, one);
    }

    // All corrections
    DrawDualAxisInclusiveSummary(results, sJESvar, sJESvarE, sJESnom, sJESnomE, drJES, drJESE,
                                 "c_sum_jes_all", "R_{JES} = <R>_{data} / <R>_{MC}",
                                 "DR (R_{JES}^{var}/R_{JES}^{nom})",
                                 summaryDir + "Summary_Scale_AllVariants.pdf", wp_label);
    DrawDualAxisInclusiveSummary(results, sJERvar, sJERvarE, sJERnom, sJERnomE, drJER, drJERE,
                                 "c_sum_jer_all", "R_{JER} = #sigma_{insitu}^{data} / #sigma_{insitu}^{MC}",
                                 "DR (R_{JER}^{var}/R_{JER}^{nom})",
                                 summaryDir + "Summary_Resolution_AllVariants.pdf", wp_label);

    // Subset from EASYBJETS_SUMMARY_VARIANTS: comma-separated legend labels
    // (VariantRatio::label, not the key), e.g.
    //   export EASYBJETS_SUMMARY_VARIANTS='bJR10v01,#mu-in-jet'
    {
        const char* envSub = std::getenv("EASYBJETS_SUMMARY_VARIANTS");
        if (envSub && *envSub) {
            std::vector<std::string> want;
            std::string s(envSub), tok;
            std::stringstream ss(s);
            while (std::getline(ss, tok, ',')) {
                size_t a = tok.find_first_not_of(" \t");
                size_t b = tok.find_last_not_of(" \t");
                if (a != std::string::npos) want.push_back(tok.substr(a, b - a + 1));
            }
            std::vector<size_t> sel;
            std::string tag;
            for (size_t k = 0; k < results.size(); ++k) {
                if (!results[k].valid) continue;
                for (const auto& w : want) {
                    if (results[k].label == w) {
                        sel.push_back(k);
                        tag += (tag.empty() ? "" : "_") + results[k].label;
                        break;
                    }
                }
            }
            if (!sel.empty()) {
                DrawDualAxisInclusiveSummary(
                    results, sJESvar, sJESvarE, sJESnom, sJESnomE, drJES, drJESE, "c_sum_jes_subset",
                    "R_{JES} = <R>_{data} / <R>_{MC}", "DR (R_{JES}^{var}/R_{JES}^{nom})",
                    summaryDir + "Summary_Scale_Subset_" + tag + ".pdf", wp_label, sel);
                DrawDualAxisInclusiveSummary(
                    results, sJERvar, sJERvarE, sJERnom, sJERnomE, drJER, drJERE, "c_sum_jer_subset",
                    "R_{JER} = #sigma_{insitu}^{data} / #sigma_{insitu}^{MC}",
                    "DR (R_{JER}^{var}/R_{JER}^{nom})",
                    summaryDir + "Summary_Resolution_Subset_" + tag + ".pdf", wp_label, sel);
                Log::info("Ratio_Response_Inclusive", "Custom-subset summary plots written for: " + tag);
            } else {
                Log::error("Ratio_Response_Inclusive",
                           std::string("EASYBJETS_SUMMARY_VARIANTS matched no variants: ") + envSub);
            }
        }
    }
}

// Two-axis summary, x = correction.
//   left axis : data/MC ratio (R_JES or R_JER), nominal in black at x = 0,
//               each correction filled at x - 0.10
//   right axis: DR correction / nominal (corr-aware), same colour with a
//               shifted marker style at x + 0.10
inline void DrawDualAxisInclusiveSummary(const std::vector<VariantRatio>& results,
                                         std::function<double(const VariantRatio&)> getSingleVar,
                                         std::function<double(const VariantRatio&)> getSingleVarErr,
                                         std::function<double(const VariantRatio&)> getSingleNom,
                                         std::function<double(const VariantRatio&)> getSingleNomErr,
                                         std::function<double(const VariantRatio&)> getDR,
                                         std::function<double(const VariantRatio&)> getDRErr,
                                         const std::string& canvasName, const std::string& leftTitle,
                                         const std::string& rightTitle, const std::string& outputPath,
                                         const std::string& wpLabel, const std::vector<size_t>& selected) {
    if (results.empty()) return;

    std::vector<size_t> idxs;
    if (selected.empty()) {
        for (size_t i = 0; i < results.size(); ++i)
            if (results[i].valid) idxs.push_back(i);
    } else {
        for (size_t k : selected)
            if (k < results.size() && results[k].valid) idxs.push_back(k);
    }
    if (idxs.empty()) return;

    const int nVar = (int)idxs.size();

    // Colour and marker by position in the plot
    static const int kPalette[] = {kRed + 1,  kBlue + 1,   kGreen + 2,  kMagenta + 1, kOrange + 7,
                                   kCyan + 2, kViolet - 1, kSpring + 4, kAzure + 1,   kPink + 7};
    const int nPal = (int)(sizeof(kPalette) / sizeof(kPalette[0]));
    static const int kMarker[] = {20, 21, 22, 23, 29, 33, 34, 47, 43, 45};
    const int nMk = (int)(sizeof(kMarker) / sizeof(kMarker[0]));

    TCanvas* c = new TCanvas(canvasName.c_str(), canvasName.c_str(), 900, 650);
    c->SetLeftMargin(0.13);
    c->SetBottomMargin(0.18);
    c->SetTopMargin(0.08);
    c->SetRightMargin(0.14);

    double nomS = getSingleNom(results[idxs[0]]);
    double nomSErr = std::abs(getSingleNomErr(results[idxs[0]]));

    // Left axis range, with extra room at the top for the legend
    double sMin = nomS - nomSErr, sMax = nomS + nomSErr;
    for (size_t k : idxs) {
        double v = getSingleVar(results[k]), e = std::abs(getSingleVarErr(results[k]));
        if (std::isfinite(v) && std::isfinite(e)) {
            sMin = std::min(sMin, v - e);
            sMax = std::max(sMax, v + e);
        }
    }
    if (sMin > sMax) {
        sMin = 0.9;
        sMax = 1.1;
    }
    double sRange = sMax - sMin;
    if (sRange < 1e-6) sRange = 0.1;
    double sLo = sMin - 0.25 * sRange;
    double sHi = sMax + 0.45 * sRange;

    // Right axis range
    double dMin = 1e9, dMax = -1e9;
    for (size_t k : idxs) {
        double v = getDR(results[k]), e = std::abs(getDRErr(results[k]));
        if (std::isfinite(v) && std::isfinite(e)) {
            dMin = std::min(dMin, v - e);
            dMax = std::max(dMax, v + e);
        }
    }
    if (dMin > dMax) {
        dMin = 0.95;
        dMax = 1.05;
    }
    double dRange = dMax - dMin;
    if (dRange < 1e-6) dRange = 0.1;
    double dLo = dMin - 0.25 * dRange;
    double dHi = dMax + 0.25 * dRange;

    auto mapDR = [&](double dr) {
        return sLo + (dr - dLo) / (dHi - dLo) * (sHi - sLo);
    };

    TH1F* frame = new TH1F((canvasName + "_frame").c_str(), "", nVar + 1, -0.5, nVar + 0.5);
    frame->SetMinimum(sLo);
    frame->SetMaximum(sHi);
    frame->GetYaxis()->SetTitle(leftTitle.c_str());
    frame->GetYaxis()->SetTitleSize(0.038);
    frame->GetYaxis()->SetTitleOffset(1.45);
    frame->GetYaxis()->SetLabelSize(0.035);
    frame->GetXaxis()->SetLabelSize(0.0); // x labels drawn below with TLatex, rotated
    frame->SetStats(0);
    frame->Draw("AXIS");

    // Dotted line at data/MC = 1
    if (sLo < 1.0 && sHi > 1.0) {
        TLine* ref = new TLine(-0.5, 1.0, nVar + 0.5, 1.0);
        ref->SetLineStyle(3);
        ref->SetLineColor(kGray + 1);
        ref->Draw("same");
    }

    TLegend* leg = new TLegend(0.50, 0.68, 0.86, 0.90);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.024);
    leg->SetNColumns(2);
    leg->SetColumnSeparation(0.05);

    // Nominal data/MC
    {
        double xs[1] = {0.0}, ys[1] = {nomS}, exs[1] = {0.0}, eys[1] = {nomSErr};
        TGraphErrors* g = new TGraphErrors(1, xs, ys, exs, eys);
        g->SetMarkerStyle(21);
        g->SetMarkerSize(1.6);
        g->SetMarkerColor(kBlack);
        g->SetLineColor(kBlack);
        g->SetLineWidth(2);
        g->Draw("P same");
        leg->AddEntry(g, "Nominal data/MC", "lp");
        leg->AddEntry((TObject*)nullptr, "", ""); // empty cell so each correction's pair stays on one row
    }

    TLatex xlab;
    xlab.SetTextSize(0.030);
    xlab.SetTextAlign(22);
    xlab.SetTextAngle(35);
    double xLabelY = sLo - 0.06 * (sHi - sLo);

    xlab.DrawLatex(0.0, xLabelY, "Nominal");

    for (int slot = 0; slot < nVar; ++slot) {
        size_t k = idxs[slot];
        const VariantRatio& r = results[k];
        int col = kPalette[slot % nPal];
        int mk = kMarker[slot % nMk];
        double xc = slot + 1;

        // data/MC, left axis
        {
            double xs[1] = {xc - 0.10}, ys[1] = {getSingleVar(r)}, exs[1] = {0.0},
                   eys[1] = {std::abs(getSingleVarErr(r))};
            TGraphErrors* gS = new TGraphErrors(1, xs, ys, exs, eys);
            gS->SetMarkerStyle(mk);
            gS->SetMarkerSize(1.5);
            gS->SetMarkerColor(col);
            gS->SetLineColor(col);
            gS->SetLineWidth(2);
            gS->Draw("P same");
            leg->AddEntry(gS, (r.label + " data/MC").c_str(), "lp");
        }
        // DR, mapped onto the left axis; the error is scaled the same way
        {
            double xs[1] = {xc + 0.10}, ys[1] = {mapDR(getDR(r))}, exs[1] = {0.0},
                   eys[1] = {std::abs(getDRErr(r)) / (dHi - dLo) * (sHi - sLo)};
            TGraphErrors* gD = new TGraphErrors(1, xs, ys, exs, eys);
            gD->SetMarkerStyle(mk + 4);
            gD->SetMarkerSize(1.5);
            gD->SetMarkerColor(col);
            gD->SetLineColor(col);
            gD->SetLineWidth(1);
            gD->SetLineStyle(2);
            gD->Draw("P same");
            leg->AddEntry(gD, (r.label + " DR").c_str(), "lp");
        }
        xlab.DrawLatex(xc, xLabelY, r.label.c_str());
    }

    TGaxis* axR = new TGaxis(nVar + 0.5, sLo, nVar + 0.5, sHi, dLo, dHi, 510, "+L");
    axR->SetTitle(rightTitle.c_str());
    axR->SetLabelSize(0.035);
    axR->SetTitleSize(0.038);
    axR->SetTitleOffset(1.4);
    axR->Draw();

    leg->Draw();
    AtlasUtils::DrawAtlasLabel(0.16, 0.88);
    TLatex wpTx;
    wpTx.SetNDC();
    wpTx.SetTextFont(42);
    wpTx.SetTextSize(0.032);
    wpTx.DrawLatex(0.16, 0.81, Form("WP: %s", wpLabel.c_str()));

    c->SaveAs(outputPath.c_str());
    delete frame;
    delete c;
    delete leg;
    delete axR;
}

// Entry point called by the channel, variant = "ALL" (or empty) runs
// every correction of the channel; otherwise give one key such as
// "DB_bJER10v00". Loops over AnalysisWPConfig::StandardWPs().
inline void Ratio_Response_Inclusive_Run(const std::string& variant, const std::string& channel) {
    gROOT->SetBatch(kTRUE);
    TH1::AddDirectory(kFALSE);
    AtlasUtils::LoadAtlasStyle();

    std::vector<std::string> variantKeys;
    if (variant == "ALL" || variant.empty()) {
        for (const auto& k : GetPrimaryVariantKeys(channel))
            if (k != channel) variantKeys.push_back(k);
    } else {
        variantKeys.push_back(variant);
    }

    if (variantKeys.empty()) {
        Log::error("Ratio_Response_Inclusive", "No variants to process.");
        return;
    }

    Log::info("Ratio_Response_Inclusive", "Processing variants:");
    for (const auto& k : variantKeys) Log::info("  ", k);

    for (const auto& wp : AnalysisWPConfig::StandardWPs()) {
        ProcessWP(variantKeys, wp, channel);
    }

    Log::info("Ratio_Response_Inclusive", "All done.");
}

#endif // RATIO_RESPONSE_INCLUSIVE_H