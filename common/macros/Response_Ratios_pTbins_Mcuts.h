/**
 * ============================================================================
 *  Response_Ratios_pTbins_Mcuts.h - b-JES and b-JER double ratios per pT and
 *                                   mass bin (correction/nominal)
 * ============================================================================
 * ===========================================
 * @author:    Thiziri AMEZZA
 * @email:     thiziri.amezza@cern.ch
 * @date:      10 June 2026 
 * ===========================================
 *
 * Reads the response fit results of one correction (bJER10v00, bJER10v01, MuonCorr) 
 * and of the nominal calibration, and compares them cell by cell in
 * (mass bin, pT bin). Shared by the DB (Z(bb)+gamma, ZbbGamma) and MJB
 * (Z(bb)+jets, ZbbJets) channels.
 *
 * Two mass modes, picked by the third argument of Ratio_Response_Run:
 * 
 *   "cuts"       (default) three mass bins 50-80, 80-110, >=110 GeV, read from
 *                the Response_Fit_pTbins_Mcuts output
 * 
 *   "inclusive"  one "MassInclusive" bin, read from the Response_Fit_pTbins
 *                output; results go to a separate folder (R2pTbin_MassInclusive)
 * 
 * pT bins: DB 200-300, 300-450 GeV; MJB 450-600, 600-900, 900-1200 GeV 
 *          (the last bin has no upper limit in the selection; 450 / 1200 is for display).
 *
 * Definitions, with R = pT(jet)/pT(reference) and mu, sigma from the fits:
 * 
 *   R_JES            = mu_data / mu_combMC
 *   sigma_insitu     = sqrt(sigma_data^2   - sigma_truth^2)
 *   sigma^MC_insitu  = sqrt(sigma_combMC^2 - sigma_truth^2)
 *   R_JER            = sigma_insitu / sigma^MC_insitu
 *   DR(X)            = X_correction / X_nominal   (the "double ratio")
 * 
 * sigma_truth comes from the pT(truth)/pT(reference) fit of the same variant.
 * The correction/nominal ratios of mu and sigma of each sample (data, signal, 
 * background, combined MC) are also drawn.
 *
 * Errors of the double ratios:``
 *
 *   uncorrelated   All fit errors (data and MC, correction and nominal) are
 *                  treated as independent and added in quadrature. This ignores
 *                  that correction and nominal use the same events, so it is the
 *                  largest error: a conservative upper bound.
 * 
 *  corr-aware     Correction and nominal fluctuate together (same events), so
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
 * How to run (from ZbbGamma/macros or ZbbJets/macros, after the matching fit):
 * 
 *   export EASYBJETS_RUN=run3          same as for the fit
 *   root 'Response_Ratios_pTbins_Mcuts.cxx()'                 3 mass bins, all corrections
 *   root 'Response_Ratios_pTbins_Mcuts.cxx("DB_bJER10v00")'   one correction (MJB_... in ZbbJets)
 *   root 'Response_Ratios_pTbins.cxx()'                       inclusive mass, all corrections
 * 
 * Then run Response_Ratios_Combined.cxx and Response_Summary_Combined.cxx (after the 3-mass-bin mode) 
 * or Response_Summary_Combined_pTbins.cxx (after the inclusive-mass mode), in both channels, with the same options.
 *
 * Input:  <PlotDir>/<variant>/Response_Fit_pTbins_Mcuts_<DB|MJB>/<wp>/<mass dir>/
 *           {data,mc,bkg,combined_MC,truth_ref}_fit_parameters_ptbin<i>.txt
 *         for the nominal and the correction; <mass dir> is 50-80GeV,
 *         80-110GeV, MassGE110GeV or MassInclusive. <PlotDir> is ../Run3Plots
 *         or ../Run2Plots, plus the fit-model suffix if any.
 * Output: <PlotDir>/Response/Ratios/<DB|MJB>_Ratio_<variant key>/<wp>/R2pTbin/
 *         (R2pTbin_MassInclusive/ in inclusive mode):
 *           heat maps of the per-sample mu and sigma ratios and of DR(R_JES),
 *           DR(sigma_insitu), DR(R_JER) (uncorrelated and _CorrAware),
 *           Summary_*_vs_pT_<wp>.pdf, *_ErrorComparison_vs_pT_<wp>.pdf,
 *           Summary_Ratios_<wp>.txt   readable table
 *           Summary_Ratios_<wp>.csv   read by Response_Ratios_Combined.h (for
 *                                     Response_Ratios_Combined.cxx and
 *                                     Response_Summary_Combined.cxx) and by
 *                                     Response_Summary_Combined_Inclusive.h (for
 *                                     Response_Summary_Combined_pTbins.cxx)
 *
 * Options (environment variables):
 *   EASYBJETS_RUN                 run2 | run3 (default run3): ../Run2Plots or ../Run3Plots
 *   EASYBJETS_FIT_MODEL           unset = Gaussian; bukin, cb, dscb, dgauss, all add a suffix
 *                                 to the plot folder, so it must match the fit step
 *   EASYBJETS_CORR_METHOD         rho1 (default) | file | bootstrap: source of the "full" rho
 *   EASYBJETS_CORRELATIONS        correlation file (default ./correlations.txt)
 *   EASYBJETS_HEATMAP_JES_RANGE   "lo,hi" colour range of the DR(R_JES) heat maps (default 0.90,1.10)
 *   EASYBJETS_HEATMAP_JER_RANGE   "lo,hi" colour range of the DR(sigma_insitu) and DR(R_JER)
 *                                 heat maps (default 0.60,1.40)
 * EASYBJETS_BINNING does not change any folder here; it only matters for the fit.
 *
 * ==========================================================================
 * @version:  1.0
 * ==========================================================================
 */

#ifndef RATIO_RESPONSE_H
#define RATIO_RESPONSE_H

#include "VariantConfig.h"
#include "CommonUtils.h" // also provides FitUtils::ComputeJERRatio
#include "CorrelationConfig.h"

#include <fstream>

#include "TPaletteAxis.h"


// Binning. The pT bins depend on the channel (DB 2 bins, MJB 3 bins), so they
// are reached through the CurrentNPtBins()/CurrentPtBins() helpers. Arrays
// are sized with MAX_PT_BINS and loops only run over CurrentNPtBins().
// The channel is set once by Ratio_Response_Run and kept in CurrentChannel();
// each macro runs in its own session, so this global state is safe.
constexpr int MAX_PT_BINS = 3;

inline std::string& CurrentChannel() {
    static std::string ch = "DB";
    return ch;
}

inline int CurrentNPtBins() {
    return (CurrentChannel() == "MJB") ? BinConfig::nPtBins_3 : BinConfig::nPtBins_DB;
}

inline const double* CurrentPtBins() {
    return (CurrentChannel() == "MJB") ? BinConfig::ptBins_3 : BinConfig::ptBins_DB;
}

// Same edges as CurrentPtBins(), but the open-ended last edge (1e9) is
// replaced by 450 (DB) or 1200 (MJB). This macro only reads fit results and
// never selects events, so every pT edge here is for display and should come
// from this helper.
inline const double* CurrentPtBinsDisplay() {
    return (CurrentChannel() == "MJB") ? BinConfig::ptBins_3_display : BinConfig::ptBins_DB_display;
}

inline const char* CurrentPtBinLabel(int i) {
    return (CurrentChannel() == "MJB") ? BinConfig::PtBinLabel_3(i) : BinConfig::PtBinLabel_DB(i);
}
// Y-axis padding of the ratio plots: the MJB ratios spread more than the DB
// ones, so they get more room.
inline double RatioPadFactor() { return (CurrentChannel() == "MJB") ? 0.55 : 0.15; }
inline double RatioYLoClamp() { return (CurrentChannel() == "MJB") ? 0.98 : 0.95; }
inline double RatioYHiClamp() { return (CurrentChannel() == "MJB") ? 1.02 : 1.05; }
// Z range used when a heatmap has no filled cell.
inline double EmptyYMin() { return (CurrentChannel() == "MJB") ? 0.85 : 0.9; }
inline double EmptyYMax() { return (CurrentChannel() == "MJB") ? 1.3 : 1.1; }


// Mass bins, chosen at run time like in the fit macro: "cuts" gives 3 bins,
// "inclusive" one bin read from the "MassInclusive" fit folder. nMassCuts is
// a run-time count and the arrays are sized with MAX_MASS_CUTS.
// ConfigureMassBins(), called from Ratio_Response_Run, fills them.
constexpr int MAX_MASS_CUTS = 6;
int nMassCuts = 3;
const char* massCutLabels[MAX_MASS_CUTS] = {"50-80 GeV", "80-110 GeV", "#geq110 GeV"};
// Folder names. They must match fileSystemMassCutLabels in
// Response_Fit_pTbins_Mcuts.h (input) and MassTag() in Response_Summary_Combined.h
// (reads the CSV written here).
const char* massCutDirNames[MAX_MASS_CUTS] = {"50-80GeV", "80-110GeV", "MassGE110GeV"};
double massCutLow[MAX_MASS_CUTS] = {50, 80, 110};
double massCutHigh[MAX_MASS_CUTS] = {80, 110, -1};

// false = 3 mass cuts, true = one inclusive-mass bin.
inline bool& MassInclusiveMode() {
    static bool m = false;
    return m;
}

inline void ConfigureMassBins() {
    if (MassInclusiveMode()) {
        nMassCuts = 1;
        massCutLow[0] = -1;
        massCutHigh[0] = -1;
        massCutLabels[0] = "Inclusive mass";
        massCutDirNames[0] = "MassInclusive";
    } else {
        nMassCuts = 3;
        massCutLow[0] = 50;
        massCutLow[1] = 80;
        massCutLow[2] = 110;
        massCutHigh[0] = 80;
        massCutHigh[1] = 110;
        massCutHigh[2] = -1;
        massCutLabels[0] = "50-80 GeV";
        massCutLabels[1] = "80-110 GeV";
        massCutLabels[2] = "#geq110 GeV";
        massCutDirNames[0] = "50-80GeV";
        massCutDirNames[1] = "80-110GeV";
        massCutDirNames[2] = "MassGE110GeV";
    }
}


// Physics formulas

// sigma_insitu = sqrt(sigma_data^2 - sigma_truth^2)  (data side).
// The second argument is sigma_truth from the pT(truth)/pT(ref) fit, not the
// sigma of the reco MC response. Returns 0 when the difference is negative.
double CalculateResolution(double data_sigma, double truth_sigma) {
    double sigma_sq_diff = data_sigma * data_sigma - truth_sigma * truth_sigma;
    if (sigma_sq_diff <= 0) return 0.0;
    return std::sqrt(sigma_sq_diff);
}

// Error on sigma_insitu, with sigma_data and sigma_truth taken as independent:
//   d(sigma_insitu)^2 = (sigma_data/sigma_insitu * d sigma_data)^2
//                     + (sigma_truth/sigma_insitu * d sigma_truth)^2
double CalculateResolutionError(double data_sigma, double data_sigma_err, double truth_sigma,
                                double truth_sigma_err) {
    double A = data_sigma * data_sigma - truth_sigma * truth_sigma;
    if (A <= 0) return 0.0;

    double denom = std::sqrt(A);
    double dfd_data = data_sigma / denom;
    double dfd_truth = -truth_sigma / denom;

    double error_sq = dfd_data * dfd_data * data_sigma_err * data_sigma_err
                      + dfd_truth * dfd_truth * truth_sigma_err * truth_sigma_err;
    return std::sqrt(error_sq);
}

// A/B, or 0 when B is 0 (an empty cell then shows up as 0, not inf).
inline double SafeRatio(double A, double B) { return (B != 0) ? A / B : 0.0; }

// Error on A/B with A and B uncorrelated:
//   d(A/B) = |A/B| * sqrt((dA/A)^2 + (dB/B)^2)
// Variant and nominal come from the same events, so for the double ratios
// this overestimates the error.
inline double SafeRatioErr(double A, double A_err, double B, double B_err) {
    if (B == 0 || A == 0) return 0.0;
    double ratio = A / B;
    return std::abs(ratio) * std::sqrt(std::pow(A_err / A, 2) + std::pow(B_err / B, 2));
}


// Correlation-aware error on DR(sigma_insitu) = sigma_insitu_var / sigma_insitu_nom.
// sigma_data comes from the same data events for variant and nominal, so its
// fit error is assumed to cancel in the ratio and is dropped. The sigma_truth
// fit errors are treated as uncorrelated between variant and nominal (each
// fit has its own window and binning) and added in quadrature:
//   d sigma_insitu / d sigma_truth = -sigma_truth / sigma_insitu
//   (dDR/DR)^2 = (sigma_truth_var * d sigma_truth_var / sigma_insitu_var^2)^2
//              + (sigma_truth_nom * d sigma_truth_nom / sigma_insitu_nom^2)^2
inline double CorrelationAwareJERRatioErr(double jer_var, double sigma_data_var, double sigma_truth_var,
                                          double sigma_truth_var_err, double jer_nom, double sigma_data_nom,
                                          double sigma_truth_nom, double sigma_truth_nom_err) {
    if (jer_var == 0 || jer_nom == 0) return 0.0;
    double ratio = jer_var / jer_nom;

    double diff_var = sigma_data_var * sigma_data_var - sigma_truth_var * sigma_truth_var;
    double diff_nom = sigma_data_nom * sigma_data_nom - sigma_truth_nom * sigma_truth_nom;
    if (diff_var <= 0 || diff_nom <= 0) return 0.0;

    double denom_var = std::sqrt(diff_var);
    double denom_nom = std::sqrt(diff_nom);

    double dJER_dtruth_var = sigma_truth_var / denom_var;
    double dJER_dtruth_nom = sigma_truth_nom / denom_nom;

    double rel_var = dJER_dtruth_var * sigma_truth_var_err / std::abs(jer_var);
    double rel_nom = dJER_dtruth_nom * sigma_truth_nom_err / std::abs(jer_nom);

    return std::abs(ratio) * std::sqrt(rel_var * rel_var + rel_nom * rel_nom);
}


// Correlation-aware error on DR(R_JER) = R_JER_var / R_JER_nom, with
// R_JER = sqrt((s_d^2 - s_t^2) / (s_m^2 - s_t^2)) (s_d data, s_m combined MC,
// s_t truth sigma). Same idea as above: the s_d and s_m fit errors are taken
// to cancel between variant and nominal, only the s_t errors are kept.
//   dR/ds_t = s_t (s_d^2 - s_m^2) / (sigma_insitu * sigma^MC_insitu^3)
//   dR/R    = |s_t (s_d^2 - s_m^2)| * ds_t / ((s_d^2 - s_t^2)(s_m^2 - s_t^2))
//   (dDR/DR)^2 = (dR_var/R_var)^2 + (dR_nom/R_nom)^2
inline double RJERRelErrFromSigmaTruth(double sigma_d, double sigma_m, double sigma_t, double sigma_t_err) {
    double u = sigma_d * sigma_d - sigma_t * sigma_t; // sigma_insitu^2
    double v = sigma_m * sigma_m - sigma_t * sigma_t; // (sigma^MC_insitu)^2
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


// Correlation key of one (mass bin, pT bin) cell, as written by
// Bootstrap_Correlations: the BinConfig pT label, plus the mass bin unless the
// run has no mass split (CorrelationConfig::RespBinLabel).
inline std::string RespCorrLabel(int massIdx, int ptBin) {
    const std::string ptLab = (CurrentChannel() == "MJB") ? BinConfig::PtBinFileLabel_3(ptBin)
                                                          : BinConfig::PtBinFileLabel_DB(ptBin);
    return CorrelationConfig::RespBinLabel(ptLab, massCutDirNames[massIdx]);
}

// Reader for the *_fit_parameters_ptbin<i>.txt files written by
// FitOutputUtils::SaveFitParameters in Response_Fit_pTbins_Mcuts.h.
struct FitParams {
    double mean = 0;
    double mean_err = 0;
    double sigma = 0;
    double sigma_err = 0;
};

// Reads the "Mean:" and "Sigma:" lines ("value +/- error" or "value +/- error").
bool ReadFitParams(const std::string& filepath, FitParams& params) {
    std::ifstream infile(filepath);
    if (!infile.is_open()) {
        std::cerr << "Warning: Could not open fit file: " << filepath << std::endl;
        return false;
    }

    std::string line;
    while (std::getline(infile, line)) {
        if (line.substr(0, 5) == "Mean:") {
            size_t colon = line.find(':');
            if (colon == std::string::npos) continue;
            std::string rest = line.substr(colon + 1);

            size_t sep = rest.find("±");
            if (sep == std::string::npos) sep = rest.find("+/-");
            if (sep == std::string::npos) {
                sscanf(rest.c_str(), " %lf", &params.mean);
                continue;
            }

            std::string val_str = rest.substr(0, sep);
            size_t after_sep = sep;
            if (rest[sep] == '+')
                after_sep += 3;
            else
                after_sep += 2;
            while (after_sep < rest.size() && !std::isdigit(rest[after_sep]) && rest[after_sep] != '-')
                after_sep++;
            std::string err_str = rest.substr(after_sep);

            params.mean = std::atof(val_str.c_str());
            params.mean_err = std::atof(err_str.c_str());
        } else if (line.substr(0, 6) == "Sigma:") {
            size_t colon = line.find(':');
            if (colon == std::string::npos) continue;
            std::string rest = line.substr(colon + 1);

            size_t sep = rest.find("±");
            if (sep == std::string::npos) sep = rest.find("+/-");
            if (sep == std::string::npos) {
                sscanf(rest.c_str(), " %lf", &params.sigma);
                continue;
            }

            std::string val_str = rest.substr(0, sep);
            size_t after_sep = sep;
            if (rest[sep] == '+')
                after_sep += 3;
            else
                after_sep += 2;
            while (after_sep < rest.size() && !std::isdigit(rest[after_sep]) && rest[after_sep] != '-')
                after_sep++;
            std::string err_str = rest.substr(after_sep);

            params.sigma = std::atof(val_str.c_str());
            params.sigma_err = std::atof(err_str.c_str());
        }
    }
    infile.close();
    return true;
}


// Drawing helpers

// Symmetric z range around 1 that covers every non-empty cell.
void SetDynamicRatioPaletteRange(TH2F* hist) {
    double minVal = 1e18, maxVal = -1e18;
    bool hasContent = false;
    for (int i = 1; i <= hist->GetNbinsX(); ++i) {
        for (int j = 1; j <= hist->GetNbinsY(); ++j) {
            double content = hist->GetBinContent(i, j);
            if (content != 0) {
                minVal = std::min(minVal, content);
                maxVal = std::max(maxVal, content);
                hasContent = true;
            }
        }
    }
    if (!hasContent) {
        hist->SetMinimum(EmptyYMin());
        hist->SetMaximum(EmptyYMax());
        return;
    }

    double center = 1.0;
    double maxAbsDist = std::max(std::abs(center - minVal), std::abs(maxVal - center));
    if (maxAbsDist < 0.01) maxAbsDist = 0.01;
    double paddedExtent = maxAbsDist * 1.02;
    hist->SetMinimum(center - paddedExtent);
    hist->SetMaximum(center + paddedExtent);
}

// Fixed colour range "lo,hi" for a double-ratio heatmap, read from envName.
// Corrections and Run periods are separate jobs, so the range has to be fixed
// for their heatmaps to be comparable. Cells outside the range take the end
// colour. Used with EASYBJETS_HEATMAP_JES_RANGE (default 0.90,1.10) and
// EASYBJETS_HEATMAP_JER_RANGE (default 0.60,1.40).
inline void GetHeatmapRange(const char* envName, double defLo, double defHi, double& lo, double& hi) {
    lo = defLo;
    hi = defHi;
    if (const char* e = std::getenv(envName)) {
        double a = 0, b = 0;
        if (std::sscanf(e, "%lf,%lf", &a, &b) == 2 && b > a) {
            lo = a;
            hi = b;
        } else
            std::cerr << "[Heatmap] WARN: bad " << envName << "='" << e
                      << "' (want \"lo,hi\"); using default\n";
    }
}

// Heatmap of one quantity, mass bin on x and pT bin on y. zmax > zmin gives a
// fixed colour range, otherwise SetDynamicRatioPaletteRange is used.
void Draw2DHeatmap(const std::string& canvasName, const std::string& title, const std::string& paletteLabel,
                   const std::string& outputPath, double values[MAX_MASS_CUTS][MAX_PT_BINS],
                   double errors[MAX_MASS_CUTS][MAX_PT_BINS], double zmin = 0.0, double zmax = 0.0) {

    TCanvas* canvas = new TCanvas(canvasName.c_str(), title.c_str(), 800, 600);
    canvas->SetTopMargin(0.12);
    canvas->SetRightMargin(0.15);
    canvas->SetLeftMargin(0.12);
    canvas->SetBottomMargin(0.12);
    canvas->cd();

    TString mapName = TString(canvasName.c_str()) + "_map";
    TH2F* map = new TH2F(mapName, ";Mass Bin;p_{T} Bin", nMassCuts, -0.5, nMassCuts - 0.5, CurrentNPtBins(),
                         -0.5, CurrentNPtBins() - 0.5);

    for (int m = 0; m < nMassCuts; m++) {
        for (int bin = 0; bin < CurrentNPtBins(); bin++) {
            map->SetBinContent(m + 1, bin + 1, values[m][bin]);
            map->SetBinError(m + 1, bin + 1, errors[m][bin]);
        }
    }

    map->SetStats(0);
    map->SetMarkerSize(1.2);

    for (int m = 0; m < nMassCuts; m++) map->GetXaxis()->SetBinLabel(m + 1, massCutLabels[m]);
    for (int bin = 0; bin < CurrentNPtBins(); bin++) {
        TString label =
            Form("%d-%d GeV", (int)CurrentPtBinsDisplay()[bin], (int)CurrentPtBinsDisplay()[bin + 1]);
        map->GetYaxis()->SetBinLabel(bin + 1, label.Data());
    }

    map->GetXaxis()->SetLabelSize(0.025);
    map->GetYaxis()->SetLabelSize(0.025);
    map->GetXaxis()->SetTitleSize(0.035);
    map->GetYaxis()->SetTitleSize(0.035);
    map->GetZaxis()->SetTitleSize(0.03);
    map->SetTitle("");

    gStyle->SetNumberContours(100);
    if (zmax > zmin) {
        map->SetMinimum(zmin);
        map->SetMaximum(zmax);
    } else
        SetDynamicRatioPaletteRange(map);
    gStyle->SetPaintTextFormat("4.3f");

    map->Draw("COLZ TEXT");
    gPad->Update();

    TPaletteAxis* palette = (TPaletteAxis*)map->GetListOfFunctions()->FindObject("palette");
    if (palette) {
        palette->SetX1NDC(0.86);
        palette->SetX2NDC(0.89);
        palette->SetY1NDC(0.12);
        palette->SetY2NDC(0.88);
        palette->SetLabelSize(0.03);
        palette->SetTitleOffset(1.3);

        TLatex ptext;
        ptext.SetNDC();
        ptext.SetTextSize(0.035);
        ptext.SetTextAlign(22);
        ptext.SetTextAngle(90);
        ptext.DrawLatex(0.97, 0.5, paletteLabel.c_str());
    }

    AtlasUtils::DrawAtlasLabel(0.18, 0.85);
    canvas->SaveAs(outputPath.c_str());

    delete map;
    delete canvas;
}


// Ratio vs pT, one marker series per mass bin (shifted by -20/0/+20 GeV so they do not overlap). 
// y_lo_force/y_hi_force fix the y range when given.
void MakeRatioSummaryPlot_vs_pT(const std::string& outputDir, const std::string& canvasName,
                                const std::string& yTitle, const std::string& saveName,
                                double data[MAX_MASS_CUTS][MAX_PT_BINS],
                                double errs[MAX_MASS_CUTS][MAX_PT_BINS],
                                double y_lo_force = std::numeric_limits<double>::quiet_NaN(),
                                double y_hi_force = std::numeric_limits<double>::quiet_NaN()) {
    Int_t colors[MAX_MASS_CUTS] = {kRed, kBlue, kGreen + 2};
    Int_t markers[MAX_MASS_CUTS] = {20, 21, 22};

    double ptCenters[MAX_PT_BINS];
    for (int i = 0; i < CurrentNPtBins(); i++)
        ptCenters[i] = (CurrentPtBinsDisplay()[i] + CurrentPtBinsDisplay()[i + 1]) / 2.0;
    double offsets[MAX_MASS_CUTS] = {-20.0, 0.0, 20.0};

    // Y range from the points. An end that falls inside [0.95, 1.05] is moved
    // to RatioYLoClamp()/RatioYHiClamp() so flat ratios do not fill the frame.
    double ymin_val = 1e6, ymax_val = -1e6;
    for (int m = 0; m < nMassCuts; m++) {
        for (int bin = 0; bin < CurrentNPtBins(); bin++) {
            double lo = data[m][bin] - errs[m][bin];
            double hi = data[m][bin] + errs[m][bin];
            if (lo < ymin_val) ymin_val = lo;
            if (hi > ymax_val) ymax_val = hi;
        }
    }
    double range = ymax_val - ymin_val;
    double padding = range * RatioPadFactor();

    double y_lo = ymin_val - padding;
    double y_hi = ymax_val + padding;
    if (y_lo > 0.95) y_lo = RatioYLoClamp();
    if (y_hi < 1.05) y_hi = RatioYHiClamp();

    // The R_JER plots pass a fixed [0.45, 2.05].
    if (!std::isnan(y_lo_force)) y_lo = y_lo_force;
    if (!std::isnan(y_hi_force)) y_hi = y_hi_force;


    TCanvas* canvas = new TCanvas(canvasName.c_str(), canvasName.c_str(), 800, 600);
    canvas->SetLeftMargin(0.15);
    canvas->SetBottomMargin(0.15);
    canvas->SetTopMargin(0.12);
    canvas->SetRightMargin(0.1);

    TH1F* frame = new TH1F((canvasName + "_frame").c_str(), "", CurrentNPtBins(), CurrentPtBinsDisplay());
    frame->SetMinimum(y_lo);
    frame->SetMaximum(y_hi);
    frame->GetXaxis()->SetTitle("p_{T} [GeV]");
    frame->GetYaxis()->SetTitle(yTitle.c_str());
    frame->GetXaxis()->SetTitleSize(0.04);
    frame->GetYaxis()->SetTitleSize(0.04);
    frame->GetXaxis()->SetLabelSize(0.04);
    frame->GetYaxis()->SetLabelSize(0.04);
    frame->GetYaxis()->SetTitleOffset(1.55);
    frame->SetStats(0);
    frame->GetXaxis()->SetRangeUser(CurrentPtBinsDisplay()[0], CurrentPtBinsDisplay()[CurrentNPtBins()]);
    frame->GetXaxis()->SetNdivisions(0, 0, 0, kFALSE);
    frame->GetXaxis()->SetLabelSize(0);
    frame->Draw("AXIS");

    double yrange = y_hi - y_lo;
    TLatex* tick_latex = new TLatex();
    tick_latex->SetTextSize(0.04);
    tick_latex->SetTextAlign(23);
    tick_latex->SetTextColor(kBlack);
    for (int i = 0; i <= CurrentNPtBins(); i++) {
        TLine* tick =
            new TLine(CurrentPtBinsDisplay()[i], y_lo, CurrentPtBinsDisplay()[i], y_lo + 0.01 * yrange);
        tick->SetLineColor(kBlack);
        tick->Draw("same");
        tick_latex->DrawLatex(CurrentPtBinsDisplay()[i], y_lo - 0.05 * yrange,
                              Form("%.0f", CurrentPtBinsDisplay()[i]));
    }

    TLatex* title_latex = new TLatex();
    title_latex->SetTextSize(0.045);
    title_latex->SetTextAlign(22);
    title_latex->SetNDC();
    title_latex->DrawLatex(0.85, 0.04, "p_{T} [GeV]");

    TLine* refline = new TLine(CurrentPtBinsDisplay()[0], 1.0, CurrentPtBinsDisplay()[CurrentNPtBins()], 1.0);
    refline->SetLineStyle(2);
    refline->SetLineColor(kBlack);
    refline->SetLineWidth(2);
    refline->Draw("same");

    gStyle->SetLegendBorderSize(1);
    TLegend* legend = new TLegend(0.18, 0.65, 0.5, 0.78);
    legend->SetBorderSize(1);
    legend->SetLineColor(kBlack);
    legend->SetLineWidth(1);
    legend->SetFillStyle(0);
    legend->SetTextSize(0.028);

    TLegend* procLeg = new TLegend(0.62, 0.75, 0.89, 0.89);
    procLeg->SetBorderSize(0);
    procLeg->SetFillStyle(0);
    procLeg->SetTextSize(0.035);
    procLeg->AddEntry((TObject*)0, SampleConfig::DisplayNames(CurrentChannel())[SampleConfig::kZbb].c_str(),
                      "");

    for (int m = 0; m < nMassCuts; m++) {
        double x_pts[MAX_PT_BINS], y_pts[MAX_PT_BINS], x_err[MAX_PT_BINS], y_err[MAX_PT_BINS];
        for (int bin = 0; bin < CurrentNPtBins(); bin++) {
            x_pts[bin] = ptCenters[bin] + offsets[m];
            y_pts[bin] = data[m][bin];
            x_err[bin] = 0;
            y_err[bin] = std::abs(errs[m][bin]);
        }

        TGraphErrors* gr = new TGraphErrors(CurrentNPtBins(), x_pts, y_pts, x_err, y_err);
        gr->SetMarkerColor(colors[m]);
        gr->SetLineColor(colors[m]);
        gr->SetMarkerStyle(markers[m]);
        gr->SetMarkerSize(1.2);
        gr->SetLineWidth(2);
        gr->Draw("P same");

        legend->AddEntry(gr, Form("Mass Bin: %s", massCutLabels[m]), "p");
    }

    legend->Draw();
    procLeg->Draw();
    AtlasUtils::DrawAtlasLabel(0.18, 0.85);

    canvas->SaveAs((outputDir + saveName).c_str());

    delete frame;
    delete canvas;
    delete tick_latex;
    delete title_latex;
    delete refline;
    delete legend;
    delete procLeg;
}


// Same central values drawn twice per mass bin, with the uncorrelated errors
// and, 8 GeV to the right, the correlation-aware errors, to compare the two.
// The y-axis range is set to the maximum of the two sets of errors.
void MakeErrorComparisonPlot_vs_pT(const std::string& outputDir, const std::string& canvasName,
                                   const std::string& yTitle, const std::string& saveName,
                                   double values[MAX_MASS_CUTS][MAX_PT_BINS],
                                   double errs_uncorr[MAX_MASS_CUTS][MAX_PT_BINS],
                                   double errs_corr[MAX_MASS_CUTS][MAX_PT_BINS],
                                   double y_lo_force = std::numeric_limits<double>::quiet_NaN(),
                                   double y_hi_force = std::numeric_limits<double>::quiet_NaN()) {
    Int_t colors_mass[MAX_MASS_CUTS] = {kRed, kBlue, kGreen + 2};
    Int_t markers_mass[MAX_MASS_CUTS] = {20, 21, 22};

    double ptCenters[MAX_PT_BINS];
    for (int i = 0; i < CurrentNPtBins(); i++)
        ptCenters[i] = (CurrentPtBinsDisplay()[i] + CurrentPtBinsDisplay()[i + 1]) / 2.0;
    double offsets[MAX_MASS_CUTS] = {-20.0, 0.0, 20.0};

    // Y range from the uncorrelated errors, which are normally the larger ones.
    double ymin_val = 1e6, ymax_val = -1e6;
    for (int m = 0; m < nMassCuts; m++) {
        for (int bin = 0; bin < CurrentNPtBins(); bin++) {
            double lo = values[m][bin] - std::abs(errs_uncorr[m][bin]);
            double hi = values[m][bin] + std::abs(errs_uncorr[m][bin]);
            if (lo < ymin_val) ymin_val = lo;
            if (hi > ymax_val) ymax_val = hi;
        }
    }
    double range = ymax_val - ymin_val;
    double y_lo = ymin_val - range * RatioPadFactor();
    double y_hi = ymax_val + range * RatioPadFactor();
    if (y_lo > 0.95) y_lo = RatioYLoClamp();
    if (y_hi < 1.05) y_hi = RatioYHiClamp();

    TCanvas* ccomp = new TCanvas(canvasName.c_str(), "Error-treatment comparison", 900, 650);
    ccomp->SetLeftMargin(0.15);
    ccomp->SetBottomMargin(0.15);
    ccomp->SetTopMargin(0.12);
    ccomp->SetRightMargin(0.1);

    TH1F* frame = new TH1F((canvasName + "_frame").c_str(), "", CurrentNPtBins(), CurrentPtBinsDisplay());
    frame->SetMinimum(y_lo);
    frame->SetMaximum(y_hi);
    frame->GetXaxis()->SetTitle("p_{T} [GeV]");
    frame->GetYaxis()->SetTitle(yTitle.c_str());
    frame->GetXaxis()->SetTitleSize(0.04);
    frame->GetYaxis()->SetTitleSize(0.04);
    frame->GetXaxis()->SetLabelSize(0.04);
    frame->GetYaxis()->SetLabelSize(0.04);
    frame->GetYaxis()->SetTitleOffset(1.55);
    frame->SetStats(0);
    frame->GetXaxis()->SetRangeUser(CurrentPtBinsDisplay()[0], CurrentPtBinsDisplay()[CurrentNPtBins()]);
    frame->GetXaxis()->SetNdivisions(0, 0, 0, kFALSE);
    frame->GetXaxis()->SetLabelSize(0);
    frame->Draw("AXIS");

    double yrange = y_hi - y_lo;
    TLatex* tick_latex = new TLatex();
    tick_latex->SetTextSize(0.04);
    tick_latex->SetTextAlign(23);
    for (int i = 0; i <= CurrentNPtBins(); i++) {
        TLine* tick =
            new TLine(CurrentPtBinsDisplay()[i], y_lo, CurrentPtBinsDisplay()[i], y_lo + 0.01 * yrange);
        tick->SetLineColor(kBlack);
        tick->Draw("same");
        tick_latex->DrawLatex(CurrentPtBinsDisplay()[i], y_lo - 0.05 * yrange,
                              Form("%.0f", CurrentPtBinsDisplay()[i]));
    }
    TLatex* title_lat = new TLatex();
    title_lat->SetTextSize(0.045);
    title_lat->SetTextAlign(22);
    title_lat->SetNDC();
    title_lat->DrawLatex(0.85, 0.04, "p_{T} [GeV]");

    TLine* refline = new TLine(CurrentPtBinsDisplay()[0], 1.0, CurrentPtBinsDisplay()[CurrentNPtBins()], 1.0);
    refline->SetLineStyle(2);
    refline->SetLineColor(kBlack);
    refline->SetLineWidth(2);
    refline->Draw("same");

    TLegend* leg = new TLegend(0.15, 0.65, 0.60, 0.87);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.030);

    TLegend* procLeg = new TLegend(0.62, 0.75, 0.89, 0.89);
    procLeg->SetBorderSize(0);
    procLeg->SetFillStyle(0);
    procLeg->SetTextSize(0.035);
    procLeg->AddEntry((TObject*)0, SampleConfig::DisplayNames(CurrentChannel())[SampleConfig::kZbb].c_str(),
                      "");

    for (int m = 0; m < nMassCuts; m++) {
        double x_pts[MAX_PT_BINS], y_pts[MAX_PT_BINS], x_err[MAX_PT_BINS];
        double y_err_uncorr[MAX_PT_BINS], y_err_corr[MAX_PT_BINS];
        for (int bin = 0; bin < CurrentNPtBins(); bin++) {
            x_pts[bin] = ptCenters[bin] + offsets[m];
            y_pts[bin] = values[m][bin];
            x_err[bin] = 0;
            y_err_uncorr[bin] = std::abs(errs_uncorr[m][bin]);
            y_err_corr[bin] = std::abs(errs_corr[m][bin]);
        }

        // Uncorrelated: thin dashed error bars.
        TGraphErrors* gr_uncorr = new TGraphErrors(CurrentNPtBins(), x_pts, y_pts, x_err, y_err_uncorr);
        gr_uncorr->SetMarkerColor(colors_mass[m]);
        gr_uncorr->SetLineColor(colors_mass[m]);
        gr_uncorr->SetMarkerStyle(markers_mass[m]);
        gr_uncorr->SetMarkerSize(1.0);
        gr_uncorr->SetLineWidth(1);
        gr_uncorr->SetLineStyle(2);
        gr_uncorr->Draw("P same");

        // Corr-aware: open marker, shifted right so the bars do not overlap.
        double x_pts_shifted[MAX_PT_BINS];
        for (int bin = 0; bin < CurrentNPtBins(); bin++) x_pts_shifted[bin] = x_pts[bin] + 8.0;

        TGraphErrors* gr_corr = new TGraphErrors(CurrentNPtBins(), x_pts_shifted, y_pts, x_err, y_err_corr);
        gr_corr->SetMarkerColor(colors_mass[m]);
        gr_corr->SetLineColor(colors_mass[m]);
        gr_corr->SetMarkerStyle(markers_mass[m] + 4);
        gr_corr->SetMarkerSize(1.2);
        gr_corr->SetLineWidth(2);
        gr_corr->Draw("P same");

        leg->AddEntry(gr_uncorr, Form("%s - uncorr.", massCutLabels[m]), "p");
        leg->AddEntry(gr_corr, Form("%s - corr-aware", massCutLabels[m]), "p");
    }

    leg->Draw();
    procLeg->Draw();
    AtlasUtils::DrawAtlasLabel(0.18, 0.85);

    ccomp->SaveAs((outputDir + saveName).c_str());

    delete frame;
    delete ccomp;
    delete tick_latex;
    delete title_lat;
    delete refline;
    delete leg;
    delete procLeg;
}


// CSV summary, one row per (mass bin, pT bin) cell, written next to the text table. 
// Response_Ratios_Combined.h (ReadCSV) reads it to put DB and MJB on one 200-1200 GeV axis,
// and Response_Summary_Combined_Inclusive.h reads the inclusive-mass version.
//
// Warning: the readers use column positions, not the header names. 
// Add new columns only at the end, never in between.
//
// Columns, counted from 0:
//    0-7    channel, wp, variant, mass_bin (folder name), pt_bin_index,
//           pt_lo, pt_hi, pt_center (display edges, center = midpoint)
//    8-9    rjes_dr, rjes_dr_err                           (uncorrelated)
//   10-12   sig_insitu_dr, _err_unc, _err_corr
//   13-15   rjer_dr, _err_unc, _err_corr
//   16-23   single data/MC ratios: rjes_nom, _err, rjes_var, _err,
//           rjer_nom, _err, rjer_var, _err
//   24-31   b_JES inputs: mu_data and mu_combMC, nominal then variant (+ errors)
//   32-39   b_JER inputs: sigma_insitu and sigma^MC_insitu, nominal then variant
//   40-43   rjes_dr_err_corr, rjes_dr_err_full, sig_insitu_dr_err_full,
//           rjer_dr_err_full
//
// ReadCSV needs at least 16 columns, and uses 16-23, 24-39 and 40-43 only
// when the row is long enough. The "full" errors (41-43) use the correlations
// measured by the bootstrap (correlations.txt); if none is found, rho = 1.
inline void GenerateSummaryCSV(
    const std::string& outputFilePath, const std::string& channel, const std::string& wp,
    const std::string& variantLabel, double rjes_dr[MAX_MASS_CUTS][MAX_PT_BINS],
    double rjes_dr_errs[MAX_MASS_CUTS][MAX_PT_BINS], double sig_insitu_dr[MAX_MASS_CUTS][MAX_PT_BINS],
    double sig_insitu_dr_errs[MAX_MASS_CUTS][MAX_PT_BINS],
    double sig_insitu_dr_errs_corr[MAX_MASS_CUTS][MAX_PT_BINS], double rjer_dr[MAX_MASS_CUTS][MAX_PT_BINS],
    double rjer_dr_errs[MAX_MASS_CUTS][MAX_PT_BINS], double rjer_dr_errs_corr[MAX_MASS_CUTS][MAX_PT_BINS],
    // single data/MC ratios (columns 16-23)
    double rjes_nom[MAX_MASS_CUTS][MAX_PT_BINS], double rjes_nom_err[MAX_MASS_CUTS][MAX_PT_BINS],
    double rjes_var[MAX_MASS_CUTS][MAX_PT_BINS], double rjes_var_err[MAX_MASS_CUTS][MAX_PT_BINS],
    double rjer_nom[MAX_MASS_CUTS][MAX_PT_BINS], double rjer_nom_err[MAX_MASS_CUTS][MAX_PT_BINS],
    double rjer_var[MAX_MASS_CUTS][MAX_PT_BINS], double rjer_var_err[MAX_MASS_CUTS][MAX_PT_BINS],
    // b_JES / b_JER data and MC values (columns 24-39)
    double bjes_data_nom[MAX_MASS_CUTS][MAX_PT_BINS], double bjes_data_nom_err[MAX_MASS_CUTS][MAX_PT_BINS],
    double bjes_mc_nom[MAX_MASS_CUTS][MAX_PT_BINS], double bjes_mc_nom_err[MAX_MASS_CUTS][MAX_PT_BINS],
    double bjes_data_var[MAX_MASS_CUTS][MAX_PT_BINS], double bjes_data_var_err[MAX_MASS_CUTS][MAX_PT_BINS],
    double bjes_mc_var[MAX_MASS_CUTS][MAX_PT_BINS], double bjes_mc_var_err[MAX_MASS_CUTS][MAX_PT_BINS],
    double bjer_data_nom[MAX_MASS_CUTS][MAX_PT_BINS], double bjer_data_nom_err[MAX_MASS_CUTS][MAX_PT_BINS],
    double bjer_mc_nom[MAX_MASS_CUTS][MAX_PT_BINS], double bjer_mc_nom_err[MAX_MASS_CUTS][MAX_PT_BINS],
    double bjer_data_var[MAX_MASS_CUTS][MAX_PT_BINS], double bjer_data_var_err[MAX_MASS_CUTS][MAX_PT_BINS],
    double bjer_mc_var[MAX_MASS_CUTS][MAX_PT_BINS], double bjer_mc_var_err[MAX_MASS_CUTS][MAX_PT_BINS],
    // correlated double-ratio errors (columns 40-43)
    double rjes_dr_errs_corr[MAX_MASS_CUTS][MAX_PT_BINS],
    double rjes_dr_errs_full[MAX_MASS_CUTS][MAX_PT_BINS],
    double sig_insitu_dr_errs_full[MAX_MASS_CUTS][MAX_PT_BINS],
    double rjer_dr_errs_full[MAX_MASS_CUTS][MAX_PT_BINS]) {
    std::ofstream out(outputFilePath.c_str());
    if (!out.is_open()) {
        std::cerr << "[Ratio_Response] WARN: cannot open CSV summary at " << outputFilePath << "\n";
        return;
    }
    // Header line; the readers skip it.
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
        // columns 40-43
        << "rjes_dr_err_corr,rjes_dr_err_full,sig_insitu_dr_err_full,rjer_dr_err_full\n";

    // Display edges of pT bins, not midpoints.
    out << std::setprecision(8);
    for (int m = 0; m < nMassCuts; ++m) {
        for (int b = 0; b < CurrentNPtBins(); ++b) {
            double lo = CurrentPtBinsDisplay()[b];
            double hi = CurrentPtBinsDisplay()[b + 1];
            double mid = 0.5 * (lo + hi);
            // The folder name ("50-80GeV", ...) has no spaces or LaTeX, and the
            // readers map it back to a display label (MassTag()).
            std::string mlabel = massCutDirNames[m];
            // A comma would shift every later column.
            for (auto& ch : mlabel)
                if (ch == ',') ch = ';';
            out << channel << "," << wp << "," << variantLabel << "," << mlabel << "," << b << "," << lo
                << "," << hi << "," << mid << "," << rjes_dr[m][b] << "," << rjes_dr_errs[m][b] << ","
                << sig_insitu_dr[m][b] << "," << sig_insitu_dr_errs[m][b] << ","
                << sig_insitu_dr_errs_corr[m][b] << "," << rjer_dr[m][b] << "," << rjer_dr_errs[m][b] << ","
                << rjer_dr_errs_corr[m][b] << "," << rjes_nom[m][b] << "," << rjes_nom_err[m][b] << ","
                << rjes_var[m][b] << "," << rjes_var_err[m][b] << "," << rjer_nom[m][b] << ","
                << rjer_nom_err[m][b] << "," << rjer_var[m][b] << "," << rjer_var_err[m][b] << ","
                << bjes_data_nom[m][b] << "," << bjes_data_nom_err[m][b] << "," << bjes_mc_nom[m][b] << ","
                << bjes_mc_nom_err[m][b] << "," << bjes_data_var[m][b] << "," << bjes_data_var_err[m][b]
                << "," << bjes_mc_var[m][b] << "," << bjes_mc_var_err[m][b] << "," << bjer_data_nom[m][b]
                << "," << bjer_data_nom_err[m][b] << "," << bjer_mc_nom[m][b] << "," << bjer_mc_nom_err[m][b]
                << "," << bjer_data_var[m][b] << "," << bjer_data_var_err[m][b] << "," << bjer_mc_var[m][b]
                << "," << bjer_mc_var_err[m][b] << "," << rjes_dr_errs_corr[m][b] << ","
                << rjes_dr_errs_full[m][b] << "," << sig_insitu_dr_errs_full[m][b] << ","
                << rjer_dr_errs_full[m][b] << "\n";
        }
    }
    out.close();
    std::cout << "[Ratio_Response] CSV summary saved: " << outputFilePath << "\n";
}


// Readable summary table (Summary_Ratios_<wp>.txt).
void GenerateSummaryTable(
    const std::string& outputFilePath, const std::string& labelWP_display, const std::string& variantLabel,
    double ratio_data_means[MAX_MASS_CUTS][MAX_PT_BINS],
    double ratio_data_mean_errs[MAX_MASS_CUTS][MAX_PT_BINS],
    double ratio_mc_means[MAX_MASS_CUTS][MAX_PT_BINS], double ratio_mc_mean_errs[MAX_MASS_CUTS][MAX_PT_BINS],
    double ratio_bkg_means[MAX_MASS_CUTS][MAX_PT_BINS],
    double ratio_bkg_mean_errs[MAX_MASS_CUTS][MAX_PT_BINS],
    double ratio_combined_means[MAX_MASS_CUTS][MAX_PT_BINS],
    double ratio_combined_mean_errs[MAX_MASS_CUTS][MAX_PT_BINS],
    double ratio_data_sigmas[MAX_MASS_CUTS][MAX_PT_BINS],
    double ratio_data_sigma_errs[MAX_MASS_CUTS][MAX_PT_BINS],
    double ratio_mc_sigmas[MAX_MASS_CUTS][MAX_PT_BINS],
    double ratio_mc_sigma_errs[MAX_MASS_CUTS][MAX_PT_BINS],
    double ratio_bkg_sigmas[MAX_MASS_CUTS][MAX_PT_BINS],
    double ratio_bkg_sigma_errs[MAX_MASS_CUTS][MAX_PT_BINS],
    double ratio_combined_sigmas[MAX_MASS_CUTS][MAX_PT_BINS],
    double ratio_combined_sigma_errs[MAX_MASS_CUTS][MAX_PT_BINS], double rjes_dr[MAX_MASS_CUTS][MAX_PT_BINS],
    double rjes_dr_errs[MAX_MASS_CUTS][MAX_PT_BINS], double sig_insitu_dr[MAX_MASS_CUTS][MAX_PT_BINS],
    double sig_insitu_dr_errs[MAX_MASS_CUTS][MAX_PT_BINS],
    double sig_insitu_dr_errs_corr[MAX_MASS_CUTS][MAX_PT_BINS], double rjer_dr[MAX_MASS_CUTS][MAX_PT_BINS],
    double rjer_dr_errs[MAX_MASS_CUTS][MAX_PT_BINS], double rjer_dr_errs_corr[MAX_MASS_CUTS][MAX_PT_BINS],
    double R_JER_var[MAX_MASS_CUTS][MAX_PT_BINS], double R_JER_var_err[MAX_MASS_CUTS][MAX_PT_BINS],
    double R_JER_nom[MAX_MASS_CUTS][MAX_PT_BINS], double R_JER_nom_err[MAX_MASS_CUTS][MAX_PT_BINS]) {
    std::ofstream out(outputFilePath.c_str());
    if (!out.is_open()) {
        std::cerr << "Error: Could not open summary file: " << outputFilePath << std::endl;
        return;
    }
    // "R_DB" or "R_MJB", only used in the text.
    const std::string nom_branch_R = GetVariant(CurrentChannel()).branch_R;

    out << "========================================================================================================================\n";
    out << "Summary of " << nom_branch_R << " Ratios (" << variantLabel
        << " / Nominal) for Working Point: " << labelWP_display << "\n";
    out << "========================================================================================================================\n\n";

    out << "Definitions:\n";
    out << "  R_JES     = <" << nom_branch_R << ">_data / <" << nom_branch_R << ">_combMC\n";
    out << "  sigma_insitu    = sqrt(sigma^2_data - sigma^2_truth)              (data-side)\n";
    out << "  sigma^MC_insitu = sqrt(sigma^2_combMC - sigma^2_truth)\n";
    out << "  R_JER     = sigma_insitu / sigma^MC_insitu                       (data/MC JER ratio)\n";
    out << "  DR(X)     = X_variant / X_nominal\n\n";

    out << std::fixed << std::setprecision(4);

    // Per-sample mean ratios
    out << "------------------------------------------------------------------------------------------------------------------------\n";
    out << "Per-sample Mean Ratios (<" << nom_branch_R << "^{After}> / <" << nom_branch_R << "^{Before}>)\n";
    out << "------------------------------------------------------------------------------------------------------------------------\n";
    out << std::setw(15) << "Mass Bin" << std::setw(25) << "pT Bin" << std::setw(20) << "Data Ratio"
        << std::setw(20) << "MC Ratio" << std::setw(20) << "Bkg Ratio" << std::setw(20) << "Comb. MC Ratio"
        << std::setw(25) << "DR(R_JES)" << "\n";
    out << "------------------------------------------------------------------------------------------------------------------------\n";
    for (int m = 0; m < nMassCuts; m++) {
        for (int bin = 0; bin < CurrentNPtBins(); bin++) {
            out << std::setw(15) << massCutLabels[m] << std::setw(25) << CurrentPtBinLabel(bin)
                << std::setw(20)
                << Form("%.4f +/- %.4f", ratio_data_means[m][bin], ratio_data_mean_errs[m][bin])
                << std::setw(20) << Form("%.4f +/- %.4f", ratio_mc_means[m][bin], ratio_mc_mean_errs[m][bin])
                << std::setw(20) << Form("%.4f +/- %.4f", ratio_bkg_means[m][bin], ratio_bkg_mean_errs[m][bin])
                << std::setw(20)
                << Form("%.4f +/- %.4f", ratio_combined_means[m][bin], ratio_combined_mean_errs[m][bin])
                << std::setw(25) << Form("%.4f +/- %.4f", rjes_dr[m][bin], rjes_dr_errs[m][bin]) << "\n";
        }
    }
    out << "------------------------------------------------------------------------------------------------------------------------\n\n";

    // Per-sample sigma ratios
    out << "------------------------------------------------------------------------------------------------------------------------\n";
    out << "Per-sample Resolution Ratios (sigma^{After} / sigma^{Before})\n";
    out << "------------------------------------------------------------------------------------------------------------------------\n";
    out << std::setw(15) << "Mass Bin" << std::setw(25) << "pT Bin" << std::setw(20) << "Data Ratio"
        << std::setw(20) << "MC Ratio" << std::setw(20) << "Bkg Ratio" << std::setw(20) << "Comb. MC Ratio"
        << "\n";
    out << "------------------------------------------------------------------------------------------------------------------------\n";
    for (int m = 0; m < nMassCuts; m++) {
        for (int bin = 0; bin < CurrentNPtBins(); bin++) {
            out << std::setw(15) << massCutLabels[m] << std::setw(25) << CurrentPtBinLabel(bin)
                << std::setw(20)
                << Form("%.4f +/- %.4f", ratio_data_sigmas[m][bin], ratio_data_sigma_errs[m][bin])
                << std::setw(20) << Form("%.4f +/- %.4f", ratio_mc_sigmas[m][bin], ratio_mc_sigma_errs[m][bin])
                << std::setw(20)
                << Form("%.4f +/- %.4f", ratio_bkg_sigmas[m][bin], ratio_bkg_sigma_errs[m][bin])
                << std::setw(20)
                << Form("%.4f +/- %.4f", ratio_combined_sigmas[m][bin], ratio_combined_sigma_errs[m][bin])
                << "\n";
        }
    }
    out << "------------------------------------------------------------------------------------------------------------------------\n\n";

    // Double ratios of sigma_insitu and R_JER
    out << "------------------------------------------------------------------------------------------------------------------------\n";
    out << "Physics Double Ratios vs Nominal (uncorrelated / correlation-aware errors)\n";
    out << "------------------------------------------------------------------------------------------------------------------------\n";
    out << std::setw(15) << "Mass Bin" << std::setw(25) << "pT Bin" << std::setw(26) << "DR(sigma_insitu) unc"
        << std::setw(26) << "DR(sigma_insitu) cor" << std::setw(22) << "DR(R_JER) unc" << std::setw(22)
        << "DR(R_JER) cor"
        << "\n";
    out << "------------------------------------------------------------------------------------------------------------------------\n";
    for (int m = 0; m < nMassCuts; m++) {
        for (int bin = 0; bin < CurrentNPtBins(); bin++) {
            out << std::setw(15) << massCutLabels[m] << std::setw(25) << CurrentPtBinLabel(bin)
                << std::setw(26) << Form("%.4f +/- %.4f", sig_insitu_dr[m][bin], sig_insitu_dr_errs[m][bin])
                << std::setw(26)
                << Form("%.4f +/- %.4f", sig_insitu_dr[m][bin], sig_insitu_dr_errs_corr[m][bin])
                << std::setw(22) << Form("%.4f +/- %.4f", rjer_dr[m][bin], rjer_dr_errs[m][bin])
                << std::setw(22) << Form("%.4f +/- %.4f", rjer_dr[m][bin], rjer_dr_errs_corr[m][bin]) << "\n";
        }
    }
    out << "------------------------------------------------------------------------------------------------------------------------\n\n";

    // R_JER of the variant and of the nominal
    out << "------------------------------------------------------------------------------------------------------------------------\n";
    out << "R_JER (data/MC) values per cell\n";
    out << "------------------------------------------------------------------------------------------------------------------------\n";
    out << std::setw(15) << "Mass Bin" << std::setw(25) << "pT Bin" << std::setw(26) << "R_JER ("
        << variantLabel << ")" << std::setw(22) << "R_JER (Nominal)"
        << "\n";
    out << "------------------------------------------------------------------------------------------------------------------------\n";
    for (int m = 0; m < nMassCuts; m++) {
        for (int bin = 0; bin < CurrentNPtBins(); bin++) {
            out << std::setw(15) << massCutLabels[m] << std::setw(25) << CurrentPtBinLabel(bin)
                << std::setw(26) << Form("%.4f +/- %.4f", R_JER_var[m][bin], R_JER_var_err[m][bin])
                << std::setw(22) << Form("%.4f +/- %.4f", R_JER_nom[m][bin], R_JER_nom_err[m][bin]) << "\n";
        }
    }
    out << "------------------------------------------------------------------------------------------------------------------------\n";

    out.close();
    std::cout << "Summary table saved: " << outputFilePath << std::endl;
}


// One variant against the nominal, for one working point.
// "before" = nominal, "after" = variant throughout.
void process_single_wp_ratio(const std::string& variantKey, const std::string& wp_str) {

    VariantConfig cfg_nominal = GetVariant(CurrentChannel());
    VariantConfig cfg_variant = GetVariant(variantKey);

    std::string fileWP = PathBuilder::WPtoFileWP(wp_str);
    std::string labelWP = PathBuilder::WPtoLabel(wp_str);
    std::string variantLabel = cfg_variant.key;

    // Separate folder for the inclusive-mass run, so it does not overwrite the
    // 3-mass-bin CSV and plots.
    std::string r2dir = MassInclusiveMode() ? "R2pTbin_MassInclusive" : "R2pTbin";
    std::string outputDir = (RunConfig::PlotDir() + "/Response/Ratios/") + CurrentChannel() + "_Ratio_"
                            + variantLabel + "/" + fileWP + "/" + r2dir + "/";
    PathBuilder::MakeDir(outputDir);

    std::cout << "\n============================================================" << std::endl;
    std::cout << "Ratio Response: " << variantLabel << " / Nominal - WP: " << labelWP << std::endl;
    std::cout << "Reading fit parameters from Rdb_pT_bins_MCut output..." << std::endl;
    std::cout << "============================================================\n" << std::endl;

    // 1. Read the fit results of the nominal ("before") and the variant ("after").
    // Missing files leave zeros, which give zero ratios downstream.
    double data_means_before[MAX_MASS_CUTS][MAX_PT_BINS] = {},
           data_mean_errs_before[MAX_MASS_CUTS][MAX_PT_BINS] = {};
    double mc_means_before[MAX_MASS_CUTS][MAX_PT_BINS] = {},
           mc_mean_errs_before[MAX_MASS_CUTS][MAX_PT_BINS] = {};
    double bkg_means_before[MAX_MASS_CUTS][MAX_PT_BINS] = {},
           bkg_mean_errs_before[MAX_MASS_CUTS][MAX_PT_BINS] = {};
    double comb_means_before[MAX_MASS_CUTS][MAX_PT_BINS] = {},
           comb_mean_errs_before[MAX_MASS_CUTS][MAX_PT_BINS] = {};
    double data_sigmas_before[MAX_MASS_CUTS][MAX_PT_BINS] = {},
           data_sigma_errs_before[MAX_MASS_CUTS][MAX_PT_BINS] = {};
    double mc_sigmas_before[MAX_MASS_CUTS][MAX_PT_BINS] = {},
           mc_sigma_errs_before[MAX_MASS_CUTS][MAX_PT_BINS] = {};
    double bkg_sigmas_before[MAX_MASS_CUTS][MAX_PT_BINS] = {},
           bkg_sigma_errs_before[MAX_MASS_CUTS][MAX_PT_BINS] = {};
    double comb_sigmas_before[MAX_MASS_CUTS][MAX_PT_BINS] = {},
           comb_sigma_errs_before[MAX_MASS_CUTS][MAX_PT_BINS] = {};

    double data_means_after[MAX_MASS_CUTS][MAX_PT_BINS] = {},
           data_mean_errs_after[MAX_MASS_CUTS][MAX_PT_BINS] = {};
    double mc_means_after[MAX_MASS_CUTS][MAX_PT_BINS] = {},
           mc_mean_errs_after[MAX_MASS_CUTS][MAX_PT_BINS] = {};
    double bkg_means_after[MAX_MASS_CUTS][MAX_PT_BINS] = {},
           bkg_mean_errs_after[MAX_MASS_CUTS][MAX_PT_BINS] = {};
    double comb_means_after[MAX_MASS_CUTS][MAX_PT_BINS] = {},
           comb_mean_errs_after[MAX_MASS_CUTS][MAX_PT_BINS] = {};
    double data_sigmas_after[MAX_MASS_CUTS][MAX_PT_BINS] = {},
           data_sigma_errs_after[MAX_MASS_CUTS][MAX_PT_BINS] = {};
    double mc_sigmas_after[MAX_MASS_CUTS][MAX_PT_BINS] = {},
           mc_sigma_errs_after[MAX_MASS_CUTS][MAX_PT_BINS] = {};
    double bkg_sigmas_after[MAX_MASS_CUTS][MAX_PT_BINS] = {},
           bkg_sigma_errs_after[MAX_MASS_CUTS][MAX_PT_BINS] = {};
    double comb_sigmas_after[MAX_MASS_CUTS][MAX_PT_BINS] = {},
           comb_sigma_errs_after[MAX_MASS_CUTS][MAX_PT_BINS] = {};

    double truth_sigmas_before[MAX_MASS_CUTS][MAX_PT_BINS] = {},
           truth_sigma_errs_before[MAX_MASS_CUTS][MAX_PT_BINS] = {};
    double truth_sigmas_after[MAX_MASS_CUTS][MAX_PT_BINS] = {},
           truth_sigma_errs_after[MAX_MASS_CUTS][MAX_PT_BINS] = {};

    const char* fitTypes[4] = {"data", "mc", "bkg", "combined_MC"};

    int filesRead = 0, filesMissing = 0;

    for (int m = 0; m < nMassCuts; m++) {
        for (int bin = 0; bin < CurrentNPtBins(); bin++) {
            for (int t = 0; t < 4; t++) {
                std::string path_before = Form(((RunConfig::PlotDir() + "/%s/Response_Fit_pTbins_Mcuts_")
                                                + CurrentChannel() + "/%s/%s/%s_fit_parameters_ptbin%d.txt")
                                                   .c_str(),
                                               cfg_nominal.plot_subdir.c_str(), fileWP.c_str(),
                                               massCutDirNames[m], fitTypes[t], bin);
                std::string path_after = Form(((RunConfig::PlotDir() + "/%s/Response_Fit_pTbins_Mcuts_")
                                               + CurrentChannel() + "/%s/%s/%s_fit_parameters_ptbin%d.txt")
                                                  .c_str(),
                                              cfg_variant.plot_subdir.c_str(), fileWP.c_str(),
                                              massCutDirNames[m], fitTypes[t], bin);

                FitParams fp_before, fp_after;
                bool ok_before = ReadFitParams(path_before, fp_before);
                bool ok_after = ReadFitParams(path_after, fp_after);

                if (ok_before)
                    filesRead++;
                else
                    filesMissing++;
                if (ok_after)
                    filesRead++;
                else
                    filesMissing++;

                double *mean_b, *mean_e_b, *sig_b, *sig_e_b;
                double *mean_a, *mean_e_a, *sig_a, *sig_e_a;

                if (t == 0) {
                    mean_b = &data_means_before[m][bin];
                    mean_e_b = &data_mean_errs_before[m][bin];
                    sig_b = &data_sigmas_before[m][bin];
                    sig_e_b = &data_sigma_errs_before[m][bin];
                    mean_a = &data_means_after[m][bin];
                    mean_e_a = &data_mean_errs_after[m][bin];
                    sig_a = &data_sigmas_after[m][bin];
                    sig_e_a = &data_sigma_errs_after[m][bin];
                } else if (t == 1) {
                    mean_b = &mc_means_before[m][bin];
                    mean_e_b = &mc_mean_errs_before[m][bin];
                    sig_b = &mc_sigmas_before[m][bin];
                    sig_e_b = &mc_sigma_errs_before[m][bin];
                    mean_a = &mc_means_after[m][bin];
                    mean_e_a = &mc_mean_errs_after[m][bin];
                    sig_a = &mc_sigmas_after[m][bin];
                    sig_e_a = &mc_sigma_errs_after[m][bin];
                } else if (t == 2) {
                    mean_b = &bkg_means_before[m][bin];
                    mean_e_b = &bkg_mean_errs_before[m][bin];
                    sig_b = &bkg_sigmas_before[m][bin];
                    sig_e_b = &bkg_sigma_errs_before[m][bin];
                    mean_a = &bkg_means_after[m][bin];
                    mean_e_a = &bkg_mean_errs_after[m][bin];
                    sig_a = &bkg_sigmas_after[m][bin];
                    sig_e_a = &bkg_sigma_errs_after[m][bin];
                } else {
                    mean_b = &comb_means_before[m][bin];
                    mean_e_b = &comb_mean_errs_before[m][bin];
                    sig_b = &comb_sigmas_before[m][bin];
                    sig_e_b = &comb_sigma_errs_before[m][bin];
                    mean_a = &comb_means_after[m][bin];
                    mean_e_a = &comb_mean_errs_after[m][bin];
                    sig_a = &comb_sigmas_after[m][bin];
                    sig_e_a = &comb_sigma_errs_after[m][bin];
                }

                *mean_b = fp_before.mean;
                *mean_e_b = fp_before.mean_err;
                *sig_b = fp_before.sigma;
                *sig_e_b = fp_before.sigma_err;
                *mean_a = fp_after.mean;
                *mean_e_a = fp_after.mean_err;
                *sig_a = fp_after.sigma;
                *sig_e_a = fp_after.sigma_err;
            }
        }
    }

    std::cout << "Fit parameter files: " << filesRead << " read, " << filesMissing << " missing."
              << std::endl;

    // Truth fits (sigma_truth)
    int truthFilesRead = 0, truthFilesMissing = 0;
    for (int m = 0; m < nMassCuts; m++) {
        for (int bin = 0; bin < CurrentNPtBins(); bin++) {
            std::string path_truth_before =
                Form(((RunConfig::PlotDir() + "/%s/Response_Fit_pTbins_Mcuts_") + CurrentChannel()
                      + "/%s/%s/truth_ref_fit_parameters_ptbin%d.txt")
                         .c_str(),
                     cfg_nominal.plot_subdir.c_str(), fileWP.c_str(), massCutDirNames[m], bin);
            std::string path_truth_after =
                Form(((RunConfig::PlotDir() + "/%s/Response_Fit_pTbins_Mcuts_") + CurrentChannel()
                      + "/%s/%s/truth_ref_fit_parameters_ptbin%d.txt")
                         .c_str(),
                     cfg_variant.plot_subdir.c_str(), fileWP.c_str(), massCutDirNames[m], bin);

            FitParams fp_truth_before, fp_truth_after;
            if (ReadFitParams(path_truth_before, fp_truth_before)) {
                truth_sigmas_before[m][bin] = fp_truth_before.sigma;
                truth_sigma_errs_before[m][bin] = fp_truth_before.sigma_err;
                truthFilesRead++;
            } else {
                truthFilesMissing++;
            }
            if (ReadFitParams(path_truth_after, fp_truth_after)) {
                truth_sigmas_after[m][bin] = fp_truth_after.sigma;
                truth_sigma_errs_after[m][bin] = fp_truth_after.sigma_err;
                truthFilesRead++;
            } else {
                truthFilesMissing++;
            }
        }
    }
    std::cout << "Truth fit files: " << truthFilesRead << " read, " << truthFilesMissing << " missing."
              << std::endl;
    if (truthFilesMissing > 0) {
        std::cerr
            << "WARNING: Missing truth fit files - JER double ratios will use sigma_truth=0 for missing cells.\n"
            << "Run Rdb_pT_bins_MCut_" << CurrentChannel() << " with truth branches to produce these files.\n"
            << std::endl;
    }
    if (filesMissing > 0) {
        std::cerr << "\nWARNING: Missing fit files will produce zero ratios!" << std::endl;
        std::cerr << "Make sure Rdb_pT_bins_MCut_" << CurrentChannel() << ".cxx was run for BOTH \""
                  << CurrentChannel() << "\" and \"" << variantKey << "\".\n"
                  << std::endl;
    }

    // 2. Ratios variant / nominal (after/before correction)

    // Per-sample ratios variant / nominal, all with uncorrelated errors.
    double ratio_data_means[MAX_MASS_CUTS][MAX_PT_BINS] = {},
           ratio_data_mean_errs[MAX_MASS_CUTS][MAX_PT_BINS] = {};
    double ratio_mc_means[MAX_MASS_CUTS][MAX_PT_BINS] = {},
           ratio_mc_mean_errs[MAX_MASS_CUTS][MAX_PT_BINS] = {};
    double ratio_bkg_means[MAX_MASS_CUTS][MAX_PT_BINS] = {},
           ratio_bkg_mean_errs[MAX_MASS_CUTS][MAX_PT_BINS] = {};
    double ratio_comb_means[MAX_MASS_CUTS][MAX_PT_BINS] = {},
           ratio_comb_mean_errs[MAX_MASS_CUTS][MAX_PT_BINS] = {};
    double ratio_data_sigmas[MAX_MASS_CUTS][MAX_PT_BINS] = {},
           ratio_data_sigma_errs[MAX_MASS_CUTS][MAX_PT_BINS] = {};
    double ratio_mc_sigmas[MAX_MASS_CUTS][MAX_PT_BINS] = {},
           ratio_mc_sigma_errs[MAX_MASS_CUTS][MAX_PT_BINS] = {};
    double ratio_bkg_sigmas[MAX_MASS_CUTS][MAX_PT_BINS] = {},
           ratio_bkg_sigma_errs[MAX_MASS_CUTS][MAX_PT_BINS] = {};
    double ratio_comb_sigmas[MAX_MASS_CUTS][MAX_PT_BINS] = {},
           ratio_comb_sigma_errs[MAX_MASS_CUTS][MAX_PT_BINS] = {};

    // DR(R_JES) = (mu_data / mu_combMC)_var / (mu_data / mu_combMC)_nom
    double rjes_dr[MAX_MASS_CUTS][MAX_PT_BINS] = {}, rjes_dr_errs[MAX_MASS_CUTS][MAX_PT_BINS] = {};
    // corr-aware (rho = 1) and full (4x4 rho from the "resp" lines of the correlation file)
    double rjes_dr_errs_corr[MAX_MASS_CUTS][MAX_PT_BINS] = {};
    double rjes_dr_errs_full[MAX_MASS_CUTS][MAX_PT_BINS] = {};
    double sig_insitu_dr_errs_full[MAX_MASS_CUTS][MAX_PT_BINS] = {};
    double rjer_dr_errs_full[MAX_MASS_CUTS][MAX_PT_BINS] = {};

    // DR(sigma_insitu) = sigma_insitu_var / sigma_insitu_nom
    double sig_insitu_dr[MAX_MASS_CUTS][MAX_PT_BINS] = {},
           sig_insitu_dr_errs[MAX_MASS_CUTS][MAX_PT_BINS] = {};
    double sig_insitu_dr_errs_corr[MAX_MASS_CUTS][MAX_PT_BINS] = {};

    // DR(R_JER) = R_JER_var / R_JER_nom, R_JER = sigma_insitu / sigma^MC_insitu
    double rjer_dr[MAX_MASS_CUTS][MAX_PT_BINS] = {}, rjer_dr_errs[MAX_MASS_CUTS][MAX_PT_BINS] = {};
    double rjer_dr_errs_corr[MAX_MASS_CUTS][MAX_PT_BINS] = {};

    // R_JER per cell, for the text table and the CSV
    double R_JER_var_arr[MAX_MASS_CUTS][MAX_PT_BINS] = {}, R_JER_var_err_arr[MAX_MASS_CUTS][MAX_PT_BINS] = {};
    double R_JER_nom_arr[MAX_MASS_CUTS][MAX_PT_BINS] = {}, R_JER_nom_err_arr[MAX_MASS_CUTS][MAX_PT_BINS] = {};

    // Single R_JES = mu_data / mu_combMC per cell, nominal and variant, for the
    // data/MC plots of Response_Ratios_Combined.
    double R_JES_var_arr[MAX_MASS_CUTS][MAX_PT_BINS] = {}, R_JES_var_err_arr[MAX_MASS_CUTS][MAX_PT_BINS] = {};
    double R_JES_nom_arr[MAX_MASS_CUTS][MAX_PT_BINS] = {}, R_JES_nom_err_arr[MAX_MASS_CUTS][MAX_PT_BINS] = {};

    // b_JER data = sigma_insitu and MC = sigma^MC_insitu (numerator and denominator of R_JER), 
    // nominal and variant. The b_JES data/MC values are the data_means_*/comb_means_* arrays themselves.
    double bJER_data_nom_arr[MAX_MASS_CUTS][MAX_PT_BINS] = {},
           bJER_data_nom_err_arr[MAX_MASS_CUTS][MAX_PT_BINS] = {};
    double bJER_mc_nom_arr[MAX_MASS_CUTS][MAX_PT_BINS] = {},
           bJER_mc_nom_err_arr[MAX_MASS_CUTS][MAX_PT_BINS] = {};
    double bJER_data_var_arr[MAX_MASS_CUTS][MAX_PT_BINS] = {},
           bJER_data_var_err_arr[MAX_MASS_CUTS][MAX_PT_BINS] = {};
    double bJER_mc_var_arr[MAX_MASS_CUTS][MAX_PT_BINS] = {},
           bJER_mc_var_err_arr[MAX_MASS_CUTS][MAX_PT_BINS] = {};

    std::cout << "\nCalculating ratios (" << variantLabel << " / Nominal)...\n" << std::endl;

    for (int m = 0; m < nMassCuts; m++) {
        for (int bin = 0; bin < CurrentNPtBins(); bin++) {
            // Per-sample mean ratios
            ratio_data_means[m][bin] = SafeRatio(data_means_after[m][bin], data_means_before[m][bin]);
            ratio_data_mean_errs[m][bin] =
                SafeRatioErr(data_means_after[m][bin], data_mean_errs_after[m][bin],
                             data_means_before[m][bin], data_mean_errs_before[m][bin]);

            ratio_mc_means[m][bin] = SafeRatio(mc_means_after[m][bin], mc_means_before[m][bin]);
            ratio_mc_mean_errs[m][bin] = SafeRatioErr(mc_means_after[m][bin], mc_mean_errs_after[m][bin],
                                                      mc_means_before[m][bin], mc_mean_errs_before[m][bin]);

            ratio_bkg_means[m][bin] = SafeRatio(bkg_means_after[m][bin], bkg_means_before[m][bin]);
            ratio_bkg_mean_errs[m][bin] =
                SafeRatioErr(bkg_means_after[m][bin], bkg_mean_errs_after[m][bin], bkg_means_before[m][bin],
                             bkg_mean_errs_before[m][bin]);

            ratio_comb_means[m][bin] = SafeRatio(comb_means_after[m][bin], comb_means_before[m][bin]);
            ratio_comb_mean_errs[m][bin] =
                SafeRatioErr(comb_means_after[m][bin], comb_mean_errs_after[m][bin],
                             comb_means_before[m][bin], comb_mean_errs_before[m][bin]);

            // Per-sample sigma ratios
            ratio_mc_sigmas[m][bin] = SafeRatio(mc_sigmas_after[m][bin], mc_sigmas_before[m][bin]);
            ratio_mc_sigma_errs[m][bin] =
                SafeRatioErr(mc_sigmas_after[m][bin], mc_sigma_errs_after[m][bin], mc_sigmas_before[m][bin],
                             mc_sigma_errs_before[m][bin]);

            ratio_bkg_sigmas[m][bin] = SafeRatio(bkg_sigmas_after[m][bin], bkg_sigmas_before[m][bin]);
            ratio_bkg_sigma_errs[m][bin] =
                SafeRatioErr(bkg_sigmas_after[m][bin], bkg_sigma_errs_after[m][bin],
                             bkg_sigmas_before[m][bin], bkg_sigma_errs_before[m][bin]);

            ratio_comb_sigmas[m][bin] = SafeRatio(comb_sigmas_after[m][bin], comb_sigmas_before[m][bin]);
            ratio_comb_sigma_errs[m][bin] =
                SafeRatioErr(comb_sigmas_after[m][bin], comb_sigma_errs_after[m][bin],
                             comb_sigmas_before[m][bin], comb_sigma_errs_before[m][bin]);

            // Raw data sigma ratio, without the truth subtraction (that is DR(sigma_insitu)).
            ratio_data_sigmas[m][bin] = SafeRatio(data_sigmas_after[m][bin], data_sigmas_before[m][bin]);
            ratio_data_sigma_errs[m][bin] =
                SafeRatioErr(data_sigmas_after[m][bin], data_sigma_errs_after[m][bin],
                             data_sigmas_before[m][bin], data_sigma_errs_before[m][bin]);

            // DR(R_JES), uncorrelated error
            double d_over_c_before = SafeRatio(data_means_before[m][bin], comb_means_before[m][bin]);
            double d_over_c_before_err =
                SafeRatioErr(data_means_before[m][bin], data_mean_errs_before[m][bin],
                             comb_means_before[m][bin], comb_mean_errs_before[m][bin]);
            double d_over_c_after = SafeRatio(data_means_after[m][bin], comb_means_after[m][bin]);
            double d_over_c_after_err = SafeRatioErr(data_means_after[m][bin], data_mean_errs_after[m][bin],
                                                     comb_means_after[m][bin], comb_mean_errs_after[m][bin]);

            rjes_dr[m][bin] = SafeRatio(d_over_c_after, d_over_c_before);
            rjes_dr_errs[m][bin] =
                SafeRatioErr(d_over_c_after, d_over_c_after_err, d_over_c_before, d_over_c_before_err);

            // Single R_JES, nominal and variant.
            R_JES_nom_arr[m][bin] = d_over_c_before;
            R_JES_nom_err_arr[m][bin] = d_over_c_before_err;
            R_JES_var_arr[m][bin] = d_over_c_after;
            R_JES_var_err_arr[m][bin] = d_over_c_after_err;

            // DR(R_JES), corr-aware: rho = 1 between variant and nominal of the same
            // sample, each pair given one common error, the geometric mean
            // dd = sqrt(err(d_var) * err(d_nom)), dc likewise for the MC means:
            //   (dDR/DR)^2 = dd^2 (1/d_var - 1/d_nom)^2 + dc^2 (1/c_nom - 1/c_var)^2
            // No data-MC cross terms. Same formula as in Response_Ratios_Inclusive.h.
            {
                double d_var = data_means_after[m][bin], d_var_err = data_mean_errs_after[m][bin];
                double d_nom = data_means_before[m][bin], d_nom_err = data_mean_errs_before[m][bin];
                double c_var = comb_means_after[m][bin], c_var_err = comb_mean_errs_after[m][bin];
                double c_nom = comb_means_before[m][bin], c_nom_err = comb_mean_errs_before[m][bin];
                if (d_var != 0 && d_nom != 0 && c_var != 0 && c_nom != 0) {
                    double dd = std::sqrt(std::abs(d_var_err * d_nom_err));
                    double dc = std::sqrt(std::abs(c_var_err * c_nom_err));
                    double bd = (1.0 / d_var) - (1.0 / d_nom);
                    double bc = (1.0 / c_nom) - (1.0 / c_var);
                    double rel2 = bd * bd * dd * dd + bc * bc * dc * dc;
                    rjes_dr_errs_corr[m][bin] = std::abs(rjes_dr[m][bin]) * std::sqrt(rel2);

                    // Full: (dDR/DR)^2 = sum_ij c_i c_j rho_ij d_i d_j over
                    // (d_var, d_nom, c_var, c_nom), with c = (1/d_var, -1/d_nom,
                    // -1/c_var, 1/c_nom) and rho measured by Bootstrap_Correlations
                    // for this correction and cell (RespCorrLabel). Without a
                    // matching line the default matrix is used (rho = 1 for
                    // correction-nominal, 0 for data-MC).
                    CorrelationConfig::CorrelationMatrix Cresp =
                        CorrelationConfig::CorrFor(CorrelationConfig::CorrDomain::Response, CurrentChannel(),
                                                   fileWP, RespCorrLabel(m, bin), cfg_variant.path_tag);
                    double c_vec[4] = {1.0 / d_var, -1.0 / d_nom, -1.0 / c_var, 1.0 / c_nom};
                    double d_vec[4] = {d_var_err, d_nom_err, c_var_err, c_nom_err};
                    double rho[4][4] = {{1.0, Cresp.rho_mu_dd, Cresp.rho_mu_dv_mv, Cresp.rho_mu_dv_mn},
                                        {Cresp.rho_mu_dd, 1.0, Cresp.rho_mu_dn_mv, Cresp.rho_mu_dn_mn},
                                        {Cresp.rho_mu_dv_mv, Cresp.rho_mu_dn_mv, 1.0, Cresp.rho_mu_mm},
                                        {Cresp.rho_mu_dv_mn, Cresp.rho_mu_dn_mn, Cresp.rho_mu_mm, 1.0}};
                    double rel2f = 0.0;
                    for (int i = 0; i < 4; ++i)
                        for (int j = 0; j < 4; ++j)
                            rel2f += c_vec[i] * c_vec[j] * rho[i][j] * d_vec[i] * d_vec[j];
                    rel2f = std::max(0.0, rel2f);
                    rjes_dr_errs_full[m][bin] = std::abs(rjes_dr[m][bin]) * std::sqrt(rel2f);
                }
            }

            // DR(sigma_insitu), data side only
            double jer_before = CalculateResolution(data_sigmas_before[m][bin], truth_sigmas_before[m][bin]);
            double jer_before_err =
                CalculateResolutionError(data_sigmas_before[m][bin], data_sigma_errs_before[m][bin],
                                         truth_sigmas_before[m][bin], truth_sigma_errs_before[m][bin]);
            double jer_after = CalculateResolution(data_sigmas_after[m][bin], truth_sigmas_after[m][bin]);
            double jer_after_err =
                CalculateResolutionError(data_sigmas_after[m][bin], data_sigma_errs_after[m][bin],
                                         truth_sigmas_after[m][bin], truth_sigma_errs_after[m][bin]);

            sig_insitu_dr[m][bin] = SafeRatio(jer_after, jer_before);
            sig_insitu_dr_errs[m][bin] = SafeRatioErr(jer_after, jer_after_err, jer_before, jer_before_err);
            sig_insitu_dr_errs_corr[m][bin] = CorrelationAwareJERRatioErr(
                jer_after, data_sigmas_after[m][bin], truth_sigmas_after[m][bin],
                truth_sigma_errs_after[m][bin], jer_before, data_sigmas_before[m][bin],
                truth_sigmas_before[m][bin], truth_sigma_errs_before[m][bin]);

            // Full errors of the resolution double ratios: c^T Sigma c over the six
            // widths (data, combined MC, truth; correction and nominal) with the
            // bootstrap correlations of this cell (default matrix if none).
            const CorrelationConfig::CorrelationMatrix Cres =
                CorrelationConfig::CorrFor(CorrelationConfig::CorrDomain::Response, CurrentChannel(), fileWP,
                                           RespCorrLabel(m, bin), cfg_variant.path_tag);
            sig_insitu_dr_errs_full[m][bin] =
                std::abs(sig_insitu_dr[m][bin])
                * CorrelationConfig::ResolutionDRRelErrFull(
                    false, data_sigmas_after[m][bin], data_sigma_errs_after[m][bin], data_sigmas_before[m][bin],
                    data_sigma_errs_before[m][bin], comb_sigmas_after[m][bin], comb_sigma_errs_after[m][bin],
                    comb_sigmas_before[m][bin], comb_sigma_errs_before[m][bin], truth_sigmas_after[m][bin],
                    truth_sigma_errs_after[m][bin], truth_sigmas_before[m][bin], truth_sigma_errs_before[m][bin],
                    Cres);

            // DR(R_JER), with R_JER = sigma_insitu / sigma^MC_insitu from the
            // combined MC sigma. ComputeJERRatio propagates sigma_truth once,
            // since it enters both numerator and denominator of one R_JER.
            FitUtils::JERRatio rjer_before = FitUtils::ComputeJERRatio(
                data_sigmas_before[m][bin], data_sigma_errs_before[m][bin], comb_sigmas_before[m][bin],
                comb_sigma_errs_before[m][bin], truth_sigmas_before[m][bin], truth_sigma_errs_before[m][bin]);
            FitUtils::JERRatio rjer_after = FitUtils::ComputeJERRatio(
                data_sigmas_after[m][bin], data_sigma_errs_after[m][bin], comb_sigmas_after[m][bin],
                comb_sigma_errs_after[m][bin], truth_sigmas_after[m][bin], truth_sigma_errs_after[m][bin]);

            R_JER_var_arr[m][bin] = rjer_after.ratio;
            R_JER_var_err_arr[m][bin] = rjer_after.ratioErr;
            R_JER_nom_arr[m][bin] = rjer_before.ratio;
            R_JER_nom_err_arr[m][bin] = rjer_before.ratioErr;

            // b_JER data and MC values
            bJER_data_nom_arr[m][bin] = rjer_before.sig_insitu;
            bJER_data_nom_err_arr[m][bin] = rjer_before.sig_insitu_err;
            bJER_mc_nom_arr[m][bin] = rjer_before.sig_MC_insitu;
            bJER_mc_nom_err_arr[m][bin] = rjer_before.sig_MC_insitu_err;
            bJER_data_var_arr[m][bin] = rjer_after.sig_insitu;
            bJER_data_var_err_arr[m][bin] = rjer_after.sig_insitu_err;
            bJER_mc_var_arr[m][bin] = rjer_after.sig_MC_insitu;
            bJER_mc_var_err_arr[m][bin] = rjer_after.sig_MC_insitu_err;

            if (rjer_before.valid && rjer_after.valid) {
                rjer_dr[m][bin] = SafeRatio(rjer_after.ratio, rjer_before.ratio);
                rjer_dr_errs[m][bin] = SafeRatioErr(rjer_after.ratio, rjer_after.ratioErr, rjer_before.ratio,
                                                    rjer_before.ratioErr);
                rjer_dr_errs_corr[m][bin] = RJERDoubleRatioErrCorrAware(
                    rjer_after.ratio, data_sigmas_after[m][bin], comb_sigmas_after[m][bin],
                    truth_sigmas_after[m][bin], truth_sigma_errs_after[m][bin], rjer_before.ratio,
                    data_sigmas_before[m][bin], comb_sigmas_before[m][bin], truth_sigmas_before[m][bin],
                    truth_sigma_errs_before[m][bin]);
                rjer_dr_errs_full[m][bin] =
                    std::abs(rjer_dr[m][bin])
                    * CorrelationConfig::ResolutionDRRelErrFull(
                        true, data_sigmas_after[m][bin], data_sigma_errs_after[m][bin],
                        data_sigmas_before[m][bin], data_sigma_errs_before[m][bin], comb_sigmas_after[m][bin],
                        comb_sigma_errs_after[m][bin], comb_sigmas_before[m][bin], comb_sigma_errs_before[m][bin],
                        truth_sigmas_after[m][bin], truth_sigma_errs_after[m][bin], truth_sigmas_before[m][bin],
                        truth_sigma_errs_before[m][bin], Cres);
            } else {
                std::cerr << "  WARNING: R_JER invalid (negative variance) for " << massCutLabels[m] << " "
                          << CurrentPtBinLabel(bin) << " - DR(R_JER) will be 0 for this cell." << std::endl;
            }

            std::cout << "  [" << massCutLabels[m] << ", " << CurrentPtBinLabel(bin) << "]"
                      << "  DR(R_JES)=" << Form("%.4f±%.4f", rjes_dr[m][bin], rjes_dr_errs[m][bin])
                      << " (cor " << Form("%.4f", rjes_dr_errs_corr[m][bin]) << ", full "
                      << Form("%.4f", rjes_dr_errs_full[m][bin]) << ")"
                      << "  DR(sig_insitu)="
                      << Form("%.4f±%.4f", sig_insitu_dr[m][bin], sig_insitu_dr_errs[m][bin]) << " (cor "
                      << Form("%.4f", sig_insitu_dr_errs_corr[m][bin]) << ")"
                      << "  DR(R_JER)=" << Form("%.4f±%.4f", rjer_dr[m][bin], rjer_dr_errs[m][bin])
                      << " (cor " << Form("%.4f", rjer_dr_errs_corr[m][bin]) << ")" << std::endl;
        }
    }

    // 3. Heatmaps (mass bin x pT bin)
    std::string vl = cfg_variant.label;

    // Per-sample mean ratios, with their own colour range
    Draw2DHeatmap("ratio_data_mean_" + vl,
                  "Data R_{" + CurrentChannel() + "} Mean Ratio (" + vl + "/Nominal) - WP: " + labelWP,
                  "#LT R_{" + CurrentChannel() + "}^{Data, " + vl + "} #GT / #LT R_{" + CurrentChannel()
                      + "}^{Data, Nominal} #GT",
                  outputDir + "Data_" + cfg_nominal.branch_R + "_mean_ratio_" + vl + ".pdf", ratio_data_means,
                  ratio_data_mean_errs);

    Draw2DHeatmap("ratio_mc_mean_" + vl,
                  "MC R_{" + CurrentChannel() + "} Mean Ratio (" + vl + "/Nominal) - WP: " + labelWP,
                  "#LT R_{" + CurrentChannel() + "}^{MC, " + vl + "} #GT / #LT R_{" + CurrentChannel()
                      + "}^{MC, Nominal} #GT",
                  outputDir + "MC_" + cfg_nominal.branch_R + "_mean_ratio_" + vl + ".pdf", ratio_mc_means,
                  ratio_mc_mean_errs);

    Draw2DHeatmap("ratio_bkg_mean_" + vl,
                  "Background R_{" + CurrentChannel() + "} Mean Ratio (" + vl + "/Nominal) - WP: " + labelWP,
                  "#LT R_{" + CurrentChannel() + "}^{Bkg, " + vl + "} #GT / #LT R_{" + CurrentChannel()
                      + "}^{Bkg, Nominal} #GT",
                  outputDir + "Bkg_" + cfg_nominal.branch_R + "_mean_ratio_" + vl + ".pdf", ratio_bkg_means,
                  ratio_bkg_mean_errs);

    Draw2DHeatmap("ratio_comb_mean_" + vl,
                  "Combined MC R_{" + CurrentChannel() + "} Mean Ratio (" + vl + "/Nominal) - WP: " + labelWP,
                  "#LT R_{" + CurrentChannel() + "}^{CombMC, " + vl + "} #GT / #LT R_{" + CurrentChannel()
                      + "}^{CombMC, Nominal} #GT",
                  outputDir + "Combined_MC_" + cfg_nominal.branch_R + "_mean_ratio_" + vl + ".pdf",
                  ratio_comb_means, ratio_comb_mean_errs);

    // Fixed colour ranges for the double ratios, so the heatmaps of different
    // corrections and Run periods can be compared.
    double jesLo, jesHi, jerLo, jerHi;
    GetHeatmapRange("EASYBJETS_HEATMAP_JES_RANGE", 0.90, 1.10, jesLo, jesHi);
    GetHeatmapRange("EASYBJETS_HEATMAP_JER_RANGE", 0.60, 1.40, jerLo, jerHi);

    // DR(R_JES), uncorrelated errors
    Draw2DHeatmap("rjes_dr_" + vl, "b-JES Double Ratio (" + vl + "/Nominal) - WP: " + labelWP,
                  "b-JES^{" + vl + "} / b-JES^{Nominal}", outputDir + "R_JES_DoubleRatio_" + vl + ".pdf",
                  rjes_dr, rjes_dr_errs, jesLo, jesHi);

    // Per-sample sigma ratios
    Draw2DHeatmap("ratio_data_sigma_" + vl,
                  "Data R_{" + CurrentChannel() + "} Sigma Ratio (" + vl + "/Nominal) - WP: " + labelWP,
                  "#sigma^{Data, " + vl + "} / #sigma^{Data, Nominal}",
                  outputDir + "Data_" + cfg_nominal.branch_R + "_sigma_ratio_" + vl + ".pdf",
                  ratio_data_sigmas, ratio_data_sigma_errs);

    Draw2DHeatmap("ratio_mc_sigma_" + vl,
                  "MC R_{" + CurrentChannel() + "} Sigma Ratio (" + vl + "/Nominal) - WP: " + labelWP,
                  "#sigma^{MC, " + vl + "} / #sigma^{MC, Nominal}",
                  outputDir + "MC_" + cfg_nominal.branch_R + "_sigma_ratio_" + vl + ".pdf", ratio_mc_sigmas,
                  ratio_mc_sigma_errs);

    Draw2DHeatmap("ratio_bkg_sigma_" + vl,
                  "Background R_{" + CurrentChannel() + "} Sigma Ratio (" + vl + "/Nominal) - WP: " + labelWP,
                  "#sigma^{Bkg, " + vl + "} / #sigma^{Bkg, Nominal}",
                  outputDir + "Bkg_" + cfg_nominal.branch_R + "_sigma_ratio_" + vl + ".pdf", ratio_bkg_sigmas,
                  ratio_bkg_sigma_errs);

    Draw2DHeatmap("ratio_comb_sigma_" + vl,
                  "Combined MC R_{" + CurrentChannel() + "} Sigma Ratio (" + vl
                      + "/Nominal) - WP: " + labelWP,
                  "#sigma^{CombMC, " + vl + "} / #sigma^{CombMC, Nominal}",
                  outputDir + "Combined_MC_" + cfg_nominal.branch_R + "_sigma_ratio_" + vl + ".pdf",
                  ratio_comb_sigmas, ratio_comb_sigma_errs);

    // DR(sigma_insitu), uncorrelated and corr-aware errors. "COLZ TEXT" prints
    // only the values, so the two files currently look the same.
    Draw2DHeatmap("sig_insitu_dr_" + vl, "#sigma_{insitu} Double Ratio (" + vl + "/Nominal) - WP: " + labelWP,
                  "#sigma_{insitu}^{" + vl + "} / #sigma_{insitu}^{Nominal}",
                  outputDir + "sigma_insitu_DoubleRatio_" + vl + ".pdf", sig_insitu_dr, sig_insitu_dr_errs,
                  jerLo, jerHi);

    Draw2DHeatmap("sig_insitu_dr_corr_" + vl,
                  "#sigma_{insitu} Double Ratio corr-aware (" + vl + "/Nominal) - WP: " + labelWP,
                  "#sigma_{insitu}^{" + vl + "} / #sigma_{insitu}^{Nominal} ",
                  outputDir + "sigma_insitu_DoubleRatio_" + vl + "_CorrAware.pdf", sig_insitu_dr,
                  sig_insitu_dr_errs_corr, jerLo, jerHi);

    // DR(R_JER), uncorrelated and corr-aware errors
    Draw2DHeatmap("rjer_dr_" + vl, "b-JER Double Ratio (" + vl + "/Nominal) - WP: " + labelWP,
                  "b-JER^{" + vl + "} / b-JER^{Nominal}", outputDir + "R_JER_DoubleRatio_" + vl + ".pdf",
                  rjer_dr, rjer_dr_errs, jerLo, jerHi);

    Draw2DHeatmap("rjer_dr_corr_" + vl, "b-JER Double Ratio corr-aware (" + vl + "/Nominal) - WP: " + labelWP,
                  "b-JER^{" + vl + "} / b-JER^{Nominal} ",
                  outputDir + "R_JER_DoubleRatio_" + vl + "_CorrAware.pdf", rjer_dr, rjer_dr_errs_corr, jerLo,
                  jerHi);

    // 4. Double ratios vs pT

    // DR(R_JES), uncorrelated errors only; the corr-aware and full errors go to
    // the CSV and are drawn by the combined macros.
    MakeRatioSummaryPlot_vs_pT(outputDir, "rjes_summary_" + vl, "b-JES^{" + vl + "} / b-JES^{Nominal}",
                               "Summary_R_JES_vs_pT_" + fileWP + ".pdf", rjes_dr, rjes_dr_errs);

    // DR(sigma_insitu), uncorrelated then corr-aware
    MakeRatioSummaryPlot_vs_pT(
        outputDir, "sig_insitu_summary_" + vl, "#sigma_{insitu}^{" + vl + "} / #sigma_{insitu}^{Nominal}",
        "Summary_sigma_insitu_vs_pT_" + fileWP + ".pdf", sig_insitu_dr, sig_insitu_dr_errs);

    MakeRatioSummaryPlot_vs_pT(outputDir, "sig_insitu_summary_corr_" + vl,
                               "#sigma_{insitu}^{" + vl + "} / #sigma_{insitu}^{Nominal} ",
                               "Summary_sigma_insitu_vs_pT_" + fileWP + "_CorrAware.pdf", sig_insitu_dr,
                               sig_insitu_dr_errs_corr);

    // DR(R_JER), uncorrelated then corr-aware, fixed y range [0.45, 2.05]
    MakeRatioSummaryPlot_vs_pT(outputDir, "rjer_summary_" + vl, "b-JER^{" + vl + "} / b-JER^{Nominal}",
                               "Summary_R_JER_vs_pT_" + fileWP + ".pdf", rjer_dr, rjer_dr_errs, 0.45, 2.05);

    MakeRatioSummaryPlot_vs_pT(outputDir, "rjer_summary_corr_" + vl, "b-JER^{" + vl + "} / b-JER^{Nominal} ",
                               "Summary_R_JER_vs_pT_" + fileWP + "_CorrAware.pdf", rjer_dr, rjer_dr_errs_corr,
                               0.45, 2.05);

    // Uncorrelated vs corr-aware errors on the resolution double ratios
    MakeErrorComparisonPlot_vs_pT(outputDir, "sig_insitu_cmp_" + vl,
                                  "#sigma_{insitu}^{" + vl + "} / #sigma_{insitu}^{Nominal}",
                                  "sigma_insitu_DoubleRatio_ErrorComparison_vs_pT_" + fileWP + ".pdf",
                                  sig_insitu_dr, sig_insitu_dr_errs, sig_insitu_dr_errs_corr);

    MakeErrorComparisonPlot_vs_pT(outputDir, "rjer_cmp_" + vl, "b-JER^{" + vl + "} / b-JER^{Nominal}",
                                  "R_JER_DoubleRatio_ErrorComparison_vs_pT_" + fileWP + ".pdf", rjer_dr,
                                  rjer_dr_errs, rjer_dr_errs_corr, 0.45, 2.05);

    // 5. Text table and CSV
    std::string summaryPath = outputDir + "Summary_Ratios_" + fileWP + ".txt";
    GenerateSummaryTable(
        summaryPath, labelWP, variantLabel, ratio_data_means, ratio_data_mean_errs, ratio_mc_means,
        ratio_mc_mean_errs, ratio_bkg_means, ratio_bkg_mean_errs, ratio_comb_means, ratio_comb_mean_errs,
        ratio_data_sigmas, ratio_data_sigma_errs, ratio_mc_sigmas, ratio_mc_sigma_errs, ratio_bkg_sigmas,
        ratio_bkg_sigma_errs, ratio_comb_sigmas, ratio_comb_sigma_errs, rjes_dr, rjes_dr_errs, sig_insitu_dr,
        sig_insitu_dr_errs, sig_insitu_dr_errs_corr, rjer_dr, rjer_dr_errs, rjer_dr_errs_corr, R_JER_var_arr,
        R_JER_var_err_arr, R_JER_nom_arr, R_JER_nom_err_arr);

    // CSV for the combined macros; see GenerateSummaryCSV for the columns.
    std::string csvPath = outputDir + "Summary_Ratios_" + fileWP + ".csv";
    GenerateSummaryCSV(csvPath, CurrentChannel(), fileWP, variantLabel, rjes_dr, rjes_dr_errs, sig_insitu_dr,
                       sig_insitu_dr_errs, sig_insitu_dr_errs_corr, rjer_dr, rjer_dr_errs, rjer_dr_errs_corr,
                       R_JES_nom_arr, R_JES_nom_err_arr, R_JES_var_arr, R_JES_var_err_arr, R_JER_nom_arr,
                       R_JER_nom_err_arr, R_JER_var_arr, R_JER_var_err_arr,
                       // b_JES data/MC = mu_data, mu_combMC (nominal, variant)
                       data_means_before, data_mean_errs_before, comb_means_before, comb_mean_errs_before,
                       data_means_after, data_mean_errs_after, comb_means_after, comb_mean_errs_after,
                       // b_JER data/MC = sigma_insitu, sigma^MC_insitu (nominal, variant)
                       bJER_data_nom_arr, bJER_data_nom_err_arr, bJER_mc_nom_arr, bJER_mc_nom_err_arr,
                       bJER_data_var_arr, bJER_data_var_err_arr, bJER_mc_var_arr, bJER_mc_var_err_arr,
                       // columns 40-43: correlated double-ratio errors
                       rjes_dr_errs_corr, rjes_dr_errs_full, sig_insitu_dr_errs_full, rjer_dr_errs_full);

    std::cout << "\n✓ Ratio_Response complete for " << variantLabel << " (WP: " << fileWP << ")" << std::endl;
    std::cout << "  Output: " << outputDir << "\n" << std::endl;
}


// Entry point, called by the shims.
//   variant   "ALL" (or "") = every correction of the channel, or one key
//             such as "DB_bJER10v00"
//   channel   "DB" or "MJB"
//   massMode  "cuts" (3 mass bins) or "inclusive" (one MassInclusive bin)
inline void Ratio_Response_Run(const std::string& variant, const std::string& channel,
                               const std::string& massMode = "cuts") {
    // Set the global channel and mass mode before anything reads them.
    CurrentChannel() = channel;
    MassInclusiveMode() = (massMode == "inclusive");
    ConfigureMassBins();

    gROOT->SetBatch(kTRUE);
    TH1::AddDirectory(kFALSE);
    AtlasUtils::LoadAtlasStyle();

    auto workingPoints = AnalysisWPConfig::StandardWPs();

    std::vector<std::string> variants;
    if (variant == "ALL" || variant.empty()) {
        auto allKeys = GetPrimaryVariantKeys(CurrentChannel());
        for (const auto& k : allKeys) {
            if (k != CurrentChannel()) variants.push_back(k);
        }
    } else {
        variants.push_back(variant);
    }

    std::cout << "========================================================" << std::endl;
    std::cout << "        RATIO RESPONSE ANALYSIS (from Rdb fit files)    " << std::endl;
    std::cout << "========================================================" << std::endl;
    std::cout << "Variants: ";
    for (const auto& v : variants) std::cout << v << " ";
    std::cout << std::endl;
    std::cout << "Working Points: ";
    for (const auto& wp : workingPoints) std::cout << PathBuilder::WPtoFileWP(wp) << " ";
    std::cout << "\n========================================================\n" << std::endl;

    for (const auto& v : variants) {
        for (const auto& wp : workingPoints) {
            process_single_wp_ratio(v, wp);
        }
    }

    std::cout << "\n========================================================" << std::endl;
    std::cout << "  All Ratio_Response analyses complete." << std::endl;
    std::cout << "========================================================" << std::endl;
}

#endif // RATIO_RESPONSE_H
