/**
 * ============================================================================
 *  Bootstrap_Correlations.h - correlations between the correction and nominal
 *                             fit results, measured with a bootstrap
 * ============================================================================
 * ===========================================
 * @author:    Thiziri AMEZZA
 * @email:     thiziri.amezza@cern.ch
 * @date:      26 October 2026
 * ===========================================
 *
 * Runs a bootstrap of the analysis: in each replica every event gets a random
 * weight, and the fits of every observable and bin are redone for the bJR10
 * corrections (bJER10v00, bJER10v01) and the nominal. From the replicas come
 * two double-ratio errors used in the summaries:
 * 
 *   full    the correlations (rho) between the correction and nominal fit
 *              results, written to correlations.txt. The ratio macros
 *              propagate the fit errors with them (EASYBJETS_CORR_METHOD=file):
 *                (dDR/DR)^2 = sum_ij c_i c_j rho_ij (d_i/x_i)(d_j/x_j)
 *              x_i the fitted inputs (mu or sigma; data, MC, truth), d_i their
 *              fit errors, c_i = +/-1 (or the derivative, for sigma_insitu and
 *              R_JER) their weight in the double ratio. Correction and nominal
 *              enter with opposite signs, so their correlation reduces the error.
 * 
 *   bootstrap  the fitted values of every replica, saved in diagnostics/.
 *              Bootstrap_Analyze computes the double ratio in each replica,
 *                DR_r = (mu_d,v / mu_m,v)_r / (mu_d,n / mu_m,n)_r  (b-JMDP; same
 *                for the other observables),
 *              and its error is the spread of the DR_r over the replicas: half
 *              of the 16-84% interval (dr_errors.csv, read by the summaries).
 *              No rho and no fit errors are used: correction and nominal share
 *              the same events and weights in each replica, so their
 *              correlation is already in the spread.
 * 
 * Both use the bootstrap, for different things: "full" takes the size of each
 * error from the analysis fits and only the correlations from the bootstrap;
 * "bootstrap" takes everything from the replicas and uses no fit errors.
 * Comparing them checks the fit errors: if they agree, the fit errors are
 * complete; if "bootstrap" is larger, the fits miss a source of fluctuation
 * (e.g. in the mass fit, the subtracted background, which the peak fit error
 * treats as exactly known but which is refitted in every replica).
 *
 * One run covers every observable and every binning:
 *
 *   observable  double ratio of          rho used          ratio macro
 *   b-JMDP      mass peak mu             mass rho_mu_*     Mass_Ratios_PeakWidth
 *   b-JMDR      mass width sigma         mass rho_sig_*    Mass_Ratios_PeakWidth
 *   b-JES       response peak (R_JES)    resp rho_mu_*     Response_Ratios_*
 *   b-JER       sigma_insitu, R_JER      resp rho_sig_*    Response_Ratios_*
 *                                        (+ truth width)
 *
 *   binning                 line label in correlations.txt, e.g.
 *   mass, inclusive         mass <ch> <wp> <corr>:pt200to450           (PtBinConfig)
 *   mass, pT bins           mass <ch> <wp> <corr>:pt200to300
 *   resp, inclusive         resp <ch> <wp> <corr>:inclusive
 *   resp, pT bins           resp <ch> <wp> <corr>:200to300GeV          (BinConfig)
 *   resp, pT x mass bins    resp <ch> <wp> <corr>:200to300GeV_50-80GeV
 *
 * Each (observable, bin) is a "case" with its own selection. An event enters
 * every case its jet falls in, separately for the corrected and the nominal
 * jet, so all cases are filled in one loop and each gets its own rho. Each
 * ratio macro reads only the line of its correction and bin.
 *
 * Method
 *   The correction and the nominal must be measured on the same events, with
 *   the same random weight per event, otherwise their correlation is lost.
 *   Both are therefore read from the NOMINAL tree, which holds, for the
 *   selected jet, the nominal pT and mass (ljet_pt, ljet_m) and the bJR10 ones
 *   (ljet_bJR10v0X_pt, _mass). 
 *   The corrected response is: R_corr = R_nom * pT_corr / pT_nom
 *   since only the jet pT changes with the correction: the reference does not
 *   depend on the large-R jet calibration (DB: photon pT and the jet direction;
 *   MJB: the small-R recoil jets), as long as the same jet is selected.
 *   In each replica every event gets one Poisson(1) weight, shared by the
 *   correction and the nominal; each side then uses its own pT and mass for
 *   the pT bins and mass cuts, so events that move between bins are included.
 *   The histograms are filled and fitted as in the analysis:
 * 
 *     mass, data       the background function Mass_Fit chose for this variant
 *                      and bin (F-test + cross-validation) is refitted on the
 *                      sidebands of the replica and subtracted; the peak is
 *                      fitted with the iterative Gaussian of Mass_Fit (exp_2 if
 *                      Mass_Fit has not been run for that bin)
 *     mass, signal MC  iterative Gaussian of Mass_Fit
 *                      (mass histograms with the range and bins of Mass_Fit)
 *     R, R_truth       two-pass Gaussian around 1
 * 
 *   MC, as in the fits: for the mass, the Z(bb) signal alone (as Mass_Fit);
 *   for the response, the combined MC = truth-matched signal + main background
 *   weighted with its sideband SF(m) (inclusive and per pT bin, for each side),
 *   and the truth widths from truth-matched jets only
 *   (EASYBJETS_BOOT_MC=signal: signal only -> if you we want to check more).
 *   The selection follows the fits: pT bins in jet pT, mass bins 50-80, 80-110,
 *   >= 110 GeV, EASYBJETS_MASS_CUT for the response cases (not for the mass,
 *   as in Mass_Fit).
 *   The correlations are the Pearson coefficients over the replicas.
 *
 * Limitations
 *   - The minor backgrounds (ttgamma, Wqqgamma, Zqqgamma) are not in the
 *     combined MC; they are small.
 *   - The fits are quicker than in the analysis, to refit every case in every
 *     replica: the Mass_Fit background function is used but not re-chosen in
 *     each replica, and R is fitted with a two-pass Gaussian instead of the
 *     window scan of the response fits. The events and bins are the same.
 *   - Gaussian peaks only: EASYBJETS_FIT_MODEL (bukin, cb, dscb, ...) is not
 *     used, and the outputs go to Run{2,3}Plots/Bootstrap whatever the model.
 *   - MuonCorr cannot be measured (TODO, see below): it is a separate production 
 *     whose events cannot be matched to the nominal ones (no event or run number 
 *     in the slim ntuples). It is skipped and gets the default rho = 1 (vonservative
 *     choice since the muon-in-jet correction relative to nominal is almost zero). 
 * 
 *  TODO: reprocess the slimming so that the slim ntuples add the eventNumber and
 *  runNumber. This is the right way for every correction, not only MuonCorr:
 *  each variant could then be read from its own slim ntuple, with its own jet
 *  choice and selection, and its events matched to the nominal ones directly 
 *  by event and run number. A quicker fix for MuonCorr only: store the muon-corrected 
 *  jet in the nominal tree, next to the bJER jets.
 *
 * How to run (from ZbbGamma/macros or ZbbJets/macros)
 *
 *   First Mass_Fit.cxx("ALL","inclusive") and Mass_Fit.cxx("ALL","ptbins"),
 *   with the same EASYBJETS_RUN, _FIT_MODEL and _BINNING: the mass fits of the
 *   bootstrap use the background function Mass_Fit chose for each bin. The
 *   response cases need nothing beforehand: the bootstrap fits the background
 *   SF and R itself. Then:
 *
 *   export EASYBJETS_RUN=run3
 *   export EASYBJETS_DATA_DIR=$HOME/lpnatlas_data/ZbbGamma    (ZbbJets for MJB)
 *   export EASYBJETS_MASS_CUT=on               same setting as the response fits
 *   root 'Bootstrap_Correlations.cxx()'                 every bJER correction
 *   root 'Bootstrap_Correlations.cxx("DB_bJER10v01")'   one correction
 *   root 'Bootstrap_Analyze.cxx()'                      bootstrap DR errors
 *
 *   Then run the ratio macros (Mass_Ratios_PeakWidth, Response_Ratios_*, after
 *   their fits) and the summaries with EASYBJETS_CORR_METHOD=file and the same
 *   EASYBJETS_RUN; they find the files by themselves.
 *
 * Output (in the plot folder of the channel, one per run era)
 *   ../Run3Plots/Bootstrap/correlations.txt    (merged with the lines there)
 *   ../Run3Plots/Bootstrap/diagnostics/bootstrap_<corr>_<wp>_<mass|resp>_<bin>.txt
 *       fitted values of every replica, read by Bootstrap_Analyze
 *   (../Run2Plots for run2)
 *
 * Options (environment variables)
 *   EASYBJETS_BOOT_N        number of replicas (default 200)
 *   EASYBJETS_BOOT_SEED     random seed (default 42)
 *   EASYBJETS_BOOT_OUTPUT   output file (default above)
 *   EASYBJETS_BOOT_MC       MC of the response cases: combined (default) | signal
 *   EASYBJETS_RUN, EASYBJETS_DATA_DIR, EASYBJETS_MASS_CUT   as for the fits
 *
 * ==========================================================================
 * @version:  1.0
 * ==========================================================================
 */

#ifndef EASYBJETS_BOOTSTRAP_CORRELATIONS_H
#define EASYBJETS_BOOTSTRAP_CORRELATIONS_H

#include "VariantConfig.h"
#include "CommonUtils.h"
#include "BinningConfig.h"
#include "MassCutConfig.h"
#include "CorrelationConfig.h"
#include "SFUtils.h"          
#include "MassFit_Models.h"   

#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TRandom3.h>
#include <TSystem.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

// Histogram binning and fit windows of the bootstrap fits.
namespace BootstrapFit {
    // mass: range and bin count of Mass_Fit (variant mass_xmin, mass_xmax,
    // GetMaxBins of the pT bin), set where the histograms are made
    // response R and truth response: around 1, finer for the narrower truth
    constexpr int    kRBins = 50, kTruthBins = 100;
    constexpr double kRLo = 0.5, kRHi = 1.5;
    // two-pass Gaussian: wide window, peak guess, sigma limits
    constexpr double kRWideLo = 0.7, kRWideHi = 1.3, kRSigMin = 0.01, kRSigMax = 0.5;
    constexpr double kTWideLo = 0.8, kTWideHi = 1.2, kTSigMin = 0.003, kTSigMax = 0.3;
}  // namespace BootstrapFit

// MC used for the response cases: "combined" (default, signal + SF-weighted
// main background, as in the response fits) or "signal" (signal only).
// Set with EASYBJETS_BOOT_MC. The mass cases always use the signal only.
inline bool BootstrapUseCombinedMC() {
    const char* env = std::getenv("EASYBJETS_BOOT_MC");
    return !(env && std::string(env) == "signal");
}

// Current channel, set once by Bootstrap_Correlations_Run.
inline std::string& BootstrapCurrentChannel() {
    static std::string ch = "DB";
    return ch;
}

// Two-pass Gaussian fit: a wide fit, then a fit in mu +/- 2 sigma of the first.
// The ranges are arguments so the same code fits mass peaks (GeV) and the
// balance R (around 1). Returns zeros if a fit fails.
struct FitResult4 {
    double mu = 0, sigma = 0, mu_err = 0, sigma_err = 0;
};

// Defaults: the Z mass peak (GeV).
inline FitResult4 TwoPassGauss(TH1F* h, double wide_lo = 80.0, double wide_hi = 115.0,
                               double peak_guess = 90.0, double sigma_min = 0.1, double sigma_max = 40.0) {
    FitResult4 fr;
    if (!h || h->GetEntries() < 20) return fr;

    TF1 g1("gauss1", "gaus", wide_lo, wide_hi);
    g1.SetParameter(0, h->GetMaximum());
    g1.SetParameter(1, peak_guess);
    g1.SetParameter(2, h->GetRMS() * 0.5);
    g1.SetParLimits(2, sigma_min, sigma_max);
    if (h->Fit(&g1, "RQ0N") != 0) return fr;

    const double mu1 = g1.GetParameter(1);
    const double sigma1 = g1.GetParameter(2);
    if (!std::isfinite(mu1) || !std::isfinite(sigma1) || sigma1 <= 0) return fr;

    const double lo = std::max(h->GetXaxis()->GetXmin(), mu1 - 2.0 * sigma1);
    const double hi = std::min(h->GetXaxis()->GetXmax(), mu1 + 2.0 * sigma1);
    if (lo >= hi) return fr;

    TF1 g2("gauss2", "gaus", lo, hi);
    g2.SetParameter(0, g1.GetParameter(0));
    g2.SetParameter(1, mu1);
    g2.SetParameter(2, sigma1);
    g2.SetParLimits(2, sigma_min, (hi - lo) / 2);
    if (h->Fit(&g2, "RQ0N") != 0) return fr;

    const double mu2 = g2.GetParameter(1);
    const double sigma2 = g2.GetParameter(2);
    if (!std::isfinite(mu2) || !std::isfinite(sigma2) || sigma2 <= 0) return fr;

    fr.mu = mu2;
    fr.sigma = sigma2;
    fr.mu_err = g2.GetParError(1);
    fr.sigma_err = g2.GetParError(2);
    return fr;
}

// Mass peak fit as in Mass_Fit (MassFit_Signal.h, iterative mode): a first
// Gaussian fit in [80, 115] GeV, then refits in mu +/- 2 sigma until sigma
// changes by less than 2% (at most 5 iterations). Returns zeros if it fails.
inline FitResult4 IterativeMassGauss(TH1F* h) {
    FitResult4 fr = TwoPassGauss(h, 80.0, 115.0, 90.0, 1.0, 40.0);  // stage 1 + first refit
    if (!(fr.sigma > 0)) return fr;
    for (int it = 0; it < 4; ++it) {
        const double lo = std::max(h->GetXaxis()->GetXmin(), fr.mu - 2.0 * fr.sigma);
        const double hi = std::min(h->GetXaxis()->GetXmax(), fr.mu + 2.0 * fr.sigma);
        if (hi - lo < 3 * h->GetBinWidth(1)) break;  // window too narrow for this binning
        TF1 g("gauss_it", "gaus", lo, hi);
        g.SetParameters(h->GetBinContent(h->FindBin(fr.mu)), fr.mu, fr.sigma);
        g.SetParLimits(2, 1.0, 40.0);
        if (h->Fit(&g, "RQ0N") != 0) break;
        const double mu = g.GetParameter(1), sg = g.GetParameter(2);
        if (!std::isfinite(mu) || !std::isfinite(sg) || sg <= 0) break;
        const bool converged = std::abs(sg - fr.sigma) < 0.02 * fr.sigma;
        fr.mu = mu;
        fr.sigma = sg;
        fr.mu_err = g.GetParError(1);
        fr.sigma_err = g.GetParError(2);
        if (converged) break;
    }
    return fr;
}

// Background model chosen by Mass_Fit (F-test + cross-validation) for one
// variant, WP and pT bin, read from its optimal_fit_parameters.txt.
struct MassBkgModel {
    bool ok = false;
    std::string name;  // e.g. polyexp_3
    std::string type;  // PolyExp, Poly or Exp
    int order = 0;
    std::vector<double> pars;  // fitted parameters, used as starting values
};

// Reads <PlotDir>/<variant>/Mass_Fit_<ch>/<wp>/<ptbin>_<binning>/optimal_fit_parameters.txt
// (same folder as Mass_Fit with the same EASYBJETS_RUN, _FIT_MODEL, _BINNING).
inline MassBkgModel ReadMassFitOptimal(const VariantConfig& cfg, const std::string& channel,
                                       const std::string& fileWP, const std::string& ptLabel) {
    MassBkgModel m;
    const std::string path = RunConfig::PlotDir() + "/" + cfg.plot_subdir + "/Mass_Fit_" + channel + "/"
                             + fileWP + "/" + ptLabel + "_" + BinningConfig::ModeTag()
                             + "/optimal_fit_parameters.txt";
    std::ifstream in(path);
    if (!in.is_open()) return m;
    std::string line;
    bool inParam = false;
    while (std::getline(in, line)) {
        auto after = [&](const std::string& key) {
            const size_t p = line.find(key);
            if (p == std::string::npos) return std::string();
            std::string v = line.substr(p + key.size());
            v.erase(0, v.find_first_not_of(" \t"));
            v.erase(v.find_last_not_of(" \t\r") + 1);
            return v;
        };
        if (!after("Function Name:").empty()) m.name = after("Function Name:");
        else if (!after("Function Type:").empty()) m.type = after("Function Type:");
        else if (!after("Polynomial Order:").empty()) m.order = std::atoi(after("Polynomial Order:").c_str());
        else if (line.find("Parameter [") != std::string::npos) inParam = true;
        else if (inParam && !after("Value =").empty()) {
            m.pars.push_back(std::atof(after("Value =").c_str()));
            inParam = false;
        }
    }
    m.ok = !m.type.empty() && m.order > 0 && !m.pars.empty();
    return m;
}

// The Mass_Fit function of a model (MassFit_Models.h), starting from its fitted parameters. 
inline TF1* MakeMassBkg(const MassBkgModel& m, double xmin, double xmax, const std::string& name) {
    TF1* f = nullptr;
    if (m.type == "PolyExp") f = CreatePolyExpFit(m.order, xmin, xmax, name);
    else if (m.type == "Poly") f = CreatePolyFit(m.order, xmin, xmax, name);
    else if (m.type == "Exp") f = CreateExpFit(m.order, xmin, xmax, name);
    if (!f) return nullptr;
    for (int i = 0; i < f->GetNpar() && i < (int)m.pars.size(); ++i) f->SetParameter(i, m.pars[i]);
    return f;
}

// Data mass fit, as in Mass_Fit (Mass_Fit.h): the bins between the sidebands (signal region)
// are set to zero, the background function is fitted over the mass range, subtracted bin by bin 
// (clamped at 0), and the remaining Z peak is fitted like the Mass_Fit peak (IterativeMassGauss).
//   model ok  ->  the function Mass_Fit chose for this variant and bin
//   otherwise ->  exp_2 (exponential of a quadratic), as a fallback
// Returns zeros if a fit fails.
inline FitResult4 SidebandSubtractedGauss(TH1F* h, double sb1hi, double sb2lo, const MassBkgModel& model) {
    FitResult4 fr;
    if (!h || h->GetEntries() < 50) return fr;
    const double xmin = h->GetXaxis()->GetXmin(), xmax = h->GetXaxis()->GetXmax();

    // sideband histogram: signal window emptied, as in Mass_Fit
    std::unique_ptr<TH1F> hsb((TH1F*)h->Clone(Form("%s_sb", h->GetName())));
    hsb->SetDirectory(nullptr);
    const int binHigh1 = h->FindBin(sb1hi), binLow2 = h->FindBin(sb2lo);
    for (int b = binHigh1 + 1; b < binLow2; ++b) {
        hsb->SetBinContent(b, 0);
        hsb->SetBinError(b, 0);
    }

    std::unique_ptr<TF1> bkg;
    if (model.ok) {
        bkg.reset(MakeMassBkg(model, xmin, xmax, Form("%s_bkg", h->GetName())));
    } else {
        bkg.reset(CreateExpFit(2, xmin, xmax, Form("%s_bkg", h->GetName())));
        TF1 e1("bkg_expo", "expo", xmin, xmax);  // starting values for exp_2
        if (hsb->Fit(&e1, "RQ0N") == 0) bkg->SetParameters(e1.GetParameter(0), e1.GetParameter(1), 0.0);
    }
    if (!bkg || hsb->Fit(bkg.get(), "RQ0N") != 0) return fr;
    for (int i = 0; i < bkg->GetNpar(); ++i)
        if (!std::isfinite(bkg->GetParameter(i))) return fr;

    // data minus background, data errors kept
    std::unique_ptr<TH1F> sub((TH1F*)h->Clone(Form("%s_sub", h->GetName())));
    sub->SetDirectory(nullptr);
    for (int b = 1; b <= sub->GetNbinsX(); ++b)
        sub->SetBinContent(b, h->GetBinContent(b) - std::max(0.0, bkg->Eval(h->GetBinCenter(b))));

    fr = IterativeMassGauss(sub.get());
    return fr;
}

// Number of (side, mass bin) fits using the Mass_Fit model / the fallback,
// shown at the end and in the settings line of correlations.txt.
inline std::pair<int, int>& BootstrapMassModelCount() {
    static std::pair<int, int> n{0, 0};
    return n;
}

// Pearson correlation of two equally long series.
inline double PearsonCorrelation(const std::vector<double>& x, const std::vector<double>& y) {
    if (x.size() != y.size() || x.size() < 3) return 0.0;
    const double n = (double)x.size();
    double mx = 0, my = 0;
    for (size_t i = 0; i < x.size(); ++i) {
        mx += x[i];
        my += y[i];
    }
    mx /= n;
    my /= n;
    double num = 0, sx2 = 0, sy2 = 0;
    for (size_t i = 0; i < x.size(); ++i) {
        const double dx = x[i] - mx, dy = y[i] - my;
        num += dx * dy;
        sx2 += dx * dx;
        sy2 += dy * dy;
    }
    const double den = std::sqrt(sx2 * sy2);
    return (den > 1e-12) ? num / den : 0.0;
}

namespace BootstrapImpl {

// One measured case: its key in correlations.txt and the selection it uses.
struct CaseDef {
    bool resp = false;        // false = mass peak, true = response R
    std::string label;        // bin label written after "<corr>:"
    double ptLo = 0, ptHi = 1e30;
    int ptIdx = 0;             // mass: index of the bin in its PtBinConfig (bin count as Mass_Fit)
    double mLo = -1, mHi = -1; // mass bin [mLo, mHi]; -1 = no limit
    bool massCut = false;      // apply EASYBJETS_MASS_CUT
    int sfBin = -1;            // response: background SF of this pT bin (-1 = inclusive SF)
};

// One event of the nominal tree, nominal (n) and correction (c) side.
struct Event {
    float m_n = 0, m_c = 0;
    float pt_n = 0, pt_c = 0;
    float R_n = 0, R_c = 0, Rt = 0;
    float w = 1;
    bool truthOK = false;
    uint64_t mask_n = 0, mask_c = 0;  // cases this event enters, per side
};

// Does a jet with this pT and mass enter the case?
inline bool InCase(const CaseDef& k, double pt, double m) {
    if (pt < k.ptLo || pt >= k.ptHi) return false;
    if (k.mLo > 0 && m < k.mLo) return false;
    if (k.mHi > 0 && m > k.mHi) return false;
    if (k.massCut && MassCutConfig::RejectScalar(m)) return false;
    return true;
}

// All cases for one channel, with the bins of the macros that use them.
inline std::vector<CaseDef> BuildCases(const std::string& channel) {
    std::vector<CaseDef> cases;
    auto addMass = [&](const PtBinConfig::PtBinDef& d) {
        for (int i = 0; i < d.nBins; ++i) {
            const std::string lab = d.fileNames[i];
            bool dup = false;
            for (const auto& c : cases)
                if (!c.resp && c.label == lab) dup = true;
            if (dup) continue;
            CaseDef k;
            k.label = lab;
            k.ptLo = d.bins[i];
            k.ptHi = d.bins[i + 1];
            k.ptIdx = i;
            cases.push_back(k);
        }
    };
    addMass(PtBinConfig::InclusiveBins(channel));
    addMass(PtBinConfig::DefaultBins(channel));

    CaseDef incl;
    incl.resp = true;
    incl.label = "inclusive";
    incl.massCut = true;
    cases.push_back(incl);

    const bool isMJB = (channel == "MJB");
    const int nPt = isMJB ? BinConfig::nPtBins_3 : BinConfig::nPtBins_DB;
    const double* edges = isMJB ? BinConfig::ptBins_3 : BinConfig::ptBins_DB;
    const double mLo[3] = {50, 80, 110};
    const double mHi[3] = {80, 110, -1};
    const char* mDir[3] = {"50-80GeV", "80-110GeV", "MassGE110GeV"};
    for (int b = 0; b < nPt; ++b) {
        const std::string ptLab = isMJB ? BinConfig::PtBinFileLabel_3(b) : BinConfig::PtBinFileLabel_DB(b);
        CaseDef k;
        k.resp = true;
        k.massCut = true;
        k.ptLo = edges[b];
        k.ptHi = edges[b + 1];
        k.sfBin = b;
        k.label = CorrelationConfig::RespBinLabel(ptLab);
        cases.push_back(k);
        for (int m = 0; m < 3; ++m) {
            CaseDef km = k;
            km.mLo = mLo[m];
            km.mHi = mHi[m];
            km.label = CorrelationConfig::RespBinLabel(ptLab, mDir[m]);
            cases.push_back(km);
        }
    }
    return cases;
}

// Read one nominal tree into memory: the selected jet, nominal and correction.
// isBkg: main background MC; the minor backgrounds are not used. It is used 
// only in the response cases (the mass MC is the signal alone, as in Mass_Fit);
inline std::vector<Event> ReadEvents(const std::string& path, const VariantConfig& cfg_c,
                                     const VariantConfig& cfg_n, bool isMC,
                                     const std::vector<CaseDef>& cases, bool isBkg = false) {
    std::vector<Event> evs;
    std::unique_ptr<TFile> f(TFile::Open(path.c_str(), "READ"));
    if (!f || f->IsZombie()) {
        std::cerr << "[Bootstrap] WARN: cannot open " << path << "\n";
        return evs;
    }
    TTree* t = dynamic_cast<TTree*>(f->Get(cfg_n.tree_name.c_str()));
    if (!t) {
        std::cerr << "[Bootstrap] WARN: no tree '" << cfg_n.tree_name << "' in " << path << "\n";
        return evs;
    }
    for (const auto& br : {cfg_c.branch_mass, cfg_c.branch_pt}) {
        if (!Validate::HasBranch(t, br)) {
            std::cerr << "[Bootstrap] WARN: branch " << br << " missing in " << path << "\n";
            return evs;
        }
    }

    const std::string truthBr = "R_truth_" + cfg_n.channel;
    const std::string dRBr = Validate::HasBranch(t, "matched_tjet_dR") ? "matched_tjet_dR"
                                                                       : "matched_tjet_dR_sig";
    const bool hasTruth = isMC && Validate::HasBranch(t, truthBr) && Validate::HasBranch(t, dRBr);

    // Read only the branches used here.
    t->SetBranchStatus("*", 0);
    std::vector<float>*m_n = nullptr, *pt_n = nullptr, *m_c = nullptr, *pt_c = nullptr, *dR = nullptr;
    int idx = -1;
    double R = 0, Rt = -1, w = 1;
    auto on = [&](const std::string& b) { t->SetBranchStatus(b.c_str(), 1); };
    on("ljet_m");
    on("ljet_pt");
    on(cfg_c.branch_mass);
    on(cfg_c.branch_pt);
    on("selected_ljet_index");
    on(cfg_n.branch_R);
    t->SetBranchAddress("ljet_m", &m_n);
    t->SetBranchAddress("ljet_pt", &pt_n);
    t->SetBranchAddress(cfg_c.branch_mass.c_str(), &m_c);
    t->SetBranchAddress(cfg_c.branch_pt.c_str(), &pt_c);
    t->SetBranchAddress("selected_ljet_index", &idx);
    t->SetBranchAddress(cfg_n.branch_R.c_str(), &R);
    if (isMC) {
        on("weight");
        t->SetBranchAddress("weight", &w);
    }
    if (hasTruth) {
        on(truthBr);
        on(dRBr);
        t->SetBranchAddress(truthBr.c_str(), &Rt);
        t->SetBranchAddress(dRBr.c_str(), &dR);
    }

    const Long64_t n = t->GetEntries();
    evs.reserve(n);
    for (Long64_t i = 0; i < n; ++i) {
        t->GetEntry(i);
        if (idx < 0 || !m_n || idx >= (int)m_n->size() || idx >= (int)m_c->size()) continue;
        Event e;
        e.m_n = (*m_n)[idx];
        e.pt_n = (*pt_n)[idx];
        e.m_c = (*m_c)[idx];
        e.pt_c = (*pt_c)[idx];
        if (e.pt_n <= 0) continue;
        e.R_n = (float)R;
        e.R_c = (float)(R * e.pt_c / e.pt_n);  // same reference, corrected jet pT
        e.w = isMC ? (float)w : 1.0f;
        if (hasTruth && dR && idx < (int)dR->size()) {
            const double d = (*dR)[idx];
            e.truthOK = (Rt > 0 && d >= 0 && d < cfg_n.maxTruthDR);
            e.Rt = (float)Rt;
        }
        for (size_t k = 0; k < cases.size() && k < 64; ++k) {
            if (isBkg && !cases[k].resp) continue;  // mass MC = signal only
            // Response cases use only truth-matched signal MC, as the fits do.
            if (isMC && !isBkg && cases[k].resp && !e.truthOK) continue;
            if (InCase(cases[k], e.pt_n, e.m_n)) e.mask_n |= (uint64_t(1) << k);
            if (InCase(cases[k], e.pt_c, e.m_c)) e.mask_c |= (uint64_t(1) << k);
        }
        if (e.mask_n || e.mask_c) evs.push_back(e);
    }
    return evs;
}

// Fitted values of one case over the replicas. Index: 0 = correction, 1 = nominal.
struct CaseSeries {
    std::vector<double> mu_d[2], sig_d[2], mu_m[2], sig_m[2], sig_t[2];
    int nOk = 0;
};

}  // namespace BootstrapImpl

// Bootstrap of one correction at one WP: fills and fits every case in each
// replica and returns the correlation matrix per case label.
inline std::map<std::string, CorrelationConfig::CorrelationMatrix>
BootstrapOneCorrection(const std::string& channel, const std::string& wp_raw, const std::string& variant_key,
                       int nBootstrap, UInt_t seed, const std::string& diagnosticDir) {
    using namespace BootstrapImpl;
    std::map<std::string, CorrelationConfig::CorrelationMatrix> out;

    const VariantConfig& cfg_c = GetVariant(variant_key);
    const VariantConfig& cfg_n = GetVariant(channel);
    const std::string fileWP = PathBuilder::WPtoFileWP(wp_raw);
    const std::vector<CaseDef> cases = BuildCases(channel);
    if (cases.size() > 64) {
        std::cerr << "[Bootstrap] too many cases (" << cases.size() << "), only the first 64 are used\n";
    }
    const int nCases = (int)std::min<size_t>(cases.size(), 64);

    // Both sides are read from the nominal files.
    const std::string dataPath = PathBuilder::SamplePath(cfg_n, "data", wp_raw);
    const std::string mcPath =
        PathBuilder::SamplePath(cfg_n, SampleConfig::FileNames(cfg_n)[SampleConfig::kZbb], wp_raw);
    std::cout << "[Bootstrap] reading " << dataPath << "\n             and " << mcPath << "\n";
    const std::vector<Event> data = ReadEvents(dataPath, cfg_c, cfg_n, false, cases);
    const std::vector<Event> mc = ReadEvents(mcPath, cfg_c, cfg_n, true, cases);
    std::vector<Event> bkg;
    if (BootstrapUseCombinedMC()) {
        const std::string bkgPath =
            PathBuilder::SamplePath(cfg_n, SampleConfig::FileNames(cfg_n)[SampleConfig::kMultijets], wp_raw);
        std::cout << "             and " << bkgPath << "\n";
        bkg = ReadEvents(bkgPath, cfg_c, cfg_n, true, cases, true);
    }
    std::cout << "[Bootstrap] " << channel << " " << variant_key << " " << fileWP << ": " << data.size()
              << " data, " << mc.size() << " signal MC and " << bkg.size()
              << " background MC events in at least one case\n";
    if (data.size() < 50 || mc.size() < 50) {
        std::cerr << "[Bootstrap] too few events, nothing written for " << variant_key << "\n";
        return out;
    }

    // Background SF(m) per side, as in the response fits: data / background
    // fitted in the mass sidebands (after EASYBJETS_MASS_CUT), inclusive and
    // per pT bin, each side with its own pT, mass and sidebands. Fixed over
    // the replicas, like the SF in the fits.
    const bool isMJB = (channel == "MJB");
    const int nPtSF = isMJB ? BinConfig::nPtBins_3 : BinConfig::nPtBins_DB;
    const double* ptEdges = isMJB ? BinConfig::ptBins_3 : BinConfig::ptBins_DB;
    std::vector<SFUtils::LinearSF> sfIncl(2), sfPt[2];
    if (!bkg.empty()) {
        for (int s = 0; s < 2; ++s) {
            const VariantConfig& cs = (s == 0) ? cfg_c : cfg_n;
            std::vector<std::unique_ptr<TH1F>> hd, hb;
            for (int b = 0; b <= nPtSF; ++b) {  // b = nPtSF: inclusive
                hd.emplace_back(new TH1F(Form("sfd_%d_%d", s, b), "", 20, 50, 150));
                hb.emplace_back(new TH1F(Form("sfb_%d_%d", s, b), "", 20, 50, 150));
                hd.back()->Sumw2();
                hb.back()->Sumw2();
            }
            auto fillSF = [&](const std::vector<Event>& evs, std::vector<std::unique_ptr<TH1F>>& h, bool w) {
                for (const auto& e : evs) {
                    const double m = (s == 0) ? e.m_c : e.m_n;
                    const double pt = (s == 0) ? e.pt_c : e.pt_n;
                    if (MassCutConfig::RejectScalar(m)) continue;
                    const double wt = w ? e.w : 1.0;
                    h[nPtSF]->Fill(m, wt);
                    for (int b = 0; b < nPtSF; ++b)
                        if (pt >= ptEdges[b] && pt < ptEdges[b + 1]) h[b]->Fill(m, wt);
                }
            };
            fillSF(data, hd, false);
            fillSF(bkg, hb, true);
            auto fitSF = [&](int b) {
                return SFUtils::FitMassDependentSF(hd[b].get(), hb[b].get(), cs.sideband1_low,
                                                   cs.sideband1_high, cs.sideband2_low, cs.sideband2_max);
            };
            sfIncl[s] = fitSF(nPtSF);
            for (int b = 0; b < nPtSF; ++b) sfPt[s].push_back(fitSF(b));
        }
    }
    
    // Scale factor of one background event: SF(m) = p0 + p1*m, fitted on the
    // mass sidebands as in the response fits. It depends on the case k (the
    // inclusive SF, or the SF of the case's pT bin) and on the side s (0 =
    // correction, 1 = nominal: each side has its own SF, evaluated at its own
    // jet mass m). Multiplies the event's MC weight and its replica weight.
    // SF = 1 if its fit failed, the weight is 1 (no scaling).
    auto bkgWeight = [&](const CaseDef& k, int s, double m) {
        const SFUtils::LinearSF& sf =
            (k.sfBin >= 0 && k.sfBin < (int)sfPt[s].size()) ? sfPt[s][k.sfBin] : sfIncl[s];
        return SFUtils::GetEventWeight(1.0, m, sf);
    };

    // Mass_Fit background model of each mass case, per side (each side with its
    // own variant: the correction's Mass_Fit for the correction side).
    std::vector<MassBkgModel> bkgModel[2];
    for (int s = 0; s < 2; ++s) {
        const VariantConfig& cs = (s == 0) ? cfg_c : cfg_n;
        for (int k = 0; k < nCases; ++k) {
            MassBkgModel m;
            if (!cases[k].resp) {
                m = ReadMassFitOptimal(cs, channel, fileWP, cases[k].label);
                if (m.ok) {
                    ++BootstrapMassModelCount().first;
                } else {
                    ++BootstrapMassModelCount().second;
                    std::cerr << "[Bootstrap] WARN: no Mass_Fit result for " << cs.key << " " << fileWP << " "
                              << cases[k].label << " (run Mass_Fit first, same options); exp_2 used\n";
                }
            }
            bkgModel[s].push_back(m);
        }
    }

    // Histograms per case and side: data, MC and (response) MC truth.
    using namespace BootstrapFit;
    TH1::AddDirectory(false);
    std::vector<std::unique_ptr<TH1F>> hD[2], hM[2], hT[2];
    for (int s = 0; s < 2; ++s) {
        for (int k = 0; k < nCases; ++k) {
            const bool r = cases[k].resp;
            const std::string tag = Form("_%d_%d", s, k);
            if (r) {
                hD[s].emplace_back(new TH1F(("hD" + tag).c_str(), "", kRBins, kRLo, kRHi));
                hM[s].emplace_back(new TH1F(("hM" + tag).c_str(), "", kRBins, kRLo, kRHi));
                hT[s].emplace_back(new TH1F(("hT" + tag).c_str(), "", kTruthBins, kRLo, kRHi));
            } else {
                // range and bin count of Mass_Fit for this variant and pT bin
                const VariantConfig& cs = (s == 0) ? cfg_c : cfg_n;
                const int nb = std::max(5, cs.GetMaxBins(cases[k].ptIdx));
                hD[s].emplace_back(new TH1F(("hD" + tag).c_str(), "", nb, cs.mass_xmin, cs.mass_xmax));
                hM[s].emplace_back(new TH1F(("hM" + tag).c_str(), "", nb, cs.mass_xmin, cs.mass_xmax));
                hT[s].emplace_back(nullptr);
            }
            hD[s].back()->Sumw2();
            hM[s].back()->Sumw2();
            if (hT[s].back()) hT[s].back()->Sumw2();
        }
    }

    std::vector<CaseSeries> series(nCases);
    // One random stream per sample (data, signal, background). Data and MC are
    // different events, so they fluctuate independently; and the data replicas
    // stay the same whatever MC is read (EASYBJETS_BOOT_MC), so results with
    // signal-only and combined MC can be compared. Within a sample, the
    // correction and the nominal of an event share the same weight.
    TRandom3 rngData(seed), rngSig(seed + 1), rngBkg(seed + 2);

    for (int rep = 0; rep < nBootstrap; ++rep) {
        for (int s = 0; s < 2; ++s)
            for (int k = 0; k < nCases; ++k) {
                hD[s][k]->Reset();
                hM[s][k]->Reset();
                if (hT[s][k]) hT[s][k]->Reset();
            }

        // One Poisson weight per event, shared by the correction and the nominal.
        // MC = signal (+ SF-weighted background for the response cases); the
        // truth histograms take truth-matched jets only.
        auto fill = [&](const std::vector<Event>& evs, bool isMC, bool isBkg, TRandom3& rng) {
            for (const auto& e : evs) {
                const int kpois = rng.Poisson(1.0);
                if (kpois == 0) continue;
                const double w = e.w * kpois;
                for (int s = 0; s < 2; ++s) {
                    const uint64_t mask = (s == 0) ? e.mask_c : e.mask_n;
                    if (!mask) continue;
                    const double m = (s == 0) ? e.m_c : e.m_n;
                    const double R = (s == 0) ? e.R_c : e.R_n;
                    for (int k = 0; k < nCases; ++k) {
                        if (!(mask & (uint64_t(1) << k))) continue;
                        const double x = cases[k].resp ? R : m;
                        if (!isMC) {
                            hD[s][k]->Fill(x, w);
                            continue;
                        }
                        const double wk = isBkg ? w * bkgWeight(cases[k], s, m) : w;
                        hM[s][k]->Fill(x, wk);
                        if (hT[s][k] && e.truthOK) hT[s][k]->Fill(e.Rt, wk);
                    }
                }
            }
        };
        fill(data, false, false, rngData);
        fill(mc, true, false, rngSig);
        fill(bkg, true, true, rngBkg);

        for (int k = 0; k < nCases; ++k) {
            const bool r = cases[k].resp;
            FitResult4 fd[2], fm[2], ft[2];
            bool ok = true;
            for (int s = 0; s < 2 && ok; ++s) {
                TH1F* hd = hD[s][k].get();
                TH1F* hm = hM[s][k].get();
                if (r) {
                    fd[s] = TwoPassGauss(hd, kRWideLo, kRWideHi, 1.0, kRSigMin, kRSigMax);
                    fm[s] = TwoPassGauss(hm, kRWideLo, kRWideHi, 1.0, kRSigMin, kRSigMax);
                    ft[s] = TwoPassGauss(hT[s][k].get(), kTWideLo, kTWideHi, 1.0, kTSigMin, kTSigMax);
                } else {
                    // data: sideband background subtracted first, as in Mass_Fit
                    const VariantConfig& cs = (s == 0) ? cfg_c : cfg_n;
                    fd[s] = SidebandSubtractedGauss(hd, cs.sideband1_high, cs.sideband2_low, bkgModel[s][k]);
                    fm[s] = IterativeMassGauss(hm);  // signal MC: peak only, as in Mass_Fit
                }
                ok = fd[s].sigma > 0 && fm[s].sigma > 0 && (!r || ft[s].sigma > 0);
            }
            if (!ok) continue;
            CaseSeries& cs = series[k];
            for (int s = 0; s < 2; ++s) {
                cs.mu_d[s].push_back(fd[s].mu);
                cs.sig_d[s].push_back(fd[s].sigma);
                cs.mu_m[s].push_back(fm[s].mu);
                cs.sig_m[s].push_back(fm[s].sigma);
                if (r) cs.sig_t[s].push_back(ft[s].sigma);
            }
            ++cs.nOk;
        }

        if ((rep + 1) % 50 == 0 || rep + 1 == nBootstrap)
            std::cout << "  [Bootstrap] replica " << (rep + 1) << "/" << nBootstrap << "\n";
    }

    if (!diagnosticDir.empty()) gSystem->Exec(("mkdir -p " + diagnosticDir).c_str());

    for (int k = 0; k < nCases; ++k) {
        const CaseDef& cd = cases[k];
        const CaseSeries& cs = series[k];
        if (cs.nOk < 10) {
            std::cerr << "[Bootstrap] " << cd.label << ": only " << cs.nOk
                      << " replicas with all fits converged, no line written (default rho used)\n";
            continue;
        }
        // Index 0 = correction ("v" in the key names), 1 = nominal ("n").
        CorrelationConfig::CorrelationMatrix C;
        C.rho_mu_dd = PearsonCorrelation(cs.mu_d[0], cs.mu_d[1]);
        C.rho_mu_mm = PearsonCorrelation(cs.mu_m[0], cs.mu_m[1]);
        C.rho_sig_dd = PearsonCorrelation(cs.sig_d[0], cs.sig_d[1]);
        C.rho_sig_mm = PearsonCorrelation(cs.sig_m[0], cs.sig_m[1]);
        C.rho_mu_dv_mv = PearsonCorrelation(cs.mu_d[0], cs.mu_m[0]);
        C.rho_mu_dv_mn = PearsonCorrelation(cs.mu_d[0], cs.mu_m[1]);
        C.rho_mu_dn_mv = PearsonCorrelation(cs.mu_d[1], cs.mu_m[0]);
        C.rho_mu_dn_mn = PearsonCorrelation(cs.mu_d[1], cs.mu_m[1]);
        C.rho_sig_dv_mv = PearsonCorrelation(cs.sig_d[0], cs.sig_m[0]);
        C.rho_sig_dv_mn = PearsonCorrelation(cs.sig_d[0], cs.sig_m[1]);
        C.rho_sig_dn_mv = PearsonCorrelation(cs.sig_d[1], cs.sig_m[0]);
        C.rho_sig_dn_mn = PearsonCorrelation(cs.sig_d[1], cs.sig_m[1]);
        if (cd.resp) {
            C.rho_sig_tt = PearsonCorrelation(cs.sig_t[0], cs.sig_t[1]);
            C.rho_sig_mv_tv = PearsonCorrelation(cs.sig_m[0], cs.sig_t[0]);
            C.rho_sig_mv_tn = PearsonCorrelation(cs.sig_m[0], cs.sig_t[1]);
            C.rho_sig_mn_tv = PearsonCorrelation(cs.sig_m[1], cs.sig_t[0]);
            C.rho_sig_mn_tn = PearsonCorrelation(cs.sig_m[1], cs.sig_t[1]);
        }
        out[(cd.resp ? "resp " : "mass ") + cd.label] = C;

        std::cout << std::fixed << std::setprecision(3) << "  " << (cd.resp ? "resp " : "mass ") << std::left
                  << std::setw(26) << cd.label << std::right << " rho_mu_dd=" << C.rho_mu_dd
                  << " rho_mu_mm=" << C.rho_mu_mm << " rho_sig_dd=" << C.rho_sig_dd
                  << " rho_sig_mm=" << C.rho_sig_mm;
        if (cd.resp) std::cout << " rho_sig_tt=" << C.rho_sig_tt;
        std::cout << "  (" << cs.nOk << "/" << nBootstrap << " replicas)\n";

        // Diagnostic file: 8 columns per replica (the format Bootstrap_Analyze
        // reads), plus a TRUTH block for the response cases.
        if (!diagnosticDir.empty()) {
            const std::string dom = cd.resp ? "resp" : "mass";
            std::ofstream d(diagnosticDir + "bootstrap_" + variant_key + "_" + fileWP + "_" + dom + "_"
                            + cd.label + ".txt");
            d << (cd.resp ? "# RESP block: " : "# MASS block: ")
              << "mu_d_var mu_d_nom mu_m_var mu_m_nom sig_d_var sig_d_nom sig_m_var sig_m_nom\n";
            for (int i = 0; i < cs.nOk; ++i)
                d << std::setprecision(6) << cs.mu_d[0][i] << " " << cs.mu_d[1][i] << " " << cs.mu_m[0][i]
                  << " " << cs.mu_m[1][i] << " " << cs.sig_d[0][i] << " " << cs.sig_d[1][i] << " "
                  << cs.sig_m[0][i] << " " << cs.sig_m[1][i] << "\n";
            if (cd.resp) {
                d << "\n# TRUTH block: sig_t_var sig_t_nom\n";
                for (int i = 0; i < cs.nOk; ++i) d << cs.sig_t[0][i] << " " << cs.sig_t[1][i] << "\n";
            }
        }
    }
    return out;
}

// Entry point, called by Bootstrap_Correlations.cxx in each channel's macros/ folder (usage in the file header).
inline void Bootstrap_Correlations_Run(const std::string& variantOrAll, const std::string& channel) {
    BootstrapCurrentChannel() = channel;
    gROOT->SetBatch(kTRUE);
    MassCutConfig::PrintBanner("Bootstrap_Correlations");

    int nBootstrap = 200;
    UInt_t seed = 42;
    std::string outPath = CorrelationConfig::DefaultCorrelationsPath();
    if (const char* env = std::getenv("EASYBJETS_BOOT_N")) nBootstrap = std::atoi(env);
    if (const char* env = std::getenv("EASYBJETS_BOOT_SEED")) seed = (UInt_t)std::atoi(env);
    if (const char* env = std::getenv("EASYBJETS_BOOT_OUTPUT")) outPath = env;

    std::vector<std::string> variants;
    if (variantOrAll == "ALL" || variantOrAll == "all") {
        for (const auto& key : GetPrimaryVariantKeys(channel))
            if (key != channel) variants.push_back(key);
    } else {
        variants.push_back(variantOrAll);
    }

    std::cout << "\nBootstrap correlations: channel " << channel << ", " << nBootstrap
              << " replicas (seed " << seed << "), output " << outPath << "\n";

    // Existing lines, keyed by "<domain(observable)> <channel> <wp> <label>". Lines whose label has 
    // no "<corr>:" prefix come from the older bootstrap (independent weights, rho ~ 0) and are dropped.
    std::map<std::string, std::string> lines;
    int nDropped = 0;
    {
        std::ifstream in(outPath);
        std::string line;
        while (std::getline(in, line)) {
            if (line.empty() || line[0] == '#') continue;
            std::istringstream iss(line);
            std::string tok1, ch, wp, lab;
            if (!(iss >> tok1)) continue;
            std::string domain = "mass";
            if (tok1 == "mass" || tok1 == "resp") {
                domain = tok1;
                if (!(iss >> ch)) continue;
            } else {
                ch = tok1;
            }
            if (!(iss >> wp >> lab)) continue;
            if (lab.find(':') == std::string::npos) {
                ++nDropped;
                continue;
            }
            lines[domain + " " + ch + " " + wp + " " + lab] = line;
        }
    }
    if (nDropped > 0)
        std::cout << "[Bootstrap] dropped " << nDropped << " lines of " << outPath
                  << " without the \"<corr>:\" prefix in their label "
                  << " measured with independent weights)\n";

    const std::string outDir = outPath.substr(0, outPath.find_last_of('/') + 1);  // "" if no folder
    if (!outDir.empty()) gSystem->Exec(("mkdir -p " + outDir).c_str());
    const std::string diagDir = outDir + "diagnostics/";
    const auto workingPoints = AnalysisWPConfig::StandardWPs();

    auto formatLine = [](const std::string& domain, const std::string& ch, const std::string& wp,
                         const std::string& lab, const CorrelationConfig::CorrelationMatrix& C) {
        std::ostringstream o;
        o << std::fixed << std::setprecision(4) << domain << "  " << ch << "  " << wp << "  " << lab
          << "  rho_mu_dd=" << C.rho_mu_dd << "  rho_mu_mm=" << C.rho_mu_mm << "  rho_sig_dd=" << C.rho_sig_dd
          << "  rho_sig_mm=" << C.rho_sig_mm << "  rho_mu_dv_mv=" << C.rho_mu_dv_mv
          << "  rho_mu_dv_mn=" << C.rho_mu_dv_mn << "  rho_mu_dn_mv=" << C.rho_mu_dn_mv
          << "  rho_mu_dn_mn=" << C.rho_mu_dn_mn << "  rho_sig_dv_mv=" << C.rho_sig_dv_mv
          << "  rho_sig_dv_mn=" << C.rho_sig_dv_mn << "  rho_sig_dn_mv=" << C.rho_sig_dn_mv
          << "  rho_sig_dn_mn=" << C.rho_sig_dn_mn;
        if (domain == "resp")
            o << "  rho_sig_tt=" << C.rho_sig_tt << "  rho_sig_mv_tv=" << C.rho_sig_mv_tv
              << "  rho_sig_mv_tn=" << C.rho_sig_mv_tn << "  rho_sig_mn_tv=" << C.rho_sig_mn_tv
              << "  rho_sig_mn_tn=" << C.rho_sig_mn_tn;
        return o.str();
    };

    const VariantConfig& cfg_nom = GetVariant(channel);
    int iv = 0;
    for (const auto& variant : variants) {
        const VariantConfig& cfg_c = GetVariant(variant);

        // MuonCorr is a separate production with its own ntuples: its
        // events cannot be matched to the nominal ones, because the slim ntuples
        // store no event number or run number. It is skipped, and the ratio
        // macros use the default rho = 1 between MuonCorr and the nominal.
        // To measure it in the future, reprocess the slim ntuples so that they
        // keep eventNumber and runNumber (then the two productions can be matched
        // event by event), or so that the nominal tree also stores the
        // muon-corrected jet.
        if (cfg_c.path_tag == "MuonCorr" || cfg_c.branch_pt == cfg_nom.branch_pt) {
            std::cout << "[Bootstrap] " << variant << ": not in the nominal tree, cannot be matched"
                      << " event by event; skipped (the ratio macros use rho = 1 with the nominal)\n";
            continue;
        }
        for (const auto& wp : workingPoints) {
            const std::string fileWP = PathBuilder::WPtoFileWP(wp);
            const UInt_t caseSeed = seed + 1000 * (++iv);  // different random sequence per (correction, WP)
            auto res = BootstrapOneCorrection(channel, wp, variant, nBootstrap, caseSeed, diagDir);
            for (const auto& kv : res) {
                const std::string domain = kv.first.substr(0, 4);
                const std::string lab = cfg_c.path_tag + ":" + kv.first.substr(5);
                lines[domain + " " + channel + " " + fileWP + " " + lab] =
                    formatLine(domain, channel, fileWP, lab, kv.second);
            }
        }
    }

    std::ofstream out(outPath);
    out << "# Written by Bootstrap_Correlations\n"
        << "# Settings: " << (RunConfig::IsRun3() ? "run3" : "run2") << ", N=" << nBootstrap
        << ", MC=" << (BootstrapUseCombinedMC() ? "combined" : "signal")
        << ", seed=" << seed << ", EASYBJETS_MASS_CUT=" << (MassCutConfig::IsActive() ? "on" : "off")
        << ", mass background: Mass_Fit model " << BootstrapMassModelCount().first << ", exp_2 fallback "
        << BootstrapMassModelCount().second
        << " (lines kept from earlier runs may have other settings)\n"
        << "# Format: <domain> <channel> <wp> <correction>:<bin> <key=value>...\n"
        << "#   domain: mass (mass peak) or resp (response R)\n"
        << "#   correction and nominal of the same quantity: rho_mu_dd rho_mu_mm rho_sig_dd rho_sig_mm\n"
        << "#   (data d, MC m), rho_sig_tt (MC truth width, resp only)\n"
        << "#   data-MC pairs (v = correction, n = nominal): rho_mu_dv_mv ... rho_sig_dn_mn\n"
        << "#   MC reco-truth width pairs (resp only):\n"
        << "#     rho_sig_mv_tv rho_sig_mv_tn rho_sig_mn_tv rho_sig_mn_tn\n"
        << "# Read by CorrelationConfig.h.\n";
    for (const auto& kv : lines) out << kv.second << "\n";
    out.close();

    std::cout << "\nCorrelations written to " << outPath << "\n"
              << "Use them with EASYBJETS_CORR_METHOD=file (same EASYBJETS_RUN).\n";
}

#endif  // EASYBJETS_BOOTSTRAP_CORRELATIONS_H
