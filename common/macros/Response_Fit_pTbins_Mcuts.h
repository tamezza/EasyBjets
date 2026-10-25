/**
 * ============================================================================
 *  Response_Fit_pTbins_Mcuts.h - response fits in pT bins and jet-mass bins
 * ============================================================================
 * ===========================================
 * @author:    Thiziri AMEZZA
 * @email:     thiziri.amezza@cern.ch
 * @date:      08 April 2026
 * ===========================================
 *
 * Fits the balance R = pT(jet) / pT(reference) (branch R_DB or R_MJB) in pT bins
 * of the selected large-R jet, for data, the Z(bb) signal MC, the background
 * (gamma+jets for DB, multijets for MJB) and the combined MC, and computes the
 * b-JES and b-JER data/MC ratios in every cell.
 *
 * Two mass layouts share this code:
 * 
 *   cuts       Response_Fit_pTbins_Mcuts.cxx: three jet-mass bins 50-80, 80-110
 *              and >=110 GeV, one output folder each
 * 
 *   inclusive  Response_Fit_pTbins.cxx: one bin with no mass cut, folder
 *              MassInclusive
 * 
 * pT bins (Response_Fit_Common.h): DB 200-300, 300-450 GeV; 
 *                                  MJB 450-600, 600-900, 900-1200 GeV. 
 *                                  The last bin has no upper limit in the selection.
 *
 * For each variant and working point (e.g. 0p46, see AnalysisWPConfig):
 * 
 *   1. Open the data, signal, background and optional minor-background ntuples.
 * 
 *   2. Fill R per (mass bin, pT bin). The optional minor backgrounds have 
 *      no truth requirement in the produced samples. The combined-MC R
 *      (truth-matched signal + all background jets) therefore does not hold
 *      exactly the same events as the combined truth widths (truth-matched
 *      signal + truth-matched main background). If you consider the minor
 *      backgrounds in the analysis you may need to produce new ntuples 
 *      with the truth requirement.
 * 
 *   3. In each pT bin, fit a linear data/background scale factor in the jet
 *      mass sidebands and weight every background event with it.
 * 
 *   4. Fit data first (seed window taken from the histogram), then signal,
 *      background and combined MC (iterative Gaussian, or the model chosen
 *      with EASYBJETS_FIT_MODEL). The MC fits start from the data window, but
 *      each then iterates its own window (with EASYBJETS_SEED_MODE=raw the
 *      data window is not used at all).
 * 
 *   5. If the truth branches exist, fit the combined truth/ref and reco/truth
 *      distributions and compute, with the combined MC:
 *        sigma_insitu    = sqrt(sigma_data^2  - sigma_truth^2)
 *        sigma^MC_insitu = sqrt(sigma_combMC^2 - sigma_truth^2)
 *        Delta_closure   = (sigma^MC - sigma^MC_insitu) / sigma_insitu
 *        R_JES = <R>_data / <R>_combMC
 *        R_JER = sigma_insitu / sigma^MC_insitu
 *      If a truth fit fails badly, the whole cell is skipped (no R_JES, JER or
 *      plots for that pT bin and mass bin).
 * 
 *   6. Draw the fit, stack and ratio plots, the summaries vs pT and vs mass
 *      bin and the heat maps, and write text summaries and a binning
 *      recommendation.
 *
 * How to run (from ZbbGamma/macros for DB, from ZbbJets/macros for MJB):
 * 
 *   export EASYBJETS_RUN=run3 (run2 for 13 TeV)
 *   export EASYBJETS_DATA_DIR=$HOME/lpnatlas_data/ZbbGamma    (ZbbJets for MJB)
 *   export EASYBJETS_MASS_CUT=on       keep only jets with 50 < m < 150 GeV
 *                                      before the mass bins (>=110 becomes 110-150)
 *   root 'Response_Fit_pTbins_Mcuts.cxx()'                 nominal, 3 mass bins
 *   root 'Response_Fit_pTbins_Mcuts.cxx("DB_bJER10v00")'   one variant
 *   root 'Response_Fit_pTbins_Mcuts.cxx("ALL")'            all 4 variants
 *   root 'Response_Fit_pTbins.cxx("ALL")'                  all 4, mass inclusive
 *   (in ZbbJets/macros the variant keys are MJB, MJB_bJER10v00, ...)
 * 
 * Then run Response_Ratios_pTbins_Mcuts.cxx (after the mass-bin fit) or
 * Response_Ratios_pTbins.cxx (after the inclusive fit) with the same options.
 *
 * Output: <channel>/Run3Plots/<variant>/Response_Fit_pTbins_Mcuts_<DB|MJB>/0p46/
 * (Run2Plots for run2, Run{2,3}Plots_<model> for a non-Gaussian fit model)
 * 
 *   <mass bin>/   50-80GeV, 80-110GeV, MassGE110GeV, or MassInclusive:
 *                 {data,mc,bkg,combined_MC}_fit_parameters_ptbin<i>.txt (read by
 *                 the ratio macros), Fit_R_DB_*_ptbin<i>.pdf and _scan.pdf,
 *                 R_DB_stack_ratio_ptbin<i>.pdf, JER_results_ptbin<i>.txt,
 *                 R_DB_mean_vs_pT.pdf, R_DB_sigma_vs_pT.pdf,
 *                 R_DB_pT_binned_fit_summary.txt
 * 
 *   Summary/      sigma, R_JES and R_JER vs pT and heat maps,
 *                 Calibration_summary_0p46.pdf, JER_full_summary_0p46.txt
 * 
 *   WP folder     R_DB_*_vs_masscut_ptbin<i>.pdf, R_DB_*_ratio_summary.pdf,
 *                 scale_factor_fit_ptbin<i>.pdf, jer_closure_diagnostic.txt,
 *                 binning_recommendation.txt
 * 
 * File names say R_DB for MJB too. Both layouts write into the same WP folder,
 * so Summary/ and the WP-level files hold whichever layout ran last.
 *
 * Options (environment variables, read once per ROOT session):
 * 
 *   EASYBJETS_RUN                     run2 | run3 (default run3): ntuple folder and plot folder
 *   EASYBJETS_DATA_DIR                ntuple base folder (default ../../DataSets/ZbbGamma or ZbbJets)
 *   EASYBJETS_MASS_CUT                unset/off = no cut (default) | on = 50-150 GeV | lo:hi
 *   EASYBJETS_FIT_MODE                scan (default) | scan_dense | scan_log | legacy | legacy7grid
 *   EASYBJETS_FIT_NSIGMA              legacy window +/-N sigma (default 2.0, capped at 2.5 here)
 *   EASYBJETS_FIT_NSIGMA_LIST         scan windows, e.g. "1,1.5,2:1.5" (default 1,1.5; >2.5 dropped)
 *   EASYBJETS_FIT_MODEL               gauss (default) | all | bukin | cb | dscb | dgauss
 *   EASYBJETS_FIT_MODEL_NSIGMA_LO/HI  fit window of the non-Gaussian models (default 3.0 / 2.0 sigma)
 *   EASYBJETS_FIT_ERR_SCALE           pdg: errors on plots and in fit files x sqrt(chi2/ndf) if > 1
 *   EASYBJETS_SEED_MODE               truncated (default) | raw: seed mean +/- 2.5 RMS, data window ignored
 *   EASYBJETS_RESPONSE_FIT_LO         lowest R used in the fits (default 0.3; only matters above 0.5)
 *   EASYBJETS_RESPONSE_BINNING_MODE   static (default, bins from VariantConfig) | dynamic
 *   EASYBJETS_BINNING                 dynamic-mode algorithm: fixed (default) | fd | scott | fitfriendly | ...
 *   EASYBJETS_BINNING_TRUTH           dynamic mode: match (default) | auto | fixed | <algorithm>
 *   EASYBJETS_BINNING_MAX_BINS        dynamic mode: cap on the bin count (_TRUTH: truth only)
 *   EASYBJETS_BINNING_CORE_BINS       fit-friendly bins across the +/-2 sigma core (default 10)
 *   EASYBJETS_SIDEBAND1_LOW/HIGH      low mass sideband of the background SF fit (default VariantConfig)
 *   EASYBJETS_SIDEBAND2_LOW/MAX       high mass sideband of the background SF fit (default VariantConfig)
 *
 * ==========================================================================
 * @version:  1.0
 * ==========================================================================
 */
#ifndef RDB_PT_BINS_MCUT_H
#define RDB_PT_BINS_MCUT_H

#include "VariantConfig.h"
#include "CommonUtils.h"
#include "SFUtils.h"
#include "BinningConfig.h"
#include "MassCutConfig.h"
#include "Response_Fit_Common.h" // RDBFitResults + channel/pT-bin accessors (shared)
#include "Response_Fit_Core.h"   // shared response-fit core functions

#include <iomanip> // For std::setw
#include <fstream> // For file output
#include <limits>
#include "TPaletteAxis.h"

// Mass bins. Two layouts share this code:
//   cuts       3 bins 50-80, 80-110, >=110 GeV (default)
//   inclusive  1 bin with no mass cut, folder "MassInclusive" (Response_Fit_pTbins.cxx)
// The arrays are sized for MAX_MASS_CUTS and nMassCuts is set at run time by
// ConfigureMassBins(), so every loop over m < nMassCuts works in both layouts.
constexpr int MAX_MASS_CUTS = 6;
int nMassCuts = 3;
const char* massCutLabels[MAX_MASS_CUTS] = {"50-80GeV", "80-110GeV", "#geq110 GeV"}; // For display
const char* fileSystemMassCutLabels[MAX_MASS_CUTS] = {"50-80GeV", "80-110GeV",
                                                      "MassGE110GeV"}; // For file system
double massCutLow[MAX_MASS_CUTS] = {50, 80, 110};                      // Lower bound for each bin
double massCutHigh[MAX_MASS_CUTS] = {80, 110, -1}; // Upper bound for each bin (-1 means no upper bound)

// Mass layout for this session: false = 3 mass bins, true = one inclusive bin.
inline bool& MassInclusiveMode() {
    static bool m = false;
    return m;
}

// Fill the arrays above for the active layout. Safe to call more than once.
inline void ConfigureMassBins() {
    if (MassInclusiveMode()) {
        nMassCuts = 1;
        massCutLow[0] = -1;  // no lower cut
        massCutHigh[0] = -1; // no upper cut -> all masses in one bin
        massCutLabels[0] = "Inclusive mass";
        fileSystemMassCutLabels[0] = "MassInclusive";
    } else {
        nMassCuts = 3;
        massCutLow[0] = 50;
        massCutLow[1] = 80;
        massCutLow[2] = 110;
        massCutHigh[0] = 80;
        massCutHigh[1] = 110;
        massCutHigh[2] = -1;
        massCutLabels[0] = "50-80GeV";
        massCutLabels[1] = "80-110GeV";
        massCutLabels[2] = "#geq110 GeV";
        fileSystemMassCutLabels[0] = "50-80GeV";
        fileSystemMassCutLabels[1] = "80-110GeV";
        fileSystemMassCutLabels[2] = "MassGE110GeV";
    }
}


// Channel-dependent plot labels: DB gives R_{DB}, Z(bb)+gamma, gamma+jets;
// MJB gives R_{MJB}, Z(bb)+jets, multijets. They return string literals, so the
inline const char* RLab() { return (CurrentChannel() == "MJB") ? "R_{MJB}" : "R_{DB}"; }
inline const char* SigLab() {
    return (CurrentChannel() == "MJB") ? "Z(#rightarrow b#bar{b})+jets" : "Z(#rightarrow b#bar{b})+#gamma";
}
inline const char* BkgLab() { return (CurrentChannel() == "MJB") ? "Multijets" : "#gamma+jets"; }

// Create the WP output folder, its Summary/ folder and one folder per mass bin.
inline void CreateDirectories(const VariantConfig& cfg, const char* fileWP) {
    // Base directory for the specific working point
    TString baseDir =
        Form(((RunConfig::PlotDir() + "/%s/Response_Fit_pTbins_Mcuts_") + CurrentChannel() + "/%s").c_str(),
             cfg.plot_subdir.c_str(), fileWP);
    gSystem->Exec(Form("mkdir -p %s/Summary", baseDir.Data()));

    // Create a directory for each mass cut within the working point's folder
    for (int i = 0; i < nMassCuts; i++) {
        TString dirCmd = Form("mkdir -p %s/%s", baseDir.Data(), fileSystemMassCutLabels[i]);
        gSystem->Exec(dirCmd.Data());
    }
}

// Whole analysis for one variant and one working point. wp_raw is the WP as
// listed in AnalysisWPConfig, e.g. "46".
inline void ProcessWorkingPoint(const VariantConfig& cfg, const char* wp_raw) {
    Log::info(("Rdb" + CurrentChannel()).c_str(), Form("ENTER ProcessWorkingPoint: wp_raw=%s cfg.key=%s",
                                                       wp_raw ? wp_raw : "<NULL>", cfg.key.c_str()));
    // Convert raw WP to file-system format for output directories
    std::string fileWP_str = PathBuilder::WPtoFileWP(std::string(wp_raw));
    Log::info(("Rdb" + CurrentChannel()).c_str(),
              Form("After WPtoFileWP: fileWP_str=%s", fileWP_str.c_str()));
    const char* wp_label = fileWP_str.c_str(); // "0p46" - used in output paths

    std::cout << "\n--- Processing Working Point: " << wp_label << " ---" << std::endl;

    // Create directories for the current working point
    CreateDirectories(cfg, wp_label);

    // username and basePath are not used below: the ntuple paths come from
    // PathBuilder::SamplePath (EASYBJETS_DATA_DIR).
    std::string username = std::getenv("USER") ? std::getenv("USER") : "tamezza";
    std::string basePath = PathBuilder::GetEosPath();

    // Construct file paths - SamplePath internally calls WPtoFileWP, so pass raw wp
    std::string dataPath = PathBuilder::SamplePath(cfg, "data", std::string(wp_raw));
    std::string zbbPath =
        PathBuilder::SamplePath(cfg, SampleConfig::FileNames(cfg)[SampleConfig::kZbb], std::string(wp_raw));
    std::string mjPath = PathBuilder::SamplePath(cfg, SampleConfig::FileNames(cfg)[SampleConfig::kMultijets],
                                                 std::string(wp_raw));

    // Extra samples from SampleConfig (beyond data/signal/bkg)
    auto scFileNames = SampleConfig::FileNames(cfg);
    auto scDisplayNames = SampleConfig::DisplayNames(cfg);
    auto scColors = SampleConfig::Colors();
    const int nSC = SampleConfig::kNSamples;

    auto allPaths_rdb = SampleConfig::AllSamplePaths(cfg, std::string(wp_raw));

    // Open ROOT files
    SafeTFile fdata(dataPath.c_str());
    SafeTFile fmc(zbbPath.c_str());
    SafeTFile fmultijets(mjPath.c_str());

    if (!fdata || !fmc || !fmultijets) {
        std::cerr << "Error: Could not open one or more ROOT files for WP " << wp_label << ". Skipping."
                  << std::endl;
        return; 
    }

    // Minor backgrounds are optional: a missing file is only reported.
    TFile* fextra[SampleConfig::kNSamples] = {nullptr};
    bool extraLoaded[SampleConfig::kNSamples] = {false};
    for (int s = 0; s < nSC; ++s) {
        if (s == SampleConfig::kData || s == SampleConfig::kMultijets || s == SampleConfig::kZbb) continue;
        // "__unused_*" names only keep the sample list the same length for both
        // channels (SampleConfig in VariantConfig.h); there is no file to open.
        if (scFileNames[s].rfind("__unused_", 0) == 0) continue;
        fextra[s] = TFile::Open(allPaths_rdb[s].c_str(), "READ");
        if (!fextra[s] || fextra[s]->IsZombie()) {
            std::cerr << "  (" << scFileNames[s] << " not found - skipping)" << std::endl;
            if (fextra[s]) {
                delete fextra[s];
                fextra[s] = nullptr;
            }
        } else {
            extraLoaded[s] = true;
        }
    }

    TTree* tree_mc = Validate::GetTree(fmc.Get(), "nominal", "RdbAnalysis");
    TTree* tree_data = Validate::GetTree(fdata.Get(), "nominal", "RdbAnalysis");
    TTree* tree_mulijets = Validate::GetTree(fmultijets.Get(), "nominal", "RdbAnalysis");

    if (!tree_mc || !tree_data || !tree_mulijets) {
        std::cerr << "Error: Could not retrieve 'nominal' tree from one or more ROOT files for WP "
                  << wp_label << ". Skipping." << std::endl;

        for (int s = 0; s < nSC; ++s) {
            if (fextra[s]) {
                fextra[s]->Close();
                delete fextra[s];
            }
        }
        return;
    }

    // Get extra sample trees
    TTree* tree_extra[SampleConfig::kNSamples] = {nullptr};
    for (int s = 0; s < nSC; ++s) {
        if (extraLoaded[s]) {
            tree_extra[s] = (TTree*)fextra[s]->Get("nominal");
            if (!tree_extra[s]) extraLoaded[s] = false;
        }
    }

    // Create histograms for each mass cut and pT bin
    TH1F* h_mc_R_balance_ptbin[MAX_MASS_CUTS][MAX_PT_BINS];
    TH1F* h_data_R_balance_ptbin[MAX_MASS_CUTS][MAX_PT_BINS];
    TH1F* h_multijets_R_balance_ptbin[MAX_MASS_CUTS][MAX_PT_BINS];
    TH1F* h_combined_MC_R_balance_ptbin[MAX_MASS_CUTS][MAX_PT_BINS];
    TH1F* h_extra_R_balance_ptbin[SampleConfig::kNSamples][MAX_MASS_CUTS][MAX_PT_BINS];

    // Response histogram binning. Static (default): bin counts from VariantConfig.
    // Dynamic (EASYBJETS_RESPONSE_BINNING_MODE=dynamic): one binning computed from
    // a probe histogram of all signal MC and used for every (mass, pT) cell, so
    // the cells can still be compared with each other.
    BinningConfig::Result resp_bin_dyn;
    BinningConfig::Result truth_bin_dyn;
    BinningConfig::Result recotruth_bin_dyn; // dedicated, finer binning for reco/truth
    bool use_dynamic_resp = BinningConfig::ResponseDynamicEnabled();

    // Truth-binning policy (EASYBJETS_BINNING_TRUTH).
    BinningConfig::Mode truth_alg = BinningConfig::Mode::Fixed;
    BinningConfig::TruthPolicy truth_policy = BinningConfig::CurrentTruthPolicy(truth_alg);

    if (use_dynamic_resp) {
        std::cout << "[Rdb_pT] dynamic response binning enabled (reco_mode=" << BinningConfig::ModeTag()
                  << ", truth_policy=" << BinningConfig::TruthPolicyName(truth_policy);
        if (truth_policy == BinningConfig::TruthPolicy::Algorithm)
            std::cout << "(" << BinningConfig::ModeTag(truth_alg) << ")";
        std::cout << ")" << std::endl;

        const bool need_truth_probe = (truth_policy == BinningConfig::TruthPolicy::Auto)
                                      || (truth_policy == BinningConfig::TruthPolicy::Algorithm);

        TH1F* h_probe_R = new TH1F(Form("h_probe_Rdb_%s_%s", cfg.key.c_str(), wp_label),
                                   "R_response probe (reco/ref)", 100, 0.5, 1.5);
        h_probe_R->Sumw2();

        TH1F* h_probe_T = nullptr;
        if (need_truth_probe) {
            h_probe_T = new TH1F(Form("h_probe_T_Rdb_%s_%s", cfg.key.c_str(), wp_label),
                                 "R_response probe (reco/truth)", 100, 0.5, 1.5);
            h_probe_T->Sumw2();
        }

        // Single MC pass filling the probe histograms.
        bool have_truth_branch =
            FillResponseBinningHists(tree_mc, cfg, need_truth_probe, h_probe_R, h_probe_T, "[Rdb_pT]");
        // No truth branch: use the reco binning for truth too.
        if (need_truth_probe && !have_truth_branch) truth_policy = BinningConfig::TruthPolicy::Match;
        const int min_r = 10, max_r = 100;
        resp_bin_dyn = BinningConfig::ComputeResponse(h_probe_R, cfg.channel, 0.5, 1.5, min_r, max_r);

        switch (truth_policy) {
        case BinningConfig::TruthPolicy::Match:
            truth_bin_dyn = resp_bin_dyn;
            break;
        case BinningConfig::TruthPolicy::Fixed:
            // nbins = 0 tells MkRptTruth to use the per-cell bin count
            // from VariantConfig.
            truth_bin_dyn.tag = "fixed-truth";
            truth_bin_dyn.nbins = 0; // caller falls back to fallback_nb
            truth_bin_dyn.caption = "binning: fixed-truth (per-cell cfg.nBinsTruth)";
            break;
        case BinningConfig::TruthPolicy::Auto:
            truth_bin_dyn = BinningConfig::ComputeResponseWithMode(
                h_probe_T, cfg.channel, 0.5, 1.5, min_r, max_r, BinningConfig::CurrentMode(), "-truth");
            break;
        case BinningConfig::TruthPolicy::Algorithm:
            truth_bin_dyn = BinningConfig::ComputeResponseWithMode(h_probe_T, cfg.channel, 0.5, 1.5, min_r,
                                                                   max_r, truth_alg, "-truth");
            break;
        }

        // Own binning for reco/truth. Its core (sigma ~0.05) is much narrower
        // than truth/ref (sigma ~0.10), but long tails inflate its RMS, so the
        // shared truth binning would leave ~5 bins in the core (fit ndf 1-2).
        // Bin it from its core width instead; h_probe_T is the reco/truth probe.
        recotruth_bin_dyn = truth_bin_dyn; // default: same as truth/ref
        if (h_probe_T && h_probe_T->GetEntries() > 20) {
            int nrt = BinningConfig::ComputeFitFriendlyCore(h_probe_T, 0.5, 1.5, min_r, max_r);
            if (nrt > 0) {
                recotruth_bin_dyn.edges.clear();
                recotruth_bin_dyn.nbins = nrt;
                recotruth_bin_dyn.tag = "fitfriendly-recotruth";
                recotruth_bin_dyn.caption =
                    "binning: fit-friendly (reco/truth core), " + std::to_string(nrt) + " bins";
                std::cout << "[BinningConfig-Resp] channel=" << cfg.channel
                          << " mode=fitfriendly-recotruth nbins=" << nrt << " range=[0,2] (truth/ref keeps "
                          << truth_bin_dyn.nbins << ")\n";
            }
        }

        delete h_probe_R;
        if (h_probe_T) delete h_probe_T;
    }

    auto MkRpt = [&](const std::string& name, const std::string& title, int fallback_nb) -> TH1F* {
        if (use_dynamic_resp) {
            TH1F* h = BinningConfig::MakeHist(name, title, resp_bin_dyn, 0.5, 1.5);
            h->Sumw2();
            return h;
        }
        TH1F* h = new TH1F(name.c_str(), title.c_str(), fallback_nb, 0.5, 1.5);
        h->Sumw2();
        return h;
    };

    // Truth histograms (truth/ref):
    //   static mode                        fallback_nb (VariantConfig)
    //   dynamic, policy match (default)    same bins as reco
    //   dynamic, policy fixed              fallback_nb
    //   dynamic, policy auto / algorithm   truth_bin_dyn
    auto MkRptTruth = [&](const std::string& name, const std::string& title, int fallback_nb) -> TH1F* {
        if (use_dynamic_resp) {
            // Policy fixed sets nbins = 0: use fallback_nb.
            if (truth_bin_dyn.nbins > 0) {
                TH1F* h = BinningConfig::MakeHist(name, title, truth_bin_dyn, 0.5, 1.5);
                h->Sumw2();
                return h;
            }
        }
        TH1F* h = new TH1F(name.c_str(), title.c_str(), fallback_nb, 0.5, 1.5);
        h->Sumw2();
        return h;
    };

    // reco/truth histograms: like MkRptTruth, with the finer recotruth_bin_dyn.
    auto MkRptRecoTruth = [&](const std::string& name, const std::string& title, int fallback_nb) -> TH1F* {
        if (use_dynamic_resp && recotruth_bin_dyn.nbins > 0) {
            TH1F* h = BinningConfig::MakeHist(name, title, recotruth_bin_dyn, 0.5, 1.5);
            h->Sumw2();
            return h;
        }
        TH1F* h = new TH1F(name.c_str(), title.c_str(), fallback_nb, 0.5, 1.5);
        h->Sumw2();
        return h;
    };

    // Initialize histograms for each mass cut and pT bin
    for (int m = 0; m < nMassCuts; m++) {
        for (int i = 0; i < CurrentNPtBins(); i++) {
            // nBinsR_perPtBin[i], else nBinsR_perMassCut[m], else nBinsR.
            // Not used in dynamic mode.
            const int nb_r = cfg.GetNBinsR_pt(i, m);
            TString histName = Form(("h_mc_" + cfg.branch_R + "_masscut%s_ptbin%d_%s").c_str(),
                                    massCutLabels[m], i, wp_label);
            TString histTitle = Form((cfg.branch_R + " for %s (Mass Cut: %s)").c_str(), CurrentPtBinLabel(i),
                                     massCutLabels[m]);
            h_mc_R_balance_ptbin[m][i] = MkRpt(histName.Data(), histTitle.Data(), nb_r);
            h_mc_R_balance_ptbin[m][i]->SetFillColor(kRed);
            h_mc_R_balance_ptbin[m][i]->SetLineColor(kRed);

            histName = Form(("h_data_" + cfg.branch_R + "_masscut%s_ptbin%d_%s").c_str(), massCutLabels[m], i,
                            wp_label);
            h_data_R_balance_ptbin[m][i] = MkRpt(histName.Data(), histTitle.Data(), nb_r);
            h_data_R_balance_ptbin[m][i]->SetFillColor(kBlack);
            h_data_R_balance_ptbin[m][i]->SetLineColor(kBlack);

            histName = Form(("h_multijets_" + cfg.branch_R + "_masscut%s_ptbin%d_%s").c_str(),
                            massCutLabels[m], i, wp_label);
            h_multijets_R_balance_ptbin[m][i] = MkRpt(histName.Data(), histTitle.Data(), nb_r);
            h_multijets_R_balance_ptbin[m][i]->SetFillColor(kBlue);
            h_multijets_R_balance_ptbin[m][i]->SetLineColor(kBlue);

            histName = Form(("h_combined_MC_" + cfg.branch_R + "_masscut%s_ptbin%d_%s").c_str(),
                            massCutLabels[m], i, wp_label);
            histTitle = Form(("Combined MC " + cfg.branch_R + " for %s (Mass Cut: %s)").c_str(),
                             CurrentPtBinLabel(i), massCutLabels[m]);
            h_combined_MC_R_balance_ptbin[m][i] = MkRpt(histName.Data(), histTitle.Data(), nb_r);
            h_combined_MC_R_balance_ptbin[m][i]->SetFillColor(kGreen + 2);
            h_combined_MC_R_balance_ptbin[m][i]->SetLineColor(kGreen + 2);

            for (int s = 0; s < nSC; ++s) {
                if (s == SampleConfig::kData || s == SampleConfig::kMultijets || s == SampleConfig::kZbb)
                    continue;
                histName = Form("h_%s_R_DB_masscut%s_ptbin%d_%s", scFileNames[s].c_str(), massCutLabels[m], i,
                                wp_label);
                h_extra_R_balance_ptbin[s][m][i] = MkRpt(histName.Data(), histTitle.Data(), nb_r);
                h_extra_R_balance_ptbin[s][m][i]->SetFillColor(scColors[s]);
                h_extra_R_balance_ptbin[s][m][i]->SetLineColor(scColors[s]);
            }
        }
    }


    // Truth histograms for the JER closure, signal and background MC.
    TH1F* h_mc_truth_R_balance_ptbin[MAX_MASS_CUTS][MAX_PT_BINS]; // pT_truth / pT_ref  (sigma^truth)
    TH1F* h_mc_reco_truth_ptbin[MAX_MASS_CUTS]
                               [MAX_PT_BINS]; // R_response / R_truth = pT_reco / pT_truth  (sigma^MC)
    TH1F* h_bkg_truth_R_balance_ptbin[MAX_MASS_CUTS][MAX_PT_BINS];  // bkg pT_truth / pT_ref
    TH1F* h_bkg_reco_truth_ptbin[MAX_MASS_CUTS][MAX_PT_BINS];       // bkg pT_reco/pT_truth
    TH1F* h_comb_truth_R_balance_ptbin[MAX_MASS_CUTS][MAX_PT_BINS]; // combined truth (built after filling)
    TH1F* h_comb_reco_truth_ptbin[MAX_MASS_CUTS][MAX_PT_BINS];

    for (int m = 0; m < nMassCuts; m++) {
        for (int i = 0; i < CurrentNPtBins(); i++) {
            // nBinsTruth_perPtBin[i], else nBinsTruth_perMassCut[m], else
            // nBinsTruth: the truth peaks are narrower than the reco one, so they
            // get their own bin count. Dynamic mode may override it (see above).
            const int nb_truth = cfg.GetNBinsTruth_pt(i, m);
            h_mc_truth_R_balance_ptbin[m][i] = MkRptTruth(Form("h_mc_truth_R_DB_m%d_b%d_%s", m, i, wp_label),
                                                          "p_{T}^{truth}/p_{T}^{ref,db}", nb_truth);
            h_mc_reco_truth_ptbin[m][i] = MkRptRecoTruth(Form("h_mc_reco_truth_m%d_b%d_%s", m, i, wp_label),
                                                         "p_{T}^{j1}/p_{T}^{truth}", nb_truth);
            h_bkg_truth_R_balance_ptbin[m][i] =
                MkRptTruth(Form("h_bkg_truth_R_DB_m%d_b%d_%s", m, i, wp_label),
                           "p_{T}^{truth}/p_{T}^{ref,db} (bkg)", nb_truth);
            h_bkg_reco_truth_ptbin[m][i] = MkRptRecoTruth(Form("h_bkg_reco_truth_m%d_b%d_%s", m, i, wp_label),
                                                          "p_{T}^{j1}/p_{T}^{truth} (bkg)", nb_truth);
            h_comb_truth_R_balance_ptbin[m][i] = nullptr;
            h_comb_reco_truth_ptbin[m][i] = nullptr;
        }
    }

    // Declare variables for tree branches
    std::vector<float>* ljet_pt = nullptr;
    std::vector<float>* ljet_m = nullptr;
    double weight, weight_mulijets;
    double weight_extra;
    double R_response;
    int selected_jet_index = -1;

    // Set branch addresses
    tree_mc->SetBranchAddress(cfg.branch_pt.c_str(), &ljet_pt);
    tree_mulijets->SetBranchAddress(cfg.branch_pt.c_str(), &ljet_pt);

    // Jet mass
    tree_mc->SetBranchAddress(cfg.branch_mass.c_str(), &ljet_m);
    tree_mulijets->SetBranchAddress(cfg.branch_mass.c_str(), &ljet_m);

    // If the variant branch is missing in data, read the nominal ljet_pt / ljet_m.
    bool dataHasVariantPt = Validate::HasBranch(tree_data, cfg.branch_pt);
    bool dataHasVariantMass = Validate::HasBranch(tree_data, cfg.branch_mass);
    tree_data->SetBranchAddress(dataHasVariantPt ? cfg.branch_pt.c_str() : "ljet_pt", &ljet_pt);
    tree_data->SetBranchAddress(dataHasVariantMass ? cfg.branch_mass.c_str() : "ljet_m", &ljet_m);
    Log::info(("Rdb" + CurrentChannel()).c_str(),
              Form("Data pT branch:   %s", dataHasVariantPt ? cfg.branch_pt.c_str() : "ljet_pt [fallback]"));
    Log::info(
        ("Rdb" + CurrentChannel()).c_str(),
        Form("Data mass branch: %s", dataHasVariantMass ? cfg.branch_mass.c_str() : "ljet_m [fallback]"));

    tree_mc->SetBranchAddress("selected_ljet_index", &selected_jet_index);
    tree_data->SetBranchAddress("selected_ljet_index", &selected_jet_index);
    tree_mulijets->SetBranchAddress("selected_ljet_index", &selected_jet_index);

    tree_mc->SetBranchAddress(cfg.branch_R.c_str(), &R_response);
    tree_data->SetBranchAddress(cfg.branch_R.c_str(), &R_response);
    tree_mulijets->SetBranchAddress(cfg.branch_R.c_str(), &R_response);

    tree_mc->SetBranchAddress("weight", &weight);
    tree_mulijets->SetBranchAddress("weight", &weight_mulijets);

    // Set branch addresses for minor backgrounds
    for (int s = 0; s < nSC; ++s) {
        if (tree_extra[s]) {
            tree_extra[s]->SetBranchAddress(cfg.branch_pt.c_str(), &ljet_pt);
            tree_extra[s]->SetBranchAddress(cfg.branch_mass.c_str(), &ljet_m);
            tree_extra[s]->SetBranchAddress("selected_ljet_index", &selected_jet_index);
            tree_extra[s]->SetBranchAddress(cfg.branch_R.c_str(), &R_response);
            tree_extra[s]->SetBranchAddress("weight", &weight_extra);
        }
    }


    // JER truth branches (MC signal and background) + DeltaR quality branch
    double R_truth_DB = -1.0;
    double R_truth_DB_bkg = -1.0;
    bool hasTruthBranches = false;
    bool hasTruthBkg = false;
    std::vector<float>* matched_tjet_dR_sig = nullptr;
    std::vector<float>* matched_tjet_dR_bkg = nullptr;

    // The truth branch is R_truth_DB or R_truth_MJB; the variable R_truth_DB
    // holds the value for whichever channel is running.
    const std::string truthBranchName = "R_truth_" + CurrentChannel();

    if (Validate::HasBranch(tree_mc, truthBranchName.c_str())) {
        tree_mc->SetBranchAddress(truthBranchName.c_str(), &R_truth_DB);
        hasTruthBranches = true;
        Log::info(("Rdb" + CurrentChannel()).c_str(), "JER truth branches found (signal)");
    }
    if (Validate::HasBranch(tree_mc, "matched_tjet_dR")) {
        tree_mc->SetBranchAddress("matched_tjet_dR", &matched_tjet_dR_sig);
        Log::info(("Rdb" + CurrentChannel()).c_str(), "matched_tjet_dR branch found (signal)");
    }

    if (Validate::HasBranch(tree_mulijets, truthBranchName.c_str())) {
        tree_mulijets->SetBranchAddress(truthBranchName.c_str(), &R_truth_DB_bkg);
        hasTruthBkg = true;
        Log::info(("Rdb" + CurrentChannel()).c_str(), "JER truth branches found (bkg)");
    }
    if (Validate::HasBranch(tree_mulijets, "matched_tjet_dR")) {
        tree_mulijets->SetBranchAddress("matched_tjet_dR", &matched_tjet_dR_bkg);
        Log::info(("Rdb" + CurrentChannel()).c_str(), "matched_tjet_dR branch found (bkg)");
    }

    Log::info(("Rdb" + CurrentChannel()).c_str(), Form("Truth-match quality cut: dR < %.2f", cfg.maxTruthDR));

    // Signal MC. Only truth-matched jets are filled, the reco histogram
    // included, so the reco, truth/ref and reco/truth histograms hold the same
    // events (see the closure note below). Without truth branches the signal
    // histograms stay empty.
    Long64_t nentries = tree_mc->GetEntries();
    for (Long64_t i = 0; i < nentries; i++) {
        tree_mc->GetEntry(i);
        // Validate selected_jet_index first
        if (selected_jet_index < 0 || !ljet_pt || !ljet_m || selected_jet_index >= ljet_pt->size()
            || selected_jet_index >= ljet_m->size()) {
            continue; 
        }
        // Optional global mass window (EASYBJETS_MASS_CUT), applied before
        // the mass bins.
        if (MassCutConfig::Reject(ljet_m, selected_jet_index)) continue;

        double leadingJetPt = (*ljet_pt)[selected_jet_index]; 
        double jetMass = (*ljet_m)[selected_jet_index];      

        // Find the appropriate pT bin
        for (int bin = 0; bin < CurrentNPtBins(); bin++) {
            if (leadingJetPt >= CurrentPtBins()[bin] && leadingJetPt < CurrentPtBins()[bin + 1]) {
                // Apply different mass cuts
                for (int m = 0; m < nMassCuts; m++) {
                    bool passesCut = true;

                    // Apply lower mass cut if specified
                    if (massCutLow[m] > 0 && jetMass < massCutLow[m]) {
                        passesCut = false;
                    }

                    // Apply upper mass cut if specified
                    if (massCutHigh[m] > 0 && jetMass > massCutHigh[m]) {
                        passesCut = false;
                    }

                    // Fill histogram if passes cut
                    if (passesCut) {
                        // Compute dR(reco, truth) for the selected jet
                        double dR_truth = -1.0;
                        if (matched_tjet_dR_sig && selected_jet_index < (int)matched_tjet_dR_sig->size())
                            dR_truth = matched_tjet_dR_sig->at(selected_jet_index);

                        // Truth match: R_truth > 0 and dR < cfg.maxTruthDR. The closure
                        // sigma^2(R) = sigma^2(truth/ref) + sigma^2(reco/truth) only holds
                        // if the three histograms are filled with the same events.
                        bool truthOK =
                            hasTruthBranches && R_truth_DB > 0 && dR_truth >= 0 && dR_truth < cfg.maxTruthDR;

                        if (truthOK) {
                            h_mc_R_balance_ptbin[m][bin]->Fill(R_response, weight);
                            h_mc_truth_R_balance_ptbin[m][bin]->Fill(R_truth_DB, weight);
                            h_mc_reco_truth_ptbin[m][bin]->Fill(R_response / R_truth_DB, weight);
                        }
                    }
                }
                break;
            }
        }
    }

    // Data 
    nentries = tree_data->GetEntries();
    for (Long64_t i = 0; i < nentries; i++) {
        tree_data->GetEntry(i);

        // Validate selected_jet_index first
        if (selected_jet_index < 0 || !ljet_pt || !ljet_m || selected_jet_index >= ljet_pt->size()
            || selected_jet_index >= ljet_m->size()) {
            continue;
        }
        // Optional global mass window (EASYBJETS_MASS_CUT)
        if (MassCutConfig::Reject(ljet_m, selected_jet_index)) continue;

        double leadingJetPt = (*ljet_pt)[selected_jet_index]; 
        double jetMass = (*ljet_m)[selected_jet_index];      


        // Find the appropriate pT bin
        for (int bin = 0; bin < CurrentNPtBins(); bin++) {
            if (leadingJetPt >= CurrentPtBins()[bin] && leadingJetPt < CurrentPtBins()[bin + 1]) {
                // Apply different mass cuts
                for (int m = 0; m < nMassCuts; m++) {
                    bool passesCut = true;

                    // Apply lower mass cut if specified
                    if (massCutLow[m] > 0 && jetMass < massCutLow[m]) {
                        passesCut = false;
                    }

                    // Apply upper mass cut if specified
                    if (massCutHigh[m] > 0 && jetMass > massCutHigh[m]) {
                        passesCut = false;
                    }

                    // Fill histogram if passes cut
                    if (passesCut) {
                        h_data_R_balance_ptbin[m][bin]->Fill(R_response);
                    }
                }
                break;
            }
        }
    }

    // Step A: background scale factor per pT bin, SF(m) = p0 + p1 * m.
    // The jet mass of data and background is filled in each pT bin (no mass-bin
    // cut) and the data/background ratio is fitted in the sidebands of the
    // variant (e.g. DB nominal 50-75 and 110-160 GeV, MJB nominal 50-65 and
    // 110-150 GeV). The background is normalised to data where the signal is
    // small, and the SF is applied event by event in step B.
    double sf_p0[MAX_PT_BINS], sf_p1[MAX_PT_BINS];
    for (int bin = 0; bin < CurrentNPtBins(); bin++) {
        sf_p0[bin] = 1.0;
        sf_p1[bin] = 0.0;
    }

    {
        // Mass histograms per pT bin (inclusive in mass - no mass cut)
        TH1F* h_data_mass_ptbin[MAX_PT_BINS];
        TH1F* h_mj_mass_ptbin[MAX_PT_BINS];
        for (int bin = 0; bin < CurrentNPtBins(); bin++) {
            h_data_mass_ptbin[bin] =
                new TH1F(Form("h_data_mass_sf_ptbin%d_%s", bin, wp_label), "ljet_m", 20, 50, 150);
            h_mj_mass_ptbin[bin] =
                new TH1F(Form("h_mj_mass_sf_ptbin%d_%s", bin, wp_label), "ljet_m", 20, 50, 150);
            h_data_mass_ptbin[bin]->Sumw2();
            h_mj_mass_ptbin[bin]->Sumw2();
        }

        // Fill data mass per pT bin
        nentries = tree_data->GetEntries();
        for (Long64_t i = 0; i < nentries; i++) {
            tree_data->GetEntry(i);
            if (selected_jet_index < 0 || !ljet_pt || !ljet_m || selected_jet_index >= (int)ljet_pt->size()
                || selected_jet_index >= (int)ljet_m->size())
                continue;
            // With EASYBJETS_MASS_CUT=on, also set EASYBJETS_SIDEBAND2_MAX=150
            // so the high sideband stays inside the mass window.
            if (MassCutConfig::Reject(ljet_m, selected_jet_index)) continue;
            double pt = (*ljet_pt)[selected_jet_index];
            double m = (*ljet_m)[selected_jet_index];
            for (int bin = 0; bin < CurrentNPtBins(); bin++) {
                if (pt >= CurrentPtBins()[bin] && pt < CurrentPtBins()[bin + 1]) {
                    h_data_mass_ptbin[bin]->Fill(m);
                    break;
                }
            }
        }

        // Fill multijet mass per pT bin
        nentries = tree_mulijets->GetEntries();
        for (Long64_t i = 0; i < nentries; i++) {
            tree_mulijets->GetEntry(i);
            if (selected_jet_index < 0 || !ljet_pt || !ljet_m || selected_jet_index >= (int)ljet_pt->size()
                || selected_jet_index >= (int)ljet_m->size())
                continue;
            if (MassCutConfig::Reject(ljet_m, selected_jet_index)) continue;
            double pt = (*ljet_pt)[selected_jet_index];
            double m = (*ljet_m)[selected_jet_index];
            for (int bin = 0; bin < CurrentNPtBins(); bin++) {
                if (pt >= CurrentPtBins()[bin] && pt < CurrentPtBins()[bin + 1]) {
                    h_mj_mass_ptbin[bin]->Fill(m, weight_mulijets);
                    break;
                }
            }
        }

        // Linear fit of data/background in the sidebands (SFUtils)
        for (int bin = 0; bin < CurrentNPtBins(); bin++) {
            SFUtils::LinearSF sf = SFUtils::FitMassDependentSF(
                h_data_mass_ptbin[bin], h_mj_mass_ptbin[bin], cfg.sideband1_low, cfg.sideband1_high,
                cfg.sideband2_low, cfg.sideband2_max,
                Form(((RunConfig::PlotDir() + "/%s/Response_Fit_pTbins_Mcuts_") + CurrentChannel()
                      + "/%s/scale_factor_fit_ptbin%d.pdf")
                         .c_str(),
                     cfg.plot_subdir.c_str(), wp_label, bin),
                Form("%s WP %s %s", cfg.key.c_str(), wp_label, CurrentPtBinLabel(bin)));

            sf_p0[bin] = sf.p0;
            sf_p1[bin] = sf.p1;

            if (!sf.valid)
                Log::warning(("Rdb" + CurrentChannel()).c_str(),
                             Form("pT bin %d (%s): SF fit failed - using SF=1", bin, CurrentPtBinLabel(bin)));
        }

        for (int bin = 0; bin < CurrentNPtBins(); bin++) {
            delete h_data_mass_ptbin[bin];
            delete h_mj_mass_ptbin[bin];
        }
    }

    // Step B: background Response histograms, each event weighted by
    // weight x (p0[pT bin] + p1[pT bin] x jet mass).
    nentries = tree_mulijets->GetEntries();
    for (Long64_t i = 0; i < nentries; i++) {
        tree_mulijets->GetEntry(i);
        if (selected_jet_index < 0 || !ljet_pt || !ljet_m || selected_jet_index >= (int)ljet_pt->size()
            || selected_jet_index >= (int)ljet_m->size()) {
            continue;
        }
        // Optional global mass window (EASYBJETS_MASS_CUT)
        if (MassCutConfig::Reject(ljet_m, selected_jet_index)) continue;

        double leadingJetPt = (*ljet_pt)[selected_jet_index];
        double jetMass = (*ljet_m)[selected_jet_index];

        for (int bin = 0; bin < CurrentNPtBins(); bin++) {
            if (leadingJetPt >= CurrentPtBins()[bin] && leadingJetPt < CurrentPtBins()[bin + 1]) {
                // Per-event mass-dependent SF for this pT bin
                double sf = sf_p0[bin] + sf_p1[bin] * jetMass;
                double corrected_weight = weight_mulijets * sf;

                for (int m = 0; m < nMassCuts; m++) {
                    bool passesCut = true;
                    if (massCutLow[m] > 0 && jetMass < massCutLow[m]) passesCut = false;
                    if (massCutHigh[m] > 0 && jetMass > massCutHigh[m]) passesCut = false;
                    if (passesCut) {
                        h_multijets_R_balance_ptbin[m][bin]->Fill(R_response, corrected_weight);

                        // Truth-match quality cut: dR(reco, truth) < cfg.maxTruthDR
                        double dR_truth = -1.0;
                        if (matched_tjet_dR_bkg && selected_jet_index < (int)matched_tjet_dR_bkg->size())
                            dR_truth = matched_tjet_dR_bkg->at(selected_jet_index);

                        bool truthOK =
                            hasTruthBkg && R_truth_DB_bkg > 0 && dR_truth >= 0 && dR_truth < cfg.maxTruthDR;

                        if (truthOK) {
                            h_bkg_truth_R_balance_ptbin[m][bin]->Fill(R_truth_DB_bkg, corrected_weight);
                            h_bkg_reco_truth_ptbin[m][bin]->Fill(R_response / R_truth_DB_bkg,
                                                                 corrected_weight);
                        }
                    }
                }
                break;
            }
        }
    }

    // Step C: minor backgrounds, MC weight only (no SF).
    for (int s = 0; s < nSC; ++s) {
        if (!tree_extra[s]) continue;
        nentries = tree_extra[s]->GetEntries();
        for (Long64_t i = 0; i < nentries; i++) {
            tree_extra[s]->GetEntry(i);
            if (selected_jet_index < 0 || !ljet_pt || !ljet_m || selected_jet_index >= (int)ljet_pt->size()
                || selected_jet_index >= (int)ljet_m->size()) {
                continue;
            }
            // Optional global mass window (EASYBJETS_MASS_CUT)
            if (MassCutConfig::Reject(ljet_m, selected_jet_index)) continue;
            double leadingJetPt = (*ljet_pt)[selected_jet_index];
            double jetMass = (*ljet_m)[selected_jet_index];
            for (int bin = 0; bin < CurrentNPtBins(); bin++) {
                if (leadingJetPt >= CurrentPtBins()[bin] && leadingJetPt < CurrentPtBins()[bin + 1]) {
                    for (int m = 0; m < nMassCuts; m++) {
                        bool passesCut = true;
                        if (massCutLow[m] > 0 && jetMass < massCutLow[m]) passesCut = false;
                        if (massCutHigh[m] > 0 && jetMass > massCutHigh[m]) passesCut = false;
                        if (passesCut) {
                            h_extra_R_balance_ptbin[s][m][bin]->Fill(R_response, weight_extra);
                        }
                    }
                    break;
                }
            }
        }
    }

    // Step D: combined MC = signal + SF-weighted background + minor backgrounds.
    // No extra scaling: the SF is already in the background event weights.
    for (int m = 0; m < nMassCuts; m++) {
        for (int bin = 0; bin < CurrentNPtBins(); bin++) {
            h_combined_MC_R_balance_ptbin[m][bin]->Add(h_mc_R_balance_ptbin[m][bin]);
            h_combined_MC_R_balance_ptbin[m][bin]->Add(h_multijets_R_balance_ptbin[m][bin]);
            for (int s = 0; s < nSC; ++s)
                if (extraLoaded[s])
                    h_combined_MC_R_balance_ptbin[m][bin]->Add(h_extra_R_balance_ptbin[s][m][bin]);

            // Build combined truth = signal + bkg (for JER)
            h_comb_truth_R_balance_ptbin[m][bin] = (TH1F*)h_mc_truth_R_balance_ptbin[m][bin]->Clone(
                Form("h_comb_truth_R_DB_m%d_b%d_%s", m, bin, wp_label));
            if (hasTruthBkg) h_comb_truth_R_balance_ptbin[m][bin]->Add(h_bkg_truth_R_balance_ptbin[m][bin]);

            h_comb_reco_truth_ptbin[m][bin] = (TH1F*)h_mc_reco_truth_ptbin[m][bin]->Clone(
                Form("h_comb_reco_truth_m%d_b%d_%s", m, bin, wp_label));
            if (hasTruthBkg) h_comb_reco_truth_ptbin[m][bin]->Add(h_bkg_reco_truth_ptbin[m][bin]);
        }
    }

    // Create arrays to store fit results for each mass cut
    double mc_means[MAX_MASS_CUTS][MAX_PT_BINS], mc_mean_errs[MAX_MASS_CUTS][MAX_PT_BINS];
    double data_means[MAX_MASS_CUTS][MAX_PT_BINS], data_mean_errs[MAX_MASS_CUTS][MAX_PT_BINS];
    double bkg_means[MAX_MASS_CUTS][MAX_PT_BINS], bkg_mean_errs[MAX_MASS_CUTS][MAX_PT_BINS];
    double combined_means[MAX_MASS_CUTS][MAX_PT_BINS], combined_mean_errs[MAX_MASS_CUTS][MAX_PT_BINS];

    double mc_sigmas[MAX_MASS_CUTS][MAX_PT_BINS], mc_sigma_errs[MAX_MASS_CUTS][MAX_PT_BINS];
    double data_sigmas[MAX_MASS_CUTS][MAX_PT_BINS], data_sigma_errs[MAX_MASS_CUTS][MAX_PT_BINS];
    double bkg_sigmas[MAX_MASS_CUTS][MAX_PT_BINS], bkg_sigma_errs[MAX_MASS_CUTS][MAX_PT_BINS];
    double combined_sigmas[MAX_MASS_CUTS][MAX_PT_BINS], combined_sigma_errs[MAX_MASS_CUTS][MAX_PT_BINS];

    // JER quantities per (mass cut, pT bin) - filled during truth fits
    double jer_sigma_RDB_data[MAX_MASS_CUTS][MAX_PT_BINS] = {},
           jer_sigma_RDB_data_err[MAX_MASS_CUTS][MAX_PT_BINS] = {};
    double jer_sigma_RDB_MC[MAX_MASS_CUTS][MAX_PT_BINS] = {},
           jer_sigma_RDB_MC_err[MAX_MASS_CUTS][MAX_PT_BINS] = {};
    double jer_sigma_truth[MAX_MASS_CUTS][MAX_PT_BINS] = {},
           jer_sigma_truth_err[MAX_MASS_CUTS][MAX_PT_BINS] = {};
    double jer_sigma_reco_truth[MAX_MASS_CUTS][MAX_PT_BINS] = {},
           jer_sigma_reco_truth_err[MAX_MASS_CUTS][MAX_PT_BINS] = {};
    double jer_sigma_insitu[MAX_MASS_CUTS][MAX_PT_BINS] = {},
           jer_sigma_insitu_err[MAX_MASS_CUTS][MAX_PT_BINS] = {};
    double jer_sigma_MC_insitu[MAX_MASS_CUTS][MAX_PT_BINS] = {},
           jer_sigma_MC_insitu_err[MAX_MASS_CUTS][MAX_PT_BINS] = {};
    double jer_ratio[MAX_MASS_CUTS][MAX_PT_BINS] = {}, jer_ratio_err[MAX_MASS_CUTS][MAX_PT_BINS] = {};
    double jes_ratio[MAX_MASS_CUTS][MAX_PT_BINS] = {}, jes_ratio_err[MAX_MASS_CUTS][MAX_PT_BINS] = {};
    double jer_sigma_MC[MAX_MASS_CUTS][MAX_PT_BINS] = {}, jer_sigma_MC_err[MAX_MASS_CUTS][MAX_PT_BINS] = {};
    double jer_delta_closure[MAX_MASS_CUTS][MAX_PT_BINS] = {},
           jer_delta_closure_err[MAX_MASS_CUTS][MAX_PT_BINS] = {};
    bool jer_valid[MAX_MASS_CUTS][MAX_PT_BINS] = {};

    double ptBinCenters[MAX_PT_BINS];
    for (int i = 0; i < CurrentNPtBins(); i++) {
        ptBinCenters[i] = (CurrentPtBinsDisplay()[i] + CurrentPtBinsDisplay()[i + 1]) / 2;
        // The last bin uses its display edge (450 or 1200 GeV), not the open selection edge.
        if (i == CurrentNPtBins() - 1) {
            ptBinCenters[i] =
                CurrentPtBinsDisplay()[i] + (CurrentPtBinsDisplay()[i + 1] - CurrentPtBinsDisplay()[i]) / 2;
        }
    }

    // For each mass cut, perform fits and create plots
    for (int m = 0; m < nMassCuts; m++) {
        // Create a summary file for this mass cut
        TString summaryFileName = Form(((RunConfig::PlotDir() + "/%s/Response_Fit_pTbins_Mcuts_")
                                        + CurrentChannel() + "/%s/%s/R_DB_pT_binned_fit_summary.txt")
                                           .c_str(),
                                       cfg.plot_subdir.c_str(), wp_label, fileSystemMassCutLabels[m]);
        ofstream summaryFile(summaryFileName.Data());

        if (summaryFile.is_open()) {
            summaryFile << "Summary of R_response Gaussian Fits in pT Bins (Mass Cut: " << massCutLabels[m]
                        << ", WP: " << wp_label << ")" << std::endl;
            summaryFile << "----------------------------------------" << std::endl << std::endl;
        }

        // Perform fits for each pT bin with this mass cut
        for (int bin = 0; bin < CurrentNPtBins(); bin++) {
            std::string plotBase =
                Form(((RunConfig::PlotDir() + "/%s/Response_Fit_pTbins_Mcuts_") + CurrentChannel() + "/%s/%s")
                         .c_str(),
                     cfg.plot_subdir.c_str(), wp_label, fileSystemMassCutLabels[m]);

            // DoFit: iterative fit of one Response histogram.
            // Without seedLo/seedHi the seed window comes from the histogram:
            // mean +/- 2 RMS, with mean and RMS recomputed inside +/-3 RMS to limit
            // the tails (as DataDrivenSeed does). It adapts to each bin, which
            // matters when a bin has few events.
            // Data is fitted first and its final window is passed as the seed
            // window of the signal, background and combined MC fits, so all
            // widths are measured over similar ranges. Each MC fit still
            // iterates its own window from that seed, and EASYBJETS_SEED_MODE=raw
            // ignores the seed altogether.
            auto DoFit = [&](TH1F* h, const char* tag_prefix, const char* label,
                             const std::string& scanPdfPath, double seedLo = -1,
                             double seedHi = -1) -> FitUtils::GausFitResult {
                // Data-driven seed window (only if caller didn't pass one)
                if (seedLo < 0 || seedHi < 0) {
                    double m1 = h->GetMean();
                    double s1 = h->GetRMS();
                    if (s1 > 0) {
                        int saveFirst = h->GetXaxis()->GetFirst();
                        int saveLast = h->GetXaxis()->GetLast();
                        h->GetXaxis()->SetRangeUser(m1 - 3 * s1, m1 + 3 * s1);
                        double m2 = h->GetMean();
                        double s2 = h->GetRMS();
                        h->GetXaxis()->SetRange(saveFirst, saveLast); // restore the axis range
                        if (s2 > 0) {
                            m1 = m2;
                            s1 = s2;
                        }
                    }
                    const double seed_n_sigma = 2.0;
                    double dataLo = m1 - seed_n_sigma * s1;
                    double dataHi = m1 + seed_n_sigma * s1;
                    dataLo = std::max(dataLo, h->GetXaxis()->GetXmin());
                    dataHi = std::min(dataHi, h->GetXaxis()->GetXmax());
                    if (seedLo < 0) seedLo = dataLo;
                    if (seedHi < 0) seedHi = dataHi;
                }
                return FitUtils::IterativeGausFit(
                    h, tag_prefix, label, scanPdfPath, seedLo, seedHi,
                    /*maxNSigma*/ 2.5,
                    /*physicalLo*/ []() {
                        // Lowest R allowed in the fit, as in Response_Fit_Inclusive.h.
                        // The histograms start at 0.5, so it only matters above that.
                        const char* e = std::getenv("EASYBJETS_RESPONSE_FIT_LO");
                        return (e && *e) ? std::atof(e) : 0.3;
                    }());
            };

            // Data fit: seed from the histogram; its window seeds the MC fits.
            std::string data_tag = Form("data_m%d_b%d_%s", m, bin, wp_label);
            auto fit_data = DoFit(h_data_R_balance_ptbin[m][bin], data_tag.c_str(),
                                  Form("Data  %s  %s", CurrentPtBinLabel(bin), massCutLabels[m]),
                                  plotBase + Form("/Fit_R_DB_data_ptbin%d_scan.pdf", bin));

            double amp_data = fit_data.amp, ampErr_data = fit_data.ampErr;
            double Mean_data = fit_data.mean, meanErr_data = fit_data.meanErr;
            double Sigma_data = fit_data.sigma, sigmaErr_data = fit_data.sigmaErr;
            double chi_data = fit_data.chi2;
            int ndf_data = fit_data.ndf;
            double chi2NDF_data = fit_data.chi2ndf;

            data_means[m][bin] = Mean_data;
            data_mean_errs[m][bin] = meanErr_data;
            data_sigmas[m][bin] = Sigma_data;
            data_sigma_errs[m][bin] = sigmaErr_data;

            // Signal MC fit, seeded with the data window.
            std::string mc_tag = Form("mc_m%d_b%d_%s", m, bin, wp_label);
            auto fit_mc = DoFit(
                h_mc_R_balance_ptbin[m][bin], mc_tag.c_str(),
                Form("MC Z(#rightarrow b#bar{b})+Jets  %s  %s", CurrentPtBinLabel(bin), massCutLabels[m]),
                plotBase + Form("/Fit_R_DB_mc_ptbin%d_scan.pdf", bin), fit_data.fitLo, fit_data.fitHi);

            double amp_mc = fit_mc.amp, ampErr_mc = fit_mc.ampErr;
            double Mean_mc = fit_mc.mean, meanErr_mc = fit_mc.meanErr;
            double Sigma_mc = fit_mc.sigma, sigmaErr_mc = fit_mc.sigmaErr;
            double chi_mc = fit_mc.chi2;
            int ndf_mc = fit_mc.ndf;
            double chi2NDF_mc = fit_mc.chi2ndf;

            mc_means[m][bin] = Mean_mc;
            mc_mean_errs[m][bin] = meanErr_mc;
            mc_sigmas[m][bin] = Sigma_mc;
            mc_sigma_errs[m][bin] = sigmaErr_mc;

            TCanvas* mc_canvas = new TCanvas(Form("mc_canvas_m%d_b%d_%s", m, bin, wp_label),
                                             Form("MC Signal R_response for %s (Mass Cut: %s, WP: %s)",
                                                  CurrentPtBinLabel(bin), massCutLabels[m], wp_label),
                                             600, 600);
            mc_canvas->cd();
            h_mc_R_balance_ptbin[m][bin]->Draw("HIST");
            h_mc_R_balance_ptbin[m][bin]->SetMaximum(h_mc_R_balance_ptbin[m][bin]->GetMaximum() * 1.1);
            h_mc_R_balance_ptbin[m][bin]->SetMinimum(0);
            h_mc_R_balance_ptbin[m][bin]->GetXaxis()->SetTitle(RLab());
            h_mc_R_balance_ptbin[m][bin]->GetYaxis()->SetTitle("Events");

            // The fitted function is drawn over its fit window only.
            TF1* gausFit_mc_inner = FitUtils::MakeModelTF1(fit_mc, Form("gfi_mc_%s", mc_tag.c_str()));
            gausFit_mc_inner->SetLineColor(kViolet - 4);
            gausFit_mc_inner->SetLineStyle(1);
            gausFit_mc_inner->SetLineWidth(3);
            gausFit_mc_inner->Draw("same");

            TLegend* Legend1 = new TLegend(0.62, 0.78, 0.92, 0.92);
            Legend1->SetBorderSize(0);
            Legend1->SetTextSize(0.03);
            Legend1->AddEntry(h_mc_R_balance_ptbin[m][bin], SigLab(), "l");
            Legend1->AddEntry(gausFit_mc_inner, FitUtils::ModelDisplayName(fit_mc.modelName).c_str(), "l");
            TLatex text;
            text.SetTextSize(0.03);
            text.DrawLatexNDC(0.65, 0.75, Form("#mu = %.3f #pm %.3f", Mean_mc, fit_mc.meanErrScaled()));
            text.DrawLatexNDC(0.65, 0.7, Form("#sigma = %.3f #pm %.3f", Sigma_mc, fit_mc.sigmaErrScaled()));
            text.DrawLatexNDC(0.65, 0.65, Form("#chi^{2}/NDF = %.3f", chi2NDF_mc));
            text.DrawLatexNDC(0.65, 0.6, Form("window: [%.3f, %.3f]", fit_mc.fitLo, fit_mc.fitHi));
            Legend1->Draw();
            AtlasUtils::DrawAtlasLabel(0.2, 0.87, 0.15);
            text.DrawLatexNDC(0.2, 0.78, CurrentPtBinLabel(bin));
            text.DrawLatexNDC(0.2, 0.73, Form("Mass: %s", massCutLabels[m]));

            mc_canvas->SaveAs(Form(((RunConfig::PlotDir() + "/%s/Response_Fit_pTbins_Mcuts_")
                                    + CurrentChannel() + "/%s/%s/Fit_R_DB_mc_ptbin%d.pdf")
                                       .c_str(),
                                   cfg.plot_subdir.c_str(), wp_label, fileSystemMassCutLabels[m], bin));
            FitOutputUtils::SaveFitParameters(
                Form(((RunConfig::PlotDir() + "/%s/Response_Fit_pTbins_Mcuts_") + CurrentChannel()
                      + "/%s/%s/mc_fit_parameters_ptbin%d.txt")
                         .c_str(),
                     cfg.plot_subdir.c_str(), wp_label, fileSystemMassCutLabels[m], bin),
                Form("MC: Z (#rightarrow bb) + Jets (%s, Mass Cut: %s)", CurrentPtBinLabel(bin),
                     massCutLabels[m]),
                Mean_mc, fit_mc.meanErrScaled(), Sigma_mc, fit_mc.sigmaErrScaled(), amp_mc,
                fit_mc.ampErrScaled(), chi_mc, ndf_mc);
            if (summaryFile.is_open()) {
                summaryFile << "MC: Z (#rightarrow bb) + Jets (" << CurrentPtBinLabel(bin) << ")"
                            << std::endl;
                summaryFile << "Mean: " << Mean_mc << " +/- " << meanErr_mc << std::endl;
                summaryFile << "Sigma: " << Sigma_mc << " +/- " << sigmaErr_mc << std::endl;
                summaryFile << "Chi2/NDF: " << chi2NDF_mc << " (" << chi_mc << "/" << ndf_mc << ")"
                            << std::endl
                            << std::endl;
            }
            delete mc_canvas;
            delete gausFit_mc_inner;
            delete Legend1;

            // Data fit plot (the fit itself ran above).
            TCanvas* data_canvas = new TCanvas(Form("data_canvas_m%d_b%d_%s", m, bin, wp_label),
                                               Form("Data R_response for %s (Mass Cut: %s, WP: %s)",
                                                    CurrentPtBinLabel(bin), massCutLabels[m], wp_label),
                                               600, 600);
            data_canvas->cd();
            h_data_R_balance_ptbin[m][bin]->Draw("EP");
            h_data_R_balance_ptbin[m][bin]->SetMaximum(h_data_R_balance_ptbin[m][bin]->GetMaximum() * 1.1);
            h_data_R_balance_ptbin[m][bin]->GetXaxis()->SetTitle(RLab());
            h_data_R_balance_ptbin[m][bin]->GetYaxis()->SetTitle("Events");

            TF1* gausFit_data_inner = FitUtils::MakeModelTF1(fit_data, Form("gfi_data_%s", data_tag.c_str()));
            gausFit_data_inner->SetLineColor(kViolet - 4);
            gausFit_data_inner->SetLineStyle(1);
            gausFit_data_inner->SetLineWidth(3);
            gausFit_data_inner->Draw("same");

            TLegend* Legend2 = new TLegend(0.62, 0.78, 0.92, 0.92);
            Legend2->SetBorderSize(0);
            Legend2->SetTextSize(0.03);
            Legend2->AddEntry(h_data_R_balance_ptbin[m][bin], "Data", "lep");
            Legend2->AddEntry(gausFit_data_inner, FitUtils::ModelDisplayName(fit_data.modelName).c_str(),
                              "l");
            TLatex text2;
            text2.SetTextSize(0.03);
            text2.DrawLatexNDC(0.65, 0.75, Form("#mu = %.3f #pm %.3f", Mean_data, fit_data.meanErrScaled()));
            text2.DrawLatexNDC(0.65, 0.7,
                               Form("#sigma = %.3f #pm %.3f", Sigma_data, fit_data.sigmaErrScaled()));
            text2.DrawLatexNDC(0.65, 0.65, Form("#chi^{2}/NDF = %.3f", chi2NDF_data));
            text2.DrawLatexNDC(0.65, 0.6, Form("window: [%.3f, %.3f]", fit_data.fitLo, fit_data.fitHi));
            Legend2->Draw();
            AtlasUtils::DrawAtlasLabel(0.2, 0.87, 0.15);
            text2.DrawLatexNDC(0.2, 0.78, CurrentPtBinLabel(bin));
            text2.DrawLatexNDC(0.2, 0.73, Form("Mass: %s", massCutLabels[m]));

            data_canvas->SaveAs(Form(((RunConfig::PlotDir() + "/%s/Response_Fit_pTbins_Mcuts_")
                                      + CurrentChannel() + "/%s/%s/Fit_R_DB_data_ptbin%d.pdf")
                                         .c_str(),
                                     cfg.plot_subdir.c_str(), wp_label, fileSystemMassCutLabels[m], bin));
            FitOutputUtils::SaveFitParameters(
                Form(((RunConfig::PlotDir() + "/%s/Response_Fit_pTbins_Mcuts_") + CurrentChannel()
                      + "/%s/%s/data_fit_parameters_ptbin%d.txt")
                         .c_str(),
                     cfg.plot_subdir.c_str(), wp_label, fileSystemMassCutLabels[m], bin),
                Form("Data (%s, Mass Cut: %s)", CurrentPtBinLabel(bin), massCutLabels[m]), Mean_data,
                fit_data.meanErrScaled(), Sigma_data, fit_data.sigmaErrScaled(), amp_data,
                fit_data.ampErrScaled(), chi_data, ndf_data);
            if (summaryFile.is_open()) {
                summaryFile << "Data (" << CurrentPtBinLabel(bin) << ")" << std::endl;
                summaryFile << "Mean: " << Mean_data << " +/- " << meanErr_data << std::endl;
                summaryFile << "Sigma: " << Sigma_data << " +/- " << sigmaErr_data << std::endl;
                summaryFile << "Chi2/NDF: " << chi2NDF_data << " (" << chi_data << "/" << ndf_data << ")"
                            << std::endl
                            << std::endl;
            }

            delete data_canvas;
            delete gausFit_data_inner;
            delete Legend2;

            // Background fit, seeded with the data window.
            std::string bkg_tag = Form("bkg_m%d_b%d_%s", m, bin, wp_label);
            auto fit_bkg =
                DoFit(h_multijets_R_balance_ptbin[m][bin], bkg_tag.c_str(),
                      Form("Multijets  %s  %s", CurrentPtBinLabel(bin), massCutLabels[m]),
                      plotBase + Form("/Fit_R_DB_bkg_ptbin%d_scan.pdf", bin), fit_data.fitLo, fit_data.fitHi);

            double amp_bkg = fit_bkg.amp, ampErr_bkg = fit_bkg.ampErr;
            double Mean_bkg = fit_bkg.mean, meanErr_bkg = fit_bkg.meanErr;
            double Sigma_bkg = fit_bkg.sigma, sigmaErr_bkg = fit_bkg.sigmaErr;
            double chi_bkg = fit_bkg.chi2;
            int ndf_bkg = fit_bkg.ndf;
            double chi2NDF_bkg = fit_bkg.chi2ndf;

            bkg_means[m][bin] = Mean_bkg;
            bkg_mean_errs[m][bin] = meanErr_bkg;
            bkg_sigmas[m][bin] = Sigma_bkg;
            bkg_sigma_errs[m][bin] = sigmaErr_bkg;

            TCanvas* bkg_canvas = new TCanvas(Form("bkg_canvas_m%d_b%d_%s", m, bin, wp_label),
                                              Form("Background R_response for %s (Mass Cut: %s, WP: %s)",
                                                   CurrentPtBinLabel(bin), massCutLabels[m], wp_label),
                                              600, 600);
            bkg_canvas->cd();
            h_multijets_R_balance_ptbin[m][bin]->Draw("HIST");
            h_multijets_R_balance_ptbin[m][bin]->SetMaximum(h_multijets_R_balance_ptbin[m][bin]->GetMaximum()
                                                            * 1.1);
            h_multijets_R_balance_ptbin[m][bin]->GetXaxis()->SetTitle(RLab());
            h_multijets_R_balance_ptbin[m][bin]->GetYaxis()->SetTitle("Events");

            TF1* gausFit_bkg_inner = FitUtils::MakeModelTF1(fit_bkg, Form("gfi_bkg_%s", bkg_tag.c_str()));
            gausFit_bkg_inner->SetLineColor(kViolet - 4);
            gausFit_bkg_inner->SetLineStyle(1);
            gausFit_bkg_inner->SetLineWidth(3);
            gausFit_bkg_inner->Draw("same");

            TLegend* Legend3 = new TLegend(0.62, 0.78, 0.92, 0.92);
            Legend3->SetBorderSize(0);
            Legend3->SetTextSize(0.03);
            Legend3->AddEntry(h_multijets_R_balance_ptbin[m][bin], BkgLab(), "l");
            Legend3->AddEntry(gausFit_bkg_inner, FitUtils::ModelDisplayName(fit_bkg.modelName).c_str(), "l");
            TLatex text3;
            text3.SetTextSize(0.03);
            text3.DrawLatexNDC(0.65, 0.75, Form("#mu = %.3f #pm %.3f", Mean_bkg, fit_bkg.meanErrScaled()));
            text3.DrawLatexNDC(0.65, 0.7,
                               Form("#sigma = %.3f #pm %.3f", Sigma_bkg, fit_bkg.sigmaErrScaled()));
            text3.DrawLatexNDC(0.65, 0.65, Form("#chi^{2}/NDF = %.3f", chi2NDF_bkg));
            text3.DrawLatexNDC(0.65, 0.6, Form("window: [%.3f, %.3f]", fit_bkg.fitLo, fit_bkg.fitHi));
            Legend3->Draw();
            AtlasUtils::DrawAtlasLabel(0.2, 0.87, 0.15);
            text3.DrawLatexNDC(0.2, 0.78, CurrentPtBinLabel(bin));
            text3.DrawLatexNDC(0.2, 0.73, Form("Mass: %s", massCutLabels[m]));

            bkg_canvas->SaveAs(Form(((RunConfig::PlotDir() + "/%s/Response_Fit_pTbins_Mcuts_")
                                     + CurrentChannel() + "/%s/%s/Fit_R_DB_bkg_ptbin%d.pdf")
                                        .c_str(),
                                    cfg.plot_subdir.c_str(), wp_label, fileSystemMassCutLabels[m], bin));
            FitOutputUtils::SaveFitParameters(
                Form(((RunConfig::PlotDir() + "/%s/Response_Fit_pTbins_Mcuts_") + CurrentChannel()
                      + "/%s/%s/bkg_fit_parameters_ptbin%d.txt")
                         .c_str(),
                     cfg.plot_subdir.c_str(), wp_label, fileSystemMassCutLabels[m], bin),
                Form("Background: Multijets (%s, Mass Cut: %s)", CurrentPtBinLabel(bin), massCutLabels[m]),
                Mean_bkg, fit_bkg.meanErrScaled(), Sigma_bkg, fit_bkg.sigmaErrScaled(), amp_bkg,
                fit_bkg.ampErrScaled(), chi_bkg, ndf_bkg);
            if (summaryFile.is_open()) {
                summaryFile << "Background: Multijets (" << CurrentPtBinLabel(bin) << ")" << std::endl;
                summaryFile << "Mean: " << Mean_bkg << " +/- " << meanErr_bkg << std::endl;
                summaryFile << "Sigma: " << Sigma_bkg << " +/- " << sigmaErr_bkg << std::endl;
                summaryFile << "Chi2/NDF: " << chi2NDF_bkg << " (" << chi_bkg << "/" << ndf_bkg << ")"
                            << std::endl
                            << std::endl;
            }
            delete bkg_canvas;
            delete gausFit_bkg_inner;
            delete Legend3;

            // Combined MC fit, seeded with the data window. Its width is the MC
            // width used in the JER closure and in R_JER.
            std::string comb_tag = Form("combined_m%d_b%d_%s", m, bin, wp_label);
            auto fit_combined = DoFit(h_combined_MC_R_balance_ptbin[m][bin], comb_tag.c_str(),
                                      Form("Combined MC  %s  %s", CurrentPtBinLabel(bin), massCutLabels[m]),
                                      plotBase + Form("/Fit_R_DB_combined_MC_ptbin%d_scan.pdf", bin),
                                      fit_data.fitLo, fit_data.fitHi);

            double amp_combined = fit_combined.amp, ampErr_combined = fit_combined.ampErr;
            double Mean_combined = fit_combined.mean, meanErr_combined = fit_combined.meanErr;
            double Sigma_combined = fit_combined.sigma, sigmaErr_combined = fit_combined.sigmaErr;
            double chi_combined = fit_combined.chi2;
            int ndf_combined = fit_combined.ndf;
            double chi2NDF_combined = fit_combined.chi2ndf;

            combined_means[m][bin] = Mean_combined;
            combined_mean_errs[m][bin] = meanErr_combined;
            combined_sigmas[m][bin] = Sigma_combined;
            combined_sigma_errs[m][bin] = sigmaErr_combined;

            TCanvas* combined_canvas =
                new TCanvas(Form("combined_canvas_m%d_b%d_%s", m, bin, wp_label),
                            Form(("Combined MC " + cfg.branch_R + " for %s (Mass Cut: %s, WP: %s)").c_str(),
                                 CurrentPtBinLabel(bin), massCutLabels[m], wp_label),
                            600, 600);
            combined_canvas->cd();
            h_combined_MC_R_balance_ptbin[m][bin]->Draw("HIST");
            h_combined_MC_R_balance_ptbin[m][bin]->SetMaximum(
                h_combined_MC_R_balance_ptbin[m][bin]->GetMaximum() * 1.1);
            h_combined_MC_R_balance_ptbin[m][bin]->GetXaxis()->SetTitle(RLab());
            h_combined_MC_R_balance_ptbin[m][bin]->GetYaxis()->SetTitle("Events");

            TF1* gausFit_combined_inner =
                FitUtils::MakeModelTF1(fit_combined, Form("gfi_comb_%s", comb_tag.c_str()));
            gausFit_combined_inner->SetLineColor(kViolet - 4);
            gausFit_combined_inner->SetLineStyle(1);
            gausFit_combined_inner->SetLineWidth(3);
            gausFit_combined_inner->Draw("same");

            TLegend* combined_legend = new TLegend(0.62, 0.78, 0.92, 0.92);
            combined_legend->SetBorderSize(0);
            combined_legend->SetTextSize(0.03);
            combined_legend->AddEntry(h_combined_MC_R_balance_ptbin[m][bin], "Combined MC", "l");
            combined_legend->AddEntry(gausFit_combined_inner,
                                      FitUtils::ModelDisplayName(fit_combined.modelName).c_str(), "l");
            TLatex combined_text;
            combined_text.SetTextSize(0.03);
            combined_text.DrawLatexNDC(
                0.65, 0.75, Form("#mu = %.3f #pm %.3f", Mean_combined, fit_combined.meanErrScaled()));
            combined_text.DrawLatexNDC(
                0.65, 0.7, Form("#sigma = %.3f #pm %.3f", Sigma_combined, fit_combined.sigmaErrScaled()));
            combined_text.DrawLatexNDC(0.65, 0.65, Form("#chi^{2}/NDF = %.3f", chi2NDF_combined));
            combined_text.DrawLatexNDC(0.65, 0.6,
                                       Form("window: [%.3f, %.3f]", fit_combined.fitLo, fit_combined.fitHi));
            combined_legend->Draw();
            AtlasUtils::DrawAtlasLabel(0.2, 0.87, 0.15);
            combined_text.DrawLatexNDC(0.2, 0.78, CurrentPtBinLabel(bin));
            combined_text.DrawLatexNDC(0.2, 0.73, Form("Mass: %s", massCutLabels[m]));

            combined_canvas->SaveAs(Form(((RunConfig::PlotDir() + "/%s/Response_Fit_pTbins_Mcuts_")
                                          + CurrentChannel() + "/%s/%s/Fit_R_DB_combined_MC_ptbin%d.pdf")
                                             .c_str(),
                                         cfg.plot_subdir.c_str(), wp_label, fileSystemMassCutLabels[m], bin));
            FitOutputUtils::SaveFitParameters(
                Form(((RunConfig::PlotDir() + "/%s/Response_Fit_pTbins_Mcuts_") + CurrentChannel()
                      + "/%s/%s/combined_MC_fit_parameters_ptbin%d.txt")
                         .c_str(),
                     cfg.plot_subdir.c_str(), wp_label, fileSystemMassCutLabels[m], bin),
                Form("Combined MC (Signal+Bkg) (%s, Mass Cut: %s)", CurrentPtBinLabel(bin), massCutLabels[m]),
                Mean_combined, fit_combined.meanErrScaled(), Sigma_combined, fit_combined.sigmaErrScaled(),
                amp_combined, fit_combined.ampErrScaled(), chi_combined, ndf_combined);
            if (summaryFile.is_open()) {
                summaryFile << "Combined MC (Signal+Bkg) (" << CurrentPtBinLabel(bin) << ")" << std::endl;
                summaryFile << "Mean: " << Mean_combined << " +/- " << meanErr_combined << std::endl;
                summaryFile << "Sigma: " << Sigma_combined << " +/- " << sigmaErr_combined << std::endl;
                summaryFile << "Chi2/NDF: " << chi2NDF_combined << " (" << chi_combined << "/" << ndf_combined
                            << ")" << std::endl
                            << std::endl;
            }
            delete combined_canvas;
            delete gausFit_combined_inner;
            delete combined_legend;

            // Reco widths are stored even when the truth fits below fail.
            jer_sigma_RDB_data[m][bin] = Sigma_data;
            jer_sigma_RDB_data_err[m][bin] = sigmaErr_data;
            jer_sigma_RDB_MC[m][bin] = Sigma_mc;
            jer_sigma_RDB_MC_err[m][bin] = sigmaErr_mc;

            // JER. Data contains signal and background, so the MC
            // side uses the combined MC to match its composition.
            if (hasTruthBranches) {
                // Combined truth/ref, with the same DoFit as the reco fits (seed
                // from the histogram, no rebinning), so its scan PDF shows every
                // candidate window like the reco scans do.
                std::string tr_tag = Form("truth_m%d_b%d_%s", m, bin, wp_label);
                auto fit_tr =
                    DoFit(h_comb_truth_R_balance_ptbin[m][bin], tr_tag.c_str(), "Combined truth/ref",
                          plotBase + Form("/Fit_R_truth_ptbin%d_scan.pdf", bin));

                // Combined pT_reco / pT_truth (sigma^MC)
                std::string rt_tag = Form("rectruth_m%d_b%d_%s", m, bin, wp_label);
                auto fit_rt = DoFit(h_comb_reco_truth_ptbin[m][bin], rt_tag.c_str(), "Combined reco/truth",
                                    plotBase + Form("/Fit_R_recotruth_ptbin%d_scan.pdf", bin));

                if (fit_tr.valid && fit_rt.valid) {
                    double s_data = Sigma_data, ds_data = sigmaErr_data;
                    // Combined MC width, to match the data composition
                    double s_mc = Sigma_combined, ds_mc = sigmaErr_combined;
                    double s_truth = fit_tr.sigma, ds_truth = fit_tr.sigmaErr;
                    double s_rt = fit_rt.sigma, ds_rt = fit_rt.sigmaErr;

                    // In low-statistics bins the iterative fit can end on a
                    // meaningless solution. Accept only 0.5 < mean < 1.5 and
                    // 0.01 < sigma < 0.20 for both truth fits; otherwise skip the
                    // closure for this cell.
                    bool truth_sane =
                        (fit_tr.mean > 0.5 && fit_tr.mean < 1.5 && s_truth > 0.01 && s_truth < 0.20);
                    bool reco_sane = (fit_rt.mean > 0.5 && fit_rt.mean < 1.5 && s_rt > 0.01 && s_rt < 0.20);

                    if (!truth_sane || !reco_sane) {
                        Log::warning(("Rdb" + CurrentChannel()).c_str(),
                                     Form("[m%d b%d %s] Pathological truth fit skipped: "
                                          "truth(mu=%.3f sigma=%.4f) reco(mu=%.3f sigma=%.4f). "
                                          "Likely low stats; no closure for this bin.",
                                          m, bin, wp_label, fit_tr.mean, s_truth, fit_rt.mean, s_rt));
                        continue;
                    }

                    //   A = sigma^MC        = s_rt (reco/truth)
                    //   B = sigma^MC_insitu = sqrt(x1^2 - x3^2)
                    //   C = sigma_insitu    = sqrt(x2^2 - x3^2)
                    //   x1 = sigma(R) combined MC, x2 = sigma(R) data, x3 = sigma_truth
                    //   Delta_closure = (A - B) / C
                    // x3 enters both B and C, so its error is propagated once through
                    // both instead of treating B and C as independent.
                    double sig_MC = s_rt; // A
                    double dsig_MC = ds_rt;
                    double d_in = s_data * s_data - s_truth * s_truth;
                    double d_mcin = s_mc * s_mc - s_truth * s_truth;
                    double sig_in = (d_in > 0) ? std::sqrt(d_in) : 0.0;       // C
                    double sig_mcin = (d_mcin > 0) ? std::sqrt(d_mcin) : 0.0; // B

                    // Error on sigma_insitu (independent inputs)
                    double dsig_in = (sig_in > 0) ? std::sqrt(std::pow(s_data * ds_data / sig_in, 2)
                                                              + std::pow(s_truth * ds_truth / sig_in, 2))
                                                  : 0.0;

                    // Error on sigma^MC_insitu (independent inputs)
                    double dsig_mcin = (sig_mcin > 0)
                                           ? std::sqrt(std::pow(s_mc * ds_mc / sig_mcin, 2)
                                                       + std::pow(s_truth * ds_truth / sig_mcin, 2))
                                           : 0.0;

                    // Delta_closure, with sigma_truth shared between B and C
                    double delta_clos = 0, ddelta_clos = 0;
                    if (sig_in > 1e-9 && sig_mcin > 1e-9) {
                        const double A = sig_MC;
                        const double B = sig_mcin;
                        const double C = sig_in;
                        const double x1 = s_mc;
                        const double sx1 = ds_mc;
                        const double x2 = s_data;
                        const double sx2 = ds_data;
                        const double x3 = s_truth;
                        const double sx3 = ds_truth;
                        const double sA = dsig_MC;
                        const double AmB = A - B;
                        delta_clos = AmB / C;

                        // Partials:  dDelta/dA  = 1/C
                        //            dDelta/dx1 = d(-B/C)/dx1 = -(1/C) · (x1/B)
                        //                                  (B depends on x1; C does not)
                        //            dDelta/dx2 = d(A-B)/C / dx2 = -(A-B)/C^2 · (x2/C)
                        //                                  (C depends on x2; A,B do not)
                        //            dDelta/dx3 = shared term - x3 in both B and C:
                        //                    (A-B)/C^2 · (x3/C) + (1/C) · (x3/B)
                        //                  = (x3/C) · [(A-B)/C^2 + 1/B]
                        const double dA = 1.0 / C;
                        const double dx1 = -(1.0 / C) * (x1 / B);
                        const double dx2 = -AmB * x2 / (C * C * C);
                        const double dx3 = (x3 / C) * ((AmB / (C * C)) + (1.0 / B));
                        ddelta_clos = std::sqrt(dA * dA * sA * sA + dx1 * dx1 * sx1 * sx1
                                                + dx2 * dx2 * sx2 * sx2 + dx3 * dx3 * sx3 * sx3);
                    }

                    Log::info(
                        ("Rdb" + CurrentChannel()).c_str(),
                        Form(
                            "JER [m%d b%d %s]: sigma_in=%.4f±%.4f sigma^MC_in=%.4f±%.4f sigma^MC=%.4f±%.4f Deltaclos=%.3f±%.3f",
                            m, bin, wp_label, sig_in, dsig_in, sig_mcin, dsig_mcin, sig_MC, dsig_MC,
                            delta_clos, ddelta_clos));

                    // Store in arrays for summary plots
                    jer_sigma_RDB_data[m][bin] = Sigma_data;
                    jer_sigma_RDB_data_err[m][bin] = sigmaErr_data;
                    jer_sigma_RDB_MC[m][bin] = Sigma_mc;
                    jer_sigma_RDB_MC_err[m][bin] = sigmaErr_mc;
                    jer_sigma_truth[m][bin] = s_truth;
                    jer_sigma_truth_err[m][bin] = ds_truth;
                    jer_sigma_reco_truth[m][bin] = s_rt;
                    jer_sigma_reco_truth_err[m][bin] = ds_rt;
                    jer_sigma_insitu[m][bin] = sig_in;
                    jer_sigma_insitu_err[m][bin] = dsig_in;
                    jer_sigma_MC_insitu[m][bin] = sig_mcin;
                    jer_sigma_MC_insitu_err[m][bin] = dsig_mcin;
                    jer_sigma_MC[m][bin] = sig_MC;
                    jer_sigma_MC_err[m][bin] = dsig_MC;
                    // b-JES data/MC ratio per (mass cut, pT bin)
                    auto jes_r =
                        FitUtils::ComputeJESRatio(Mean_data, meanErr_data, Mean_combined, meanErr_combined);
                    jes_ratio[m][bin] = jes_r.ratio;
                    jes_ratio_err[m][bin] = jes_r.ratioErr;
                    // b-JER data/MC ratio per (mass cut, pT bin)
                    auto jer_r = FitUtils::ComputeJERRatio(s_data, ds_data, s_mc, ds_mc, s_truth, ds_truth);
                    jer_ratio[m][bin] = jer_r.ratio;
                    jer_ratio_err[m][bin] = jer_r.ratioErr;
                    jer_delta_closure[m][bin] = delta_clos;
                    jer_delta_closure_err[m][bin] = ddelta_clos;
                    jer_valid[m][bin] = true;

                    // Save truth fit parameters (combined - used for JER)
                    FitOutputUtils::SaveFitParameters(
                        Form("%s/truth_ref_fit_parameters_ptbin%d.txt", plotBase.c_str(), bin),
                        Form("Combined truth/ref (%s, %s)", CurrentPtBinLabel(bin), massCutLabels[m]),
                        fit_tr.mean, fit_tr.meanErr, fit_tr.sigma, fit_tr.sigmaErr, fit_tr.amp, fit_tr.ampErr,
                        fit_tr.chi2, fit_tr.ndf);
                    FitOutputUtils::SaveFitParameters(
                        Form("%s/reco_truth_fit_parameters_ptbin%d.txt", plotBase.c_str(), bin),
                        Form("Combined reco/truth (%s, %s)", CurrentPtBinLabel(bin), massCutLabels[m]),
                        fit_rt.mean, fit_rt.meanErr, fit_rt.sigma, fit_rt.sigmaErr, fit_rt.amp, fit_rt.ampErr,
                        fit_rt.chi2, fit_rt.ndf);

                    // Truth fit plots, in the style of the reco fit plots, to check
                    // that the chosen window sits on the peak.
                    auto DrawTruthFitOverlay = [&](TH1F* h, const FitUtils::GausFitResult& r,
                                                   const char* legendEntry, const char* xAxisTitle,
                                                   const char* savePath) {
                        if (!h || !r.valid) return;
                        TCanvas* c = new TCanvas(Form("c_truth_%s_m%d_b%d_%s", legendEntry, m, bin, wp_label),
                                                 "truth fit overlay", 600, 600);
                        c->cd();

                        h->Draw("HIST");
                        h->SetLineColor(kBlack);
                        h->SetLineWidth(2);
                        h->SetFillStyle(0);
                        h->SetMaximum(h->GetMaximum() * 1.35);
                        h->SetMinimum(0);
                        h->GetXaxis()->SetTitle(xAxisTitle);
                        h->GetYaxis()->SetTitle("Events");

                        // Fitted function over its fit window only
                        TF1* gfInner = FitUtils::MakeModelTF1(
                            r, Form("gfi_truth_%s_m%d_b%d_%s", legendEntry, m, bin, wp_label));
                        gfInner->SetLineColor(kViolet - 4);
                        gfInner->SetLineStyle(1);
                        gfInner->SetLineWidth(3);
                        gfInner->Draw("same");

                        TLegend* leg = new TLegend(0.62, 0.78, 0.92, 0.92);
                        leg->SetBorderSize(0);
                        leg->SetTextSize(0.03);
                        leg->AddEntry(h, legendEntry, "l");
                        leg->AddEntry(gfInner, "Gaussian Fit", "l");
                        leg->Draw();

                        TLatex lt;
                        lt.SetTextSize(0.03);
                        lt.SetNDC();
                        lt.DrawLatexNDC(0.65, 0.75, Form("#mu = %.4f #pm %.4f", r.mean, r.meanErrScaled()));
                        lt.DrawLatexNDC(0.65, 0.7,
                                        Form("#sigma = %.4f #pm %.4f", r.sigma, r.sigmaErrScaled()));
                        lt.DrawLatexNDC(0.65, 0.65, Form("#chi^{2}/NDF = %.3f", r.chi2ndf));
                        lt.DrawLatexNDC(0.65, 0.6, Form("window: [%.3f, %.3f]", r.fitLo, r.fitHi));

                        AtlasUtils::DrawAtlasLabel(0.2, 0.87, 0.15);
                        lt.DrawLatexNDC(0.2, 0.78, CurrentPtBinLabel(bin));
                        lt.DrawLatexNDC(0.2, 0.73, Form("Mass: %s", massCutLabels[m]));

                        c->SaveAs(savePath);
                        delete gfInner;
                        delete leg;
                        delete c;
                    };

                    DrawTruthFitOverlay(h_comb_truth_R_balance_ptbin[m][bin], fit_tr, "Combined MC truth/ref",
                                        "p_{T}^{truth} / p_{T}^{ref,DB}",
                                        Form("%s/Fit_R_truth_ptbin%d.pdf", plotBase.c_str(), bin));

                    DrawTruthFitOverlay(h_comb_reco_truth_ptbin[m][bin], fit_rt, "Combined MC reco/truth",
                                        "p_{T}^{reco} / p_{T}^{truth}",
                                        Form("%s/Fit_R_recotruth_ptbin%d.pdf", plotBase.c_str(), bin));

                    // Signal-only and background-only truth plots, to see why the
                    // combined result differs. These fits use the IterativeGausFit
                    // defaults (max 2 sigma, no R floor), unlike DoFit.
                    if (h_mc_truth_R_balance_ptbin[m][bin]
                        && h_mc_truth_R_balance_ptbin[m][bin]->GetEntries() > 20) {
                        auto [sTrLo, sTrHi] = DataDrivenSeed(h_mc_truth_R_balance_ptbin[m][bin]);
                        FitUtils::GausFitResult r_sig_tr =
                            FitUtils::IterativeGausFit(h_mc_truth_R_balance_ptbin[m][bin],
                                                       Form("sig_truth_overlay_m%d_b%d_%s", m, bin, wp_label),
                                                       "Signal truth/ref", "", sTrLo, sTrHi);
                        if (r_sig_tr.valid)
                            DrawTruthFitOverlay(
                                h_mc_truth_R_balance_ptbin[m][bin], r_sig_tr, "Signal truth/ref",
                                "p_{T}^{truth} / p_{T}^{ref,DB}",
                                Form("%s/Fit_R_truth_signal_ptbin%d.pdf", plotBase.c_str(), bin));
                    }
                    if (h_mc_reco_truth_ptbin[m][bin] && h_mc_reco_truth_ptbin[m][bin]->GetEntries() > 20) {
                        auto [sRtLo, sRtHi] = DataDrivenSeed(h_mc_reco_truth_ptbin[m][bin]);
                        FitUtils::GausFitResult r_sig_rt = FitUtils::IterativeGausFit(
                            h_mc_reco_truth_ptbin[m][bin],
                            Form("sig_recotruth_overlay_m%d_b%d_%s", m, bin, wp_label), "Signal reco/truth",
                            "", sRtLo, sRtHi);
                        if (r_sig_rt.valid)
                            DrawTruthFitOverlay(
                                h_mc_reco_truth_ptbin[m][bin], r_sig_rt, "Signal reco/truth",
                                "p_{T}^{reco} / p_{T}^{truth}",
                                Form("%s/Fit_R_recotruth_signal_ptbin%d.pdf", plotBase.c_str(), bin));
                    }
                    if (hasTruthBkg && h_bkg_truth_R_balance_ptbin[m][bin]
                        && h_bkg_truth_R_balance_ptbin[m][bin]->GetEntries() > 20) {
                        auto [bTrLo, bTrHi] = DataDrivenSeed(h_bkg_truth_R_balance_ptbin[m][bin]);
                        FitUtils::GausFitResult r_bkg_tr =
                            FitUtils::IterativeGausFit(h_bkg_truth_R_balance_ptbin[m][bin],
                                                       Form("bkg_truth_overlay_m%d_b%d_%s", m, bin, wp_label),
                                                       "Bkg truth/ref", "", bTrLo, bTrHi);
                        if (r_bkg_tr.valid)
                            DrawTruthFitOverlay(
                                h_bkg_truth_R_balance_ptbin[m][bin], r_bkg_tr, Form("%s truth/ref", BkgLab()),
                                "p_{T}^{truth} / p_{T}^{ref,DB}",
                                Form("%s/Fit_R_truth_bkg_ptbin%d.pdf", plotBase.c_str(), bin));
                    }
                    if (hasTruthBkg && h_bkg_reco_truth_ptbin[m][bin]
                        && h_bkg_reco_truth_ptbin[m][bin]->GetEntries() > 20) {
                        auto [bRtLo, bRtHi] = DataDrivenSeed(h_bkg_reco_truth_ptbin[m][bin]);
                        FitUtils::GausFitResult r_bkg_rt = FitUtils::IterativeGausFit(
                            h_bkg_reco_truth_ptbin[m][bin],
                            Form("bkg_recotruth_overlay_m%d_b%d_%s", m, bin, wp_label), "Bkg reco/truth", "",
                            bRtLo, bRtHi);
                        if (r_bkg_rt.valid)
                            DrawTruthFitOverlay(
                                h_bkg_reco_truth_ptbin[m][bin], r_bkg_rt, Form("%s reco/truth", BkgLab()),
                                "p_{T}^{reco} / p_{T}^{truth}",
                                Form("%s/Fit_R_recotruth_bkg_ptbin%d.pdf", plotBase.c_str(), bin));
                    }

                    // Save JER summary
                    {
                        std::ofstream jer_txt(Form("%s/JER_results_ptbin%d.txt", plotBase.c_str(), bin));
                        jer_txt << "JER Results -- " << cfg.key << " WP " << wp_label << " "
                                << CurrentPtBinLabel(bin) << " " << massCutLabels[m] << "\n";
                        jer_txt << std::fixed << std::setprecision(6) << "sigma_insitu:    " << sig_in
                                << " +/- " << dsig_in << "\n"
                                << "sigma_MC_insitu:      " << sig_mcin << " +/- " << dsig_mcin << "\n"
                                << "sigma_MC:        " << sig_MC << " +/- " << dsig_MC << "\n"
                                << "delta_closure:   " << delta_clos << " +/- " << ddelta_clos << "\n"
                                << "R_JER (data/MC):      " << jer_r.ratio << " +/- " << jer_r.ratioErr
                                << "\n";
                        jer_txt.close();
                    }
                } else {
                    Log::warning(("Rdb" + CurrentChannel()).c_str(),
                                 Form("JER truth fits failed for m%d b%d - skipping", m, bin));
                }
            }

            // Per-sample JER breakdown (signal, background, combined), as in
            // Response_Fit_Inclusive.h.
            if (hasTruthBranches && jer_valid[m][bin]) {
                // PerSampFit and ComputePerSampleJER are in Response_Fit_Core.h.

                // Per-sample truth fits
                auto QuickFit = [&](TH1F* h, const std::string& tag,
                                    const std::string& label) -> std::pair<double, double> {
                    if (!h || h->GetEntries() < 20) return {0, 0};
                    // Same DoFit as the reco fits
                    auto r = DoFit(h, tag.c_str(), label.c_str(), plotBase + "/" + tag + "_scan.pdf");
                    if (!r.valid) return {0, 0};
                    FitOutputUtils::SaveFitParameters(
                        Form("%s/%s_fit_parameters_ptbin%d.txt", plotBase.c_str(), tag.c_str(), bin),
                        label.c_str(), r.mean, r.meanErr, r.sigma, r.sigmaErr, r.amp, r.ampErr, r.chi2,
                        r.ndf);
                    return {r.sigma, r.sigmaErr};
                };

                auto [st_sig, ste_sig] =
                    QuickFit(h_mc_truth_R_balance_ptbin[m][bin],
                             Form("truth_sig_m%d_b%d_%s", m, bin, wp_label), "Sig truth/ref");
                auto [sr_sig, sre_sig] =
                    QuickFit(h_mc_reco_truth_ptbin[m][bin],
                             Form("reco_truth_sig_m%d_b%d_%s", m, bin, wp_label), "Sig reco/truth");
                auto [st_bkg, ste_bkg] =
                    QuickFit(h_bkg_truth_R_balance_ptbin[m][bin],
                             Form("truth_bkg_m%d_b%d_%s", m, bin, wp_label), "Bkg truth/ref");
                auto [sr_bkg, sre_bkg] =
                    QuickFit(h_bkg_reco_truth_ptbin[m][bin],
                             Form("reco_truth_bkg_m%d_b%d_%s", m, bin, wp_label), "Bkg reco/truth");
                // Combined already fitted above (fit_tr, fit_rt)
                double st_comb = jer_sigma_truth[m][bin], ste_comb = jer_sigma_truth_err[m][bin];
                double sr_comb = jer_sigma_reco_truth[m][bin], sre_comb = jer_sigma_reco_truth_err[m][bin];

                PerSampFit psf[3];
                psf[0] = ComputePerSampleJER(Sigma_data, sigmaErr_data, mc_sigmas[m][bin],
                                             mc_sigma_errs[m][bin], st_sig, ste_sig, sr_sig, sre_sig);
                psf[1] = ComputePerSampleJER(Sigma_data, sigmaErr_data, bkg_sigmas[m][bin],
                                             bkg_sigma_errs[m][bin], st_bkg, ste_bkg, sr_bkg, sre_bkg);
                psf[2] =
                    ComputePerSampleJER(Sigma_data, sigmaErr_data, combined_sigmas[m][bin],
                                        combined_sigma_errs[m][bin], st_comb, ste_comb, sr_comb, sre_comb);
                const char* psLabels[3] = {"Signal", "Bkg", "Combined"};

                // 6-marker sigma summary plot
                {
                    TCanvas* csum =
                        new TCanvas(Form("csum_m%d_b%d_%s", m, bin, wp_label), "sigma summary", 900, 600);
                    csum->SetLeftMargin(0.12);
                    csum->SetBottomMargin(0.14);
                    csum->SetGridy();
                    double ymax = 0;
                    for (int i = 0; i < 3; i++)
                        ymax = std::max(
                            {ymax, psf[i].sig_RDB + psf[i].sig_RDB_err,
                             psf[i].sig_RDB_MC + psf[i].sig_RDB_MC_err,
                             psf[i].sig_truth + psf[i].sig_truth_err, psf[i].sig_MC + psf[i].sig_MC_err,
                             psf[i].sig_insitu + psf[i].sig_insitu_err, psf[i].sig_MCi + psf[i].sig_MCi_err});
                    TH1F* fr = new TH1F(Form("fr_sum_m%d_b%d_%s", m, bin, wp_label), ";;#sigma", 3, 0, 2);
                    fr->SetMinimum(-0.05 * ymax);
                    fr->SetMaximum(ymax * 1.3);
                    fr->GetXaxis()->SetBinLabel(1, "Signal");
                    fr->GetXaxis()->SetBinLabel(2, "Background");
                    fr->GetXaxis()->SetBinLabel(3, "Combined MC");
                    fr->GetXaxis()->SetLabelSize(0.05);
                    fr->SetStats(0);
                    fr->Draw("AXIS");
                    int mkS[6] = {20, 21, 22, 23, 29, 33};
                    int mkC[6] = {kBlack, kRed, kAzure + 2, kOrange + 7, kViolet - 1, kGreen + 2};
                    std::string qLd = std::string("#sigma(") + RLab() + ")_{data}",
                                qLm = std::string("#sigma(") + RLab() + ")_{MC}";
                    const char* qL[6] = {qLd.c_str(),
                                         qLm.c_str(),
                                         "#sigma(R_{truth})_{MC}",
                                         "#sigma^{MC}(reco/truth)",
                                         "#sigma_{insitu}",
                                         "#sigma^{MC}_{insitu}"};
                    double sOff[3] = {0.333, 1.0, 1.5};
                    TLegend* leg = new TLegend(0.72, 0.60, 0.92, 0.92);
                    leg->SetBorderSize(0);
                    leg->SetFillStyle(0);
                    leg->SetTextSize(0.022);
                    for (int q = 0; q < 6; q++) {
                        double v[3], e[3];
                        for (int i = 0; i < 3; i++) {
                            switch (q) {
                            case 0:
                                v[i] = psf[i].sig_RDB;
                                e[i] = psf[i].sig_RDB_err;
                                break;
                            case 1:
                                v[i] = psf[i].sig_RDB_MC;
                                e[i] = psf[i].sig_RDB_MC_err;
                                break;
                            case 2:
                                v[i] = psf[i].sig_truth;
                                e[i] = psf[i].sig_truth_err;
                                break;
                            case 3:
                                v[i] = psf[i].sig_MC;
                                e[i] = psf[i].sig_MC_err;
                                break;
                            case 4:
                                v[i] = psf[i].sig_insitu;
                                e[i] = psf[i].sig_insitu_err;
                                break;
                            case 5:
                                v[i] = psf[i].sig_MCi;
                                e[i] = psf[i].sig_MCi_err;
                                break;
                            }
                        }
                        double xp[3], xe[3] = {0, 0, 0};
                        for (int i = 0; i < 3; i++) {
                            xp[i] = sOff[i] + (q - 2.5) * 0.04;
                        }
                        TGraphErrors* gr = new TGraphErrors(3, xp, v, xe, e);
                        gr->SetMarkerStyle(mkS[q]);
                        gr->SetMarkerColor(mkC[q]);
                        gr->SetMarkerSize(1.3);
                        gr->SetLineColor(mkC[q]);
                        gr->SetLineWidth(2);
                        gr->Draw("P same");
                        leg->AddEntry(gr, qL[q], "lp");
                    }
                    leg->Draw();
                    AtlasUtils::DrawAtlasLabel(0.18, 0.85);
                    TLatex lt;
                    lt.SetNDC();
                    lt.SetTextSize(0.028);
                    lt.DrawLatex(0.18, 0.75,
                                 Form("%s  WP: %s  %s  %s", cfg.label.c_str(), wp_label, massCutLabels[m],
                                      CurrentPtBinLabel(bin)));
                    csum->SaveAs(Form("%s/JER_sigma_summary_markers_ptbin%d.pdf", plotBase.c_str(), bin));
                    delete csum;
                    delete fr;
                }

                // 3-point JER summary with R_JES, R_JER and Delta_closure (same plot as Response_Fit_Inclusive.h)
                {
                    TCanvas* c = new TCanvas(Form("jer3_m%d_b%d_%s", m, bin, wp_label), "b-JER", 800, 600);
                    c->SetLeftMargin(0.15);

                    double vals[3] = {jer_sigma_insitu[m][bin], jer_sigma_MC_insitu[m][bin],
                                      jer_sigma_MC[m][bin]};
                    double errs[3] = {jer_sigma_insitu_err[m][bin], jer_sigma_MC_insitu_err[m][bin],
                                      jer_sigma_MC_err[m][bin]};
                    double xp[3] = {1, 2, 3}, xe[3] = {0, 0, 0};

                    TGraphErrors* g = new TGraphErrors(3, xp, vals, xe, errs);
                    g->SetTitle(
                        Form("b-JER (%s %s %s);;#sigma", cfg.key.c_str(), wp_label, CurrentPtBinLabel(bin)));
                    g->SetMarkerStyle(21);
                    g->SetMarkerSize(1.5);
                    g->SetLineWidth(2);
                    g->GetXaxis()->SetBinLabel(g->GetXaxis()->FindBin(1), "#sigma_{in situ}");
                    g->GetXaxis()->SetBinLabel(g->GetXaxis()->FindBin(2), "#sigma^{MC}_{in situ}");
                    g->GetXaxis()->SetBinLabel(g->GetXaxis()->FindBin(3), "#sigma^{MC}");
                    g->GetXaxis()->SetLabelSize(0.06);

                    double maxY =
                        *std::max_element(vals, vals + 3) + *std::max_element(errs, errs + 3) + 0.05;
                    g->GetYaxis()->SetRangeUser(0, maxY > 0.05 ? maxY : 0.5);
                    g->Draw("AP");

                    TLatex lt;
                    lt.SetNDC();
                    lt.SetTextFont(42);
                    lt.SetTextSize(0.035);
                    lt.DrawLatex(0.18, 0.80, Form("%s  %s", massCutLabels[m], CurrentPtBinLabel(bin)));
                    lt.DrawLatex(0.18, 0.75,
                                 Form("#Delta#sigma_{cl}=%.4f#pm%.4f", jer_delta_closure[m][bin],
                                      jer_delta_closure_err[m][bin]));
                    if (jes_ratio[m][bin] != 0)
                        lt.DrawLatex(0.18, 0.70,
                                     Form("R_{JES}=%.4f#pm%.4f", jes_ratio[m][bin], jes_ratio_err[m][bin]));
                    lt.DrawLatex(0.18, 0.65,
                                 Form("R_{JER}=%.4f#pm%.4f", jer_ratio[m][bin], jer_ratio_err[m][bin]));

                    AtlasUtils::DrawAtlasLabel(0.2, 0.88);
                    c->SaveAs(Form("%s/JER_summary_ptbin%d.pdf", plotBase.c_str(), bin));
                    delete c;
                }

                // Per-sample R_response overlay (reco, truth, reco/truth)
                {
                    const char* sLab[3] = {"Signal", "Bkg", "Combined"};
                    const char* sTag[3] = {"sig", "bkg", "comb"};
                    TH1F* hRDB[3] = {h_mc_R_balance_ptbin[m][bin], h_multijets_R_balance_ptbin[m][bin],
                                     h_combined_MC_R_balance_ptbin[m][bin]};
                    TH1F* hTr[3] = {h_mc_truth_R_balance_ptbin[m][bin], h_bkg_truth_R_balance_ptbin[m][bin],
                                    h_comb_truth_R_balance_ptbin[m][bin]};
                    TH1F* hRT[3] = {h_mc_reco_truth_ptbin[m][bin], h_bkg_reco_truth_ptbin[m][bin],
                                    h_comb_reco_truth_ptbin[m][bin]};
                    for (int s = 0; s < 3; s++) {
                        if (hRDB[s]->GetEntries() < 10) continue;
                        TCanvas* c = new TCanvas(Form("c_3ov_%s_m%d_b%d_%s", sTag[s], m, bin, wp_label),
                                                 sLab[s], 800, 600);
                        TH1F* hr = (TH1F*)hRDB[s]->Clone();
                        TH1F* ht = (TH1F*)hTr[s]->Clone();
                        TH1F* hrt = (TH1F*)hRT[s]->Clone();
                        if (hr->Integral() > 0) hr->Scale(1. / hr->Integral());
                        if (ht->Integral() > 0) ht->Scale(1. / ht->Integral());
                        if (hrt->Integral() > 0) hrt->Scale(1. / hrt->Integral());
                        hr->SetLineColor(kRed);
                        hr->SetLineWidth(2);
                        hr->SetFillStyle(0);
                        ht->SetLineColor(kAzure + 2);
                        ht->SetLineWidth(2);
                        ht->SetFillStyle(0);
                        hrt->SetLineColor(kOrange + 1);
                        hrt->SetLineWidth(2);
                        hrt->SetFillStyle(0);
                        double ym = std::max({hr->GetMaximum(), ht->GetMaximum(), hrt->GetMaximum()}) * 1.4;
                        hr->SetMaximum(ym);
                        hr->SetMinimum(0);
                        hr->GetXaxis()->SetTitle("Response");
                        hr->GetYaxis()->SetTitle("Normalised");
                        hr->Draw("HIST");
                        if (ht->GetEntries() > 0) ht->Draw("HIST SAME");
                        if (hrt->GetEntries() > 0) hrt->Draw("HIST SAME");
                        TLegend* l = new TLegend(0.48, 0.65, 0.92, 0.88);
                        l->SetBorderSize(0);
                        l->SetTextSize(0.023);
                        double sR = (s == 0)   ? mc_sigmas[m][bin]
                                    : (s == 1) ? bkg_sigmas[m][bin]
                                               : combined_sigmas[m][bin];
                        l->AddEntry(hr, Form("%s (#sigma=%.3f) %s", RLab(), sR, sLab[s]), "l");
                        l->AddEntry(ht, Form("R_{truth} (#sigma=%.3f)", psf[s].sig_truth), "l");
                        l->AddEntry(hrt, Form("reco/truth (#sigma=%.3f)", psf[s].sig_MC), "l");
                        l->Draw();
                        AtlasUtils::DrawAtlasLabel(0.25, 0.85);
                        c->SaveAs(Form("%s/JER_R_DB_overlay_%s_ptbin%d.pdf", plotBase.c_str(), sTag[s], bin));
                        delete c;
                    }
                }

                // Per-sample console summary
                std::cout << Form("  ┌── JER per sample [%s %s %s] ──\n", massCutLabels[m],
                                  CurrentPtBinLabel(bin), wp_label);
                for (int i = 0; i < 3; i++)
                    std::cout << Form(
                        "  │ %-10s sigma_d=%.4f sigma_MC=%.4f sigma_tr=%.4f sigma^MC=%.4f sigma_in=%.4f sigma^MCi=%.4f Deltacl=%.3f\n",
                        psLabels[i], psf[i].sig_RDB, psf[i].sig_RDB_MC, psf[i].sig_truth, psf[i].sig_MC,
                        psf[i].sig_insitu, psf[i].sig_MCi, psf[i].dC);
                std::cout << "  └──────────────────────────────────────────────────\n";

                // JER closure diagnostic: one table per variant and WP
                // (jer_closure_diagnostic.txt), to diff a variant against nominal
                // and spot bad cells. NEGVAR flags sigma_reco <= sigma_truth, which
                // cannot happen physically (sigma_reco^2 = sigma_truth^2 + sigma_det^2)
                // and leaves sigma^MC_insitu undefined, so R_JER is unstable there.
                {
                    std::string wpDir = Form(((RunConfig::PlotDir() + "/%s/Response_Fit_pTbins_Mcuts_")
                                              + CurrentChannel() + "/%s/")
                                                 .c_str(),
                                             cfg.plot_subdir.c_str(), wp_label);
                    std::string diagPath = wpDir + "jer_closure_diagnostic.txt";
                    const bool firstCell = (m == 0 && bin == 0);
                    std::ofstream diag(diagPath, firstCell ? std::ios::trunc : std::ios::app);
                    if (diag.is_open()) {
                        if (firstCell) {
                            diag << "# JER closure diagnostic  -  channel=" << CurrentChannel()
                                 << "  variant=" << cfg.key << "  WP=" << wp_label << "\n";
                            diag << "# sig_in   = sqrt(sig_data^2  - sig_tr^2)   (in-situ, data)\n";
                            diag << "# sig^MCi  = sqrt(sig_MC^2    - sig_tr^2)   (in-situ, MC reco)\n";
                            diag << "# R_JER    = sig_in / sig^MCi              (JER double-ratio term)\n";
                            diag
                                << "# dClos    = (sig_r/t - sig^MCi) / sig_in   (closure of the two MC routes)\n";
                            diag << "# FLAG: NEGVAR-MC = sig_MC<=sig_tr (sig^MCi undefined, R_JER unstable);"
                                    "  NEGVAR-data = sig_data<=sig_tr;  ok otherwise\n";
                            diag << Form("# %-9s %-22s %-12s | %8s %8s %8s %8s %8s %8s | %8s %8s | %s\n",
                                         "sample", "pTbin", "masscut", "sig_data", "sig_MC", "sig_tr",
                                         "sig_r/t", "sig_in", "sig^MCi", "R_JER", "dClos", "FLAG");
                        }
                        for (int i = 0; i < 3; ++i) {
                            const auto& p = psf[i];
                            const char* flag = (p.sig_truth <= 0)              ? "no-truth"
                                               : (p.sig_RDB_MC <= p.sig_truth) ? "NEGVAR-MC"
                                               : (p.sig_RDB <= p.sig_truth)    ? "NEGVAR-data"
                                                                               : "ok";
                            std::string rjer = (p.sig_MCi > 0)
                                                   ? std::string(Form("%8.4f", p.sig_insitu / p.sig_MCi))
                                                   : std::string("     n/a");
                            diag << Form(
                                "  %-9s %-22s %-12s | %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f | %8s %8.3f | %s\n",
                                psLabels[i], CurrentPtBinLabel(bin), massCutLabels[m], p.sig_RDB,
                                p.sig_RDB_MC, p.sig_truth, p.sig_MC, p.sig_insitu, p.sig_MCi, rjer.c_str(),
                                p.dC, flag);
                        }
                    }
                }
            }

            // Stack plot: data vs signal + background (+ minor backgrounds), combined MC overlaid
            TCanvas* stack_canvas =
                new TCanvas(Form("stack_canvas_masscut%s_bin%d_%s", massCutLabels[m], bin, wp_label),
                            Form("R_response Comparison for %s (Mass Cut: %s, WP: %s)",
                                 CurrentPtBinLabel(bin), massCutLabels[m], wp_label),
                            800, 600);
            stack_canvas->cd();

            // Clone histograms for comparison
            TH1F* h_data_clone = (TH1F*)h_data_R_balance_ptbin[m][bin]->Clone(
                Form("h_data_clone_masscut%s_bin%d_%s", massCutLabels[m], bin, wp_label));
            TH1F* h_mc_clone = (TH1F*)h_mc_R_balance_ptbin[m][bin]->Clone(
                Form("h_mc_clone_masscut%s_bin%d_%s", massCutLabels[m], bin, wp_label));
            TH1F* h_bkg_clone = (TH1F*)h_multijets_R_balance_ptbin[m][bin]->Clone(
                Form("h_bkg_clone_masscut%s_bin%d_%s", massCutLabels[m], bin, wp_label));
            TH1F* h_combined_clone = (TH1F*)h_combined_MC_R_balance_ptbin[m][bin]->Clone(
                Form("h_combined_clone_masscut%s_bin%d_%s", massCutLabels[m], bin, wp_label));

            // Style the clones
            h_data_clone->SetMarkerStyle(20);
            h_data_clone->SetMarkerColor(kBlack);
            h_data_clone->SetLineColor(kBlack);

            h_mc_clone->SetFillColor(kRed);
            h_mc_clone->SetLineColor(kRed);
            h_mc_clone->SetFillStyle(1001);

            h_bkg_clone->SetFillColor(kBlue);
            h_bkg_clone->SetLineColor(kBlue);
            h_bkg_clone->SetFillStyle(1001);

            h_combined_clone->SetFillColor(kGreen + 2);
            h_combined_clone->SetLineColor(kGreen + 2);
            h_combined_clone->SetFillStyle(1001);

            // Create stack - extra samples at bottom, multijets above, signal on top
            THStack* stack =
                new THStack(Form("stack_masscut%s_bin%d_%s", massCutLabels[m], bin, wp_label),
                            Form("%s Distribution for %s;%s;Events", RLab(), CurrentPtBinLabel(bin), RLab()));
            // Add extra samples first (bottom of stack)
            for (int idx = 0; idx < nSC; ++idx) {
                if (idx == SampleConfig::kData || idx == SampleConfig::kMultijets
                    || idx == SampleConfig::kZbb)
                    continue;
                if (extraLoaded[idx]) {
                    TH1F* h_extra_clone = (TH1F*)h_extra_R_balance_ptbin[idx][m][bin]->Clone(
                        Form("h_extra%d_clone_masscut%s_bin%d_%s", idx, massCutLabels[m], bin, wp_label));
                    h_extra_clone->SetFillStyle(1001);
                    stack->Add(h_extra_clone);
                }
            }
            stack->Add(h_bkg_clone);
            stack->Add(h_mc_clone);

            stack->Draw("HIST");
            h_combined_clone->Draw("HIST SAME");
            h_data_clone->Draw("E SAME");

            // Legend
            TLegend* stack_legend = new TLegend(0.7, 0.7, 0.95, 0.95);
            stack_legend->SetBorderSize(0);
            stack_legend->SetFillStyle(0);
            stack_legend->SetTextSize(0.025);
            stack_legend->AddEntry(h_data_clone, "Data", "lp");
            stack_legend->AddEntry(h_mc_clone, SigLab(), "f");
            stack_legend->AddEntry(h_bkg_clone, BkgLab(), "f");
            for (int idx = 0; idx < nSC; ++idx) {
                if (idx == SampleConfig::kData || idx == SampleConfig::kMultijets
                    || idx == SampleConfig::kZbb)
                    continue;
                if (extraLoaded[idx])
                    stack_legend->AddEntry(h_extra_R_balance_ptbin[idx][m][bin], scDisplayNames[idx].c_str(),
                                           "f");
            }
            stack_legend->AddEntry(h_combined_clone, "Combined MC (Signal+Bkg)", "l");
            stack_legend->Draw();

            AtlasUtils::DrawAtlasLabel(0.2, 0.85);
            TLatex text_stack;
            text_stack.SetTextSize(0.035);
            text_stack.DrawLatexNDC(0.2, 0.76, CurrentPtBinLabel(bin));

            // Add mass cut information
            text_stack.DrawLatexNDC(0.2, 0.71, Form("Mass: %s", massCutLabels[m]));

            stack_canvas->SaveAs(Form(((RunConfig::PlotDir() + "/%s/Response_Fit_pTbins_Mcuts_")
                                       + CurrentChannel() + "/%s/%s/R_DB_stack_comparison_ptbin%d.pdf")
                                          .c_str(),
                                      cfg.plot_subdir.c_str(), wp_label, fileSystemMassCutLabels[m], bin));

            delete stack_canvas;
            delete h_data_clone;
            delete h_mc_clone;
            delete h_bkg_clone;
            delete h_combined_clone;
            delete stack;
            delete stack_legend;

            // Create canvas with ratio plot
            TCanvas* combinedstack_canvas =
                new TCanvas(Form("combinedstack_canvas_masscut%s_bin%d_%s", massCutLabels[m], bin, wp_label),
                            Form("R_response Comparison for %s (Mass Cut: %s, WP: %s)",
                                 CurrentPtBinLabel(bin), massCutLabels[m], wp_label),
                            600, 600);
            combinedstack_canvas->cd();

            // Configure pads
            TPad* pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
            pad1->SetBottomMargin(0.02); // Reduced bottom margin
            pad1->Draw();
            pad1->cd(); // Switch to pad1

            // Clone histograms
            TH1F* h_data_clone_comb = (TH1F*)h_data_R_balance_ptbin[m][bin]->Clone(
                Form("h_data_clone_comb_masscut%s_bin%d_%s", massCutLabels[m], bin, wp_label));
            TH1F* h_mc_clone_comb = (TH1F*)h_mc_R_balance_ptbin[m][bin]->Clone(
                Form("h_mc_clone_comb_masscut%s_bin%d_%s", massCutLabels[m], bin, wp_label));
            TH1F* h_bkg_clone_comb = (TH1F*)h_multijets_R_balance_ptbin[m][bin]->Clone(
                Form("h_bkg_clone_comb_masscut%s_bin%d_%s", massCutLabels[m], bin, wp_label));
            TH1F* h_combined_clone_comb = (TH1F*)h_combined_MC_R_balance_ptbin[m][bin]->Clone(
                Form("h_combined_clone_comb_masscut%s_bin%d_%s", massCutLabels[m], bin, wp_label));

            // Style clones
            h_data_clone_comb->SetMarkerStyle(20);
            h_data_clone_comb->SetMarkerColor(kBlack);
            h_data_clone_comb->SetLineColor(kBlack);

            h_mc_clone_comb->SetFillColor(kRed);
            h_mc_clone_comb->SetLineColor(kRed);
            h_mc_clone_comb->SetFillStyle(1001);

            h_bkg_clone_comb->SetFillColor(kBlue);
            h_bkg_clone_comb->SetLineColor(kBlue);
            h_bkg_clone_comb->SetFillStyle(1001);

            h_combined_clone_comb->SetFillColor(kGreen + 2);
            h_combined_clone_comb->SetLineColor(kGreen + 2);
            h_combined_clone_comb->SetFillStyle(1001);

            // Create and draw stack
            THStack* stack_comb =
                new THStack(Form("stack_masscut%s_bin%d_%s", massCutLabels[m], bin, wp_label),
                            Form("%s Distribution for %s;;Events", RLab(),
                                 CurrentPtBinLabel(bin))); // Note empty x-axis title
            stack_comb->Add(h_bkg_clone_comb);
            stack_comb->Add(h_mc_clone_comb);
            stack_comb->Draw("HIST");
            h_data_clone_comb->Draw("EP SAME");

            // Data and combined MC fits, each drawn over its fit window with the
            // fitted model (Gaussian by default, or the EASYBJETS_FIT_MODEL shape).
            TF1* gausFit_data_recreate_inner = FitUtils::MakeModelTF1(
                fit_data,
                Form("gausFit_data_recreate_inner_masscut%s_bin%d_%s", massCutLabels[m], bin, wp_label));
            gausFit_data_recreate_inner->SetLineColor(kViolet - 4);
            gausFit_data_recreate_inner->SetLineStyle(1);
            gausFit_data_recreate_inner->SetLineWidth(3);
            gausFit_data_recreate_inner->Draw("same");

            TF1* gausFit_combined_recreate_inner = FitUtils::MakeModelTF1(
                fit_combined,
                Form("gausFit_combined_recreate_inner_masscut%s_bin%d_%s", massCutLabels[m], bin, wp_label));
            gausFit_combined_recreate_inner->SetLineColor(kCyan - 7);
            gausFit_combined_recreate_inner->SetLineStyle(1);
            gausFit_combined_recreate_inner->SetLineWidth(3);
            gausFit_combined_recreate_inner->Draw("same");

            // Hide x-axis labels on top pad
            stack_comb->GetXaxis()->SetLabelSize(0);
            stack_comb->GetXaxis()->SetTitleSize(0);

            // Style y-axis
            stack_comb->GetYaxis()->SetTitleSize(0.05);
            stack_comb->GetYaxis()->SetTitle("Events");
            stack_comb->GetYaxis()->SetLabelSize(0.04);

            stack_comb->SetMaximum(stack_comb->GetMaximum() * 1.5); 

            TH1F* h_mc_total = (TH1F*)stack_comb->GetStack()->Last(); 

            TGraphAsymmErrors* err_band = HistUtils::MakeUncertaintyBand(h_mc_total);
            err_band->Draw("2 SAME");

            // Legend
            TLegend* stack_legend_comb = new TLegend(0.6, 0.6, 0.95, 0.92);
            stack_legend_comb->SetBorderSize(0);  // No border
            stack_legend_comb->SetFillStyle(0);   // Transparent background
            stack_legend_comb->SetTextSize(0.03); // Readable text size
            stack_legend_comb->AddEntry(h_data_clone_comb, "Data", "lep");
            stack_legend_comb->AddEntry(h_mc_clone_comb, SigLab(), "f");
            stack_legend_comb->AddEntry(h_bkg_clone_comb, BkgLab(), "f");
            stack_legend_comb->AddEntry(gausFit_data_recreate_inner, "Data Fit", "l");
            stack_legend_comb->AddEntry(gausFit_combined_recreate_inner, "Combined MC (Sig + Bkg) Fit", "l");
            stack_legend_comb->AddEntry(err_band, "MC Stat. Unc.", "f");
            stack_legend_comb->Draw();

            // Add ATLAS label to the plot
            AtlasUtils::DrawAtlasLabel(0.2, 0.85);
            TLatex text_stack_comb;
            text_stack_comb.SetTextSize(0.035);
            text_stack_comb.DrawLatexNDC(0.2, 0.76, CurrentPtBinLabel(bin));

            // Add mass cut information
            text_stack_comb.DrawLatexNDC(0.2, 0.71, Form("Mass: %s", massCutLabels[m]));

            pad1->Modified();

            // Switch to ratio pad
            combinedstack_canvas->cd();
            TPad* pad2 = new TPad("pad2", "pad2", 0, 0.0, 1, 0.3);
            pad2->SetTopMargin(0.05);
            pad2->SetBottomMargin(0.25);
            pad2->Draw();
            pad2->cd();

            // Create ratio histogram
            TH1F* h_ratio = (TH1F*)h_data_clone_comb->Clone(
                Form("h_ratio_masscut%s_bin%d_%s", massCutLabels[m], bin, wp_label));
            h_ratio->Divide(h_combined_clone_comb);

            // Style ratio plot
            h_ratio->SetTitle("");
            h_ratio->SetMarkerStyle(21);
            h_ratio->SetMarkerSize(1.0);
            h_ratio->SetLineColor(kBlack);
            h_ratio->GetYaxis()->SetRangeUser(0.5, 1.5);
            h_ratio->GetXaxis()->SetTitle(RLab());
            h_ratio->GetYaxis()->SetTitle("Data/MC");
            h_ratio->GetYaxis()->SetTitleSize(0.12);
            h_ratio->GetYaxis()->SetTitleOffset(0.4);
            h_ratio->GetYaxis()->SetLabelSize(0.10);
            h_ratio->GetXaxis()->SetTitleSize(0.12);
            h_ratio->GetXaxis()->SetTitleOffset(0.75);
            h_ratio->GetXaxis()->SetLabelSize(0.10);
            h_ratio->GetYaxis()->SetNdivisions(505); // Fewer divisions for cleaner look
            h_ratio->GetYaxis()->CenterTitle();

            // Draw the ratio plot
            h_ratio->Draw("El same");

            // MC statistical uncertainty band
            TGraphAsymmErrors* err_band_ratio = HistUtils::MakeRatioUncertaintyBand(h_combined_clone_comb);
            err_band_ratio->Draw("2 SAME");

            // Add a horizontal reference line at y=1
            double xMin = h_ratio->GetXaxis()->GetXmin();
            double xMax = h_ratio->GetXaxis()->GetXmax();
            TLine* line = new TLine(xMin, 1.0, xMax, 1.0);
            line->SetLineColor(kBlack);
            line->SetLineWidth(2);
            line->SetLineStyle(2); // Dashed line
            line->Draw("SAME");

            // Update both pads to ensure proper display
            pad2->Update();

            // Save canvas
            combinedstack_canvas->Update();
            combinedstack_canvas->SaveAs(Form(((RunConfig::PlotDir() + "/%s/Response_Fit_pTbins_Mcuts_")
                                               + CurrentChannel() + "/%s/%s/R_DB_stack_ratio_ptbin%d.pdf")
                                                  .c_str(),
                                              cfg.plot_subdir.c_str(), wp_label, fileSystemMassCutLabels[m],
                                              bin));

            // Cleanup for combinedstack_canvas
            delete pad1;
            delete pad2;
            delete combinedstack_canvas;
            delete h_data_clone_comb;
            delete h_mc_clone_comb;
            delete h_bkg_clone_comb;
            delete h_combined_clone_comb;
            delete stack_comb;
            delete gausFit_data_recreate_inner;
            delete gausFit_combined_recreate_inner;
            delete err_band;
            delete stack_legend_comb;
            delete h_ratio;
            delete err_band_ratio;
            delete line;
        }

        // Close summary file
        if (summaryFile.is_open()) {
            summaryFile.close();
            std::cout << "Summary of all fit parameters saved to " << summaryFileName.Data() << std::endl;
        }

        // Create comparison plots of mean and sigma vs pT for this mass cut
        // Mean comparison
        TCanvas* mean_canvas = new TCanvas(
            Form("mean_canvas_masscut%s_%s", massCutLabels[m], wp_label),
            Form("R_response Mean vs pT (Mass Cut: %s, WP: %s)", massCutLabels[m], wp_label), 800, 600);
        mean_canvas->cd();

        // Create graph for each sample
        TGraphErrors* g_mc_mean =
            new TGraphErrors(CurrentNPtBins(), ptBinCenters, mc_means[m], nullptr, mc_mean_errs[m]);
        TGraphErrors* g_data_mean =
            new TGraphErrors(CurrentNPtBins(), ptBinCenters, data_means[m], nullptr, data_mean_errs[m]);
        TGraphErrors* g_bkg_mean =
            new TGraphErrors(CurrentNPtBins(), ptBinCenters, bkg_means[m], nullptr, bkg_mean_errs[m]);
        TGraphErrors* g_combined_mean = new TGraphErrors(CurrentNPtBins(), ptBinCenters, combined_means[m],
                                                         nullptr, combined_mean_errs[m]);

        // Style graphs
        g_mc_mean->SetMarkerStyle(20);
        g_mc_mean->SetMarkerColor(kRed);
        g_mc_mean->SetLineColor(kRed);
        g_mc_mean->SetMarkerSize(1.2);

        g_data_mean->SetMarkerStyle(21);
        g_data_mean->SetMarkerColor(kBlack);
        g_data_mean->SetLineColor(kBlack);
        g_data_mean->SetMarkerSize(1.2);

        g_bkg_mean->SetMarkerStyle(22);
        g_bkg_mean->SetMarkerColor(kBlue);
        g_bkg_mean->SetLineColor(kBlue);
        g_bkg_mean->SetMarkerSize(1.2);

        g_combined_mean->SetMarkerStyle(23);
        g_combined_mean->SetMarkerColor(kGreen + 2);
        g_combined_mean->SetLineColor(kGreen + 2);
        g_combined_mean->SetMarkerSize(1.2);

        // Draw graphs
        g_mc_mean->SetTitle(Form("%s Mean vs p_{T};p_{T} [GeV];#LT %s#GT", RLab(), RLab()));
        g_mc_mean->GetYaxis()->SetRangeUser(0.95, 1.05); // Adjust range as needed
        g_mc_mean->Draw("APL");
        g_data_mean->Draw("PL SAME");
        g_bkg_mean->Draw("PL SAME");
        g_combined_mean->Draw("PL SAME");

        // Add legend
        TLegend* mean_legend = new TLegend(0.6, 0.7, 0.89, 0.89);
        mean_legend->SetBorderSize(0);
        mean_legend->AddEntry(g_mc_mean, SigLab(), "pl");
        mean_legend->AddEntry(g_data_mean, "Data", "pl");
        mean_legend->AddEntry(g_bkg_mean, BkgLab(), "pl");
        mean_legend->AddEntry(g_combined_mean, "Combined MC (Signal+Bkg)", "pl");
        mean_legend->Draw();

        AtlasUtils::DrawAtlasLabel(0.18, 0.85);

        // Add mass cut information
        TLatex mean_text;
        mean_text.SetTextSize(0.035);
        mean_text.DrawLatexNDC(0.2, 0.76, Form("Mass: %s", massCutLabels[m]));

        mean_canvas->SaveAs(Form(((RunConfig::PlotDir() + "/%s/Response_Fit_pTbins_Mcuts_") + CurrentChannel()
                                  + "/%s/%s/R_DB_mean_vs_pT.pdf")
                                     .c_str(),
                                 cfg.plot_subdir.c_str(), wp_label, fileSystemMassCutLabels[m]));

        delete mean_canvas;
        delete g_mc_mean;
        delete g_data_mean;
        delete g_bkg_mean;
        delete g_combined_mean;
        delete mean_legend;

        // Sigma comparison
        TCanvas* sigma_canvas = new TCanvas(
            Form("sigma_canvas_masscut%s_%s", massCutLabels[m], wp_label),
            Form("R_response Sigma vs pT (Mass Cut: %s, WP: %s)", massCutLabels[m], wp_label), 800, 600);
        sigma_canvas->cd();

        // Create graphs
        TGraphErrors* g_mc_sigma =
            new TGraphErrors(CurrentNPtBins(), ptBinCenters, mc_sigmas[m], nullptr, mc_sigma_errs[m]);
        TGraphErrors* g_data_sigma =
            new TGraphErrors(CurrentNPtBins(), ptBinCenters, data_sigmas[m], nullptr, data_sigma_errs[m]);
        TGraphErrors* g_bkg_sigma =
            new TGraphErrors(CurrentNPtBins(), ptBinCenters, bkg_sigmas[m], nullptr, bkg_sigma_errs[m]);
        TGraphErrors* g_combined_sigma = new TGraphErrors(CurrentNPtBins(), ptBinCenters, combined_sigmas[m],
                                                          nullptr, combined_sigma_errs[m]);

        // Style graphs
        g_mc_sigma->SetMarkerStyle(20);
        g_mc_sigma->SetMarkerColor(kRed);
        g_mc_sigma->SetLineColor(kRed);
        g_mc_sigma->SetMarkerSize(1.2);

        g_data_sigma->SetMarkerStyle(21);
        g_data_sigma->SetMarkerColor(kBlack);
        g_data_sigma->SetLineColor(kBlack);
        g_data_sigma->SetMarkerSize(1.2);

        g_bkg_sigma->SetMarkerStyle(22);
        g_bkg_sigma->SetMarkerColor(kBlue);
        g_bkg_sigma->SetLineColor(kBlue);
        g_bkg_sigma->SetMarkerSize(1.2);

        g_combined_sigma->SetMarkerStyle(23);
        g_combined_sigma->SetMarkerColor(kGreen + 2);
        g_combined_sigma->SetLineColor(kGreen + 2);
        g_combined_sigma->SetMarkerSize(1.2);

        // Draw graphs
        g_mc_sigma->SetTitle(Form("%s Width vs p_{T};p_{T} [GeV];#sigma (%s)", RLab(), RLab()));
        g_mc_sigma->GetYaxis()->SetRangeUser(0.0, 0.2); // Adjust range as needed
        g_mc_sigma->Draw("APL");
        g_data_sigma->Draw("PL SAME");
        g_bkg_sigma->Draw("PL SAME");
        g_combined_sigma->Draw("PL SAME");

        // Add legend
        TLegend* sigma_legend = new TLegend(0.6, 0.7, 0.89, 0.89);
        sigma_legend->SetBorderSize(0);
        sigma_legend->AddEntry(g_mc_sigma, SigLab(), "pl");
        sigma_legend->AddEntry(g_data_sigma, "Data", "pl");
        sigma_legend->AddEntry(g_bkg_sigma, BkgLab(), "pl");
        sigma_legend->AddEntry(g_combined_sigma, "Combined MC (Signal+Bkg)", "pl");
        sigma_legend->Draw();

        AtlasUtils::DrawAtlasLabel(0.18, 0.85);

        // Add mass cut information
        TLatex sigma_text;
        sigma_text.SetTextSize(0.035);
        sigma_text.DrawLatexNDC(0.2, 0.76, Form("Mass: %s", massCutLabels[m]));

        sigma_canvas->SaveAs(Form(((RunConfig::PlotDir() + "/%s/Response_Fit_pTbins_Mcuts_")
                                   + CurrentChannel() + "/%s/%s/R_DB_sigma_vs_pT.pdf")
                                      .c_str(),
                                  cfg.plot_subdir.c_str(), wp_label, fileSystemMassCutLabels[m]));

        delete sigma_canvas;
        delete g_mc_sigma;
        delete g_data_sigma;
        delete g_bkg_sigma;
        delete g_combined_sigma;
        delete sigma_legend;

        // Create ratio plots (data/MC and background/MC)
        TCanvas* ratio_canvas = new TCanvas(
            Form("ratio_canvas_masscut%s_%s", massCutLabels[m], wp_label),
            Form("R_response Mean Ratio vs pT (Mass Cut: %s, WP: %s)", massCutLabels[m], wp_label), 800, 600);
        ratio_canvas->cd();

        // Calculate ratios
        double data_mc_ratio[MAX_PT_BINS], data_mc_ratio_err[MAX_PT_BINS];
        double bkg_mc_ratio[MAX_PT_BINS], bkg_mc_ratio_err[MAX_PT_BINS];
        double data_combined_ratio[MAX_PT_BINS], data_combined_ratio_err[MAX_PT_BINS];

        for (int i = 0; i < CurrentNPtBins(); i++) {
            // Data/MC ratio
            if (mc_means[m][i] != 0) {
                data_mc_ratio[i] = data_means[m][i] / mc_means[m][i];
                // Error propagation for ratio
                data_mc_ratio_err[i] = data_mc_ratio[i]
                                       * std::sqrt(std::pow(data_mean_errs[m][i] / data_means[m][i], 2)
                                                   + std::pow(mc_mean_errs[m][i] / mc_means[m][i], 2));
            } else {
                data_mc_ratio[i] = 0.0;
                data_mc_ratio_err[i] = 0.0;
            }

            // Background/MC ratio
            if (mc_means[m][i] != 0) {
                bkg_mc_ratio[i] = bkg_means[m][i] / mc_means[m][i];
                // Error propagation for ratio
                bkg_mc_ratio_err[i] = bkg_mc_ratio[i]
                                      * std::sqrt(std::pow(bkg_mean_errs[m][i] / bkg_means[m][i], 2)
                                                  + std::pow(mc_mean_errs[m][i] / mc_means[m][i], 2));
            } else {
                bkg_mc_ratio[i] = 0.0;
                bkg_mc_ratio_err[i] = 0.0;
            }

            // Data/Combined MC ratio
            if (combined_means[m][i] != 0) {
                data_combined_ratio[i] = data_means[m][i] / combined_means[m][i];
                // Error propagation for ratio
                data_combined_ratio_err[i] =
                    data_combined_ratio[i]
                    * std::sqrt(std::pow(data_mean_errs[m][i] / data_means[m][i], 2)
                                + std::pow(combined_mean_errs[m][i] / combined_means[m][i], 2));
            } else {
                data_combined_ratio[i] = 0.0;
                data_combined_ratio_err[i] = 0.0;
            }
        }

        // Create graphs
        TGraphErrors* g_data_mc_ratio =
            new TGraphErrors(CurrentNPtBins(), ptBinCenters, data_mc_ratio, nullptr, data_mc_ratio_err);
        TGraphErrors* g_bkg_mc_ratio =
            new TGraphErrors(CurrentNPtBins(), ptBinCenters, bkg_mc_ratio, nullptr, bkg_mc_ratio_err);
        TGraphErrors* g_data_combined_ratio = new TGraphErrors(
            CurrentNPtBins(), ptBinCenters, data_combined_ratio, nullptr, data_combined_ratio_err);

        // Style
        g_data_mc_ratio->SetMarkerStyle(21);
        g_data_mc_ratio->SetMarkerColor(kBlack);
        g_data_mc_ratio->SetLineColor(kBlack);
        g_data_mc_ratio->SetMarkerSize(1.2);

        g_bkg_mc_ratio->SetMarkerStyle(22);
        g_bkg_mc_ratio->SetMarkerColor(kBlue);
        g_bkg_mc_ratio->SetLineColor(kBlue);
        g_bkg_mc_ratio->SetMarkerSize(1.2);

        g_data_combined_ratio->SetMarkerStyle(23);
        g_data_combined_ratio->SetMarkerColor(kGreen + 2);
        g_data_combined_ratio->SetLineColor(kGreen + 2);
        g_data_combined_ratio->SetMarkerSize(1.2);

        // Draw
        g_data_mc_ratio->SetTitle(Form("Ratio of %s Mean to MC;p_{T} [GeV];Ratio to MC", RLab()));
        g_data_mc_ratio->GetYaxis()->SetRangeUser(0.95, 1.05); // Adjust range
        g_data_mc_ratio->Draw("APL");
        g_bkg_mc_ratio->Draw("PL SAME");
        g_data_combined_ratio->Draw("PL SAME");

        // First make sure we get the actual axis range
        double xmin = g_data_mc_ratio->GetXaxis()->GetXmin();
        double xmax = g_data_mc_ratio->GetXaxis()->GetXmax();

        // Add a reference line at 1
        TLine* line_ratio = new TLine(xmin, 1.0, xmax, 1.0);
        line_ratio->SetLineStyle(2);
        line_ratio->SetLineColor(kRed);
        line_ratio->Draw();

        // Legend
        TLegend* ratio_legend = new TLegend(0.6, 0.7, 0.89, 0.89);
        ratio_legend->SetBorderSize(0);
        ratio_legend->AddEntry(g_data_mc_ratio, "Data/Signal MC", "pl");
        ratio_legend->AddEntry(g_bkg_mc_ratio, "Multijets/Signal MC", "pl");
        ratio_legend->AddEntry(g_data_combined_ratio, "Data/Combined MC", "pl");
        ratio_legend->AddEntry(line_ratio, "Ratio = 1", "l");
        ratio_legend->Draw();

        AtlasUtils::DrawAtlasLabel(0.18, 0.85);

        // Add mass cut information
        TLatex ratio_text;
        ratio_text.SetTextSize(0.035);
        ratio_text.DrawLatexNDC(0.2, 0.8, Form("Mass: %s", massCutLabels[m]));

        ratio_canvas->SaveAs(Form(((RunConfig::PlotDir() + "/%s/Response_Fit_pTbins_Mcuts_")
                                   + CurrentChannel() + "/%s/%s/R_DB_mean_ratio_vs_pT.pdf")
                                      .c_str(),
                                  cfg.plot_subdir.c_str(), wp_label, fileSystemMassCutLabels[m]));

        delete ratio_canvas;
        delete g_data_mc_ratio;
        delete g_bkg_mc_ratio;
        delete g_data_combined_ratio;
        delete line_ratio;
        delete ratio_legend;

        // Create data/combined MC comparison plot
        TCanvas* data_combined_canvas = new TCanvas(
            Form("data_combined_canvas_masscut%s_%s", massCutLabels[m], wp_label),
            Form("Data vs Combined MC (Mass Cut: %s, WP: %s)", massCutLabels[m], wp_label), 800, 600);
        data_combined_canvas->cd();

        // Draw only data/combined MC ratio
        TGraphErrors* g_data_combined_ratio_only = new TGraphErrors(
            CurrentNPtBins(), ptBinCenters, data_combined_ratio, nullptr, data_combined_ratio_err);
        g_data_combined_ratio_only->SetTitle("Ratio of Data to Combined MC;p_{T} [GeV];Data/Combined MC");
        g_data_combined_ratio_only->GetYaxis()->SetRangeUser(0.95, 1.05); // Adjust range
        g_data_combined_ratio_only->SetMarkerStyle(23);
        g_data_combined_ratio_only->SetMarkerColor(kGreen + 2);
        g_data_combined_ratio_only->SetLineColor(kGreen + 2);
        g_data_combined_ratio_only->SetMarkerSize(1.2);
        g_data_combined_ratio_only->Draw("APL");

        // Add a reference line at 1
        TLine* line2 =
            new TLine(CurrentPtBinsDisplay()[0], 1.0, CurrentPtBinsDisplay()[CurrentNPtBins()], 1.0);
        line2->SetLineStyle(2);
        line2->SetLineColor(kRed);
        line2->Draw();

        // Legend
        TLegend* data_combined_legend = new TLegend(0.6, 0.7, 0.89, 0.89);
        data_combined_legend->SetBorderSize(0);
        data_combined_legend->AddEntry(g_data_combined_ratio_only, "Data/Combined MC", "pl");
        data_combined_legend->AddEntry(line2, "MC Reference", "l");
        data_combined_legend->Draw();

        AtlasUtils::DrawAtlasLabel(0.18, 0.85);

        // Add mass cut information
        TLatex data_comb_text;
        data_comb_text.SetTextSize(0.035);
        data_comb_text.DrawLatexNDC(0.2, 0.8, Form("Mass: %s", massCutLabels[m]));

        data_combined_canvas->SaveAs(Form(((RunConfig::PlotDir() + "/%s/Response_Fit_pTbins_Mcuts_")
                                           + CurrentChannel() + "/%s/%s/R_DB_data_combined_ratio_vs_pT.pdf")
                                              .c_str(),
                                          cfg.plot_subdir.c_str(), wp_label, fileSystemMassCutLabels[m]));

        delete data_combined_canvas;
        delete g_data_combined_ratio_only;
        delete line2;
        delete data_combined_legend;
    }

    // Now create comparison plots across different mass cuts
    // For each pT bin, compare results across mass cuts
    for (int bin = 0; bin < CurrentNPtBins(); bin++) {
        // Create comparison of means across mass cuts for this pT bin
        TCanvas* mass_mean_canvas = new TCanvas(
            Form("mass_mean_canvas_ptbin%d_%s", bin, wp_label),
            Form("R_response Mean vs Mass Cut (pT bin: %s, WP: %s)", CurrentPtBinLabel(bin), wp_label), 800,
            600);
        mass_mean_canvas->cd();

        double massCutX[MAX_MASS_CUTS];
        for (int m = 0; m < nMassCuts; m++) {
            massCutX[m] = m; // Use indices as x-coordinates
        }

        // Extract means for this pT bin across different mass cuts
        double mc_mean_values[MAX_MASS_CUTS], mc_mean_errors[MAX_MASS_CUTS];
        double data_mean_values[MAX_MASS_CUTS], data_mean_errors[MAX_MASS_CUTS];
        double bkg_mean_values[MAX_MASS_CUTS], bkg_mean_errors[MAX_MASS_CUTS];
        double combined_mean_values[MAX_MASS_CUTS], combined_mean_errors[MAX_MASS_CUTS];

        for (int m = 0; m < nMassCuts; m++) {
            mc_mean_values[m] = mc_means[m][bin];
            mc_mean_errors[m] = mc_mean_errs[m][bin];
            data_mean_values[m] = data_means[m][bin];
            data_mean_errors[m] = data_mean_errs[m][bin];
            bkg_mean_values[m] = bkg_means[m][bin];
            bkg_mean_errors[m] = bkg_mean_errs[m][bin];
            combined_mean_values[m] = combined_means[m][bin];
            combined_mean_errors[m] = combined_mean_errs[m][bin];
        }

        // Create graphs
        TGraphErrors* g_mc_mean_mass =
            new TGraphErrors(nMassCuts, massCutX, mc_mean_values, nullptr, mc_mean_errors);
        TGraphErrors* g_data_mean_mass =
            new TGraphErrors(nMassCuts, massCutX, data_mean_values, nullptr, data_mean_errors);
        TGraphErrors* g_bkg_mean_mass =
            new TGraphErrors(nMassCuts, massCutX, bkg_mean_values, nullptr, bkg_mean_errors);
        TGraphErrors* g_combined_mean_mass =
            new TGraphErrors(nMassCuts, massCutX, combined_mean_values, nullptr, combined_mean_errors);

        // Style graphs
        g_mc_mean_mass->SetMarkerStyle(20);
        g_mc_mean_mass->SetMarkerColor(kRed);
        g_mc_mean_mass->SetLineColor(kRed);
        g_mc_mean_mass->SetMarkerSize(1.2);

        g_data_mean_mass->SetMarkerStyle(21);
        g_data_mean_mass->SetMarkerColor(kBlack);
        g_data_mean_mass->SetLineColor(kBlack);
        g_data_mean_mass->SetMarkerSize(1.2);

        g_bkg_mean_mass->SetMarkerStyle(22);
        g_bkg_mean_mass->SetMarkerColor(kBlue);
        g_bkg_mean_mass->SetLineColor(kBlue);
        g_bkg_mean_mass->SetMarkerSize(1.2);

        g_combined_mean_mass->SetMarkerStyle(23);
        g_combined_mean_mass->SetMarkerColor(kGreen + 2);
        g_combined_mean_mass->SetLineColor(kGreen + 2);
        g_combined_mean_mass->SetMarkerSize(1.2);

        // Draw graphs
        g_mc_mean_mass->SetTitle(
            Form("%s Mean vs Mass Cut for %s;Mass Cut;#LT %s#GT", RLab(), CurrentPtBinLabel(bin), RLab()));
        g_mc_mean_mass->GetYaxis()->SetRangeUser(0.95, 1.05); // Adjust range as needed
        g_mc_mean_mass->GetXaxis()->SetLimits(-0.5, nMassCuts - 0.5);
        g_mc_mean_mass->Draw("APL");
        g_data_mean_mass->Draw("PL SAME");
        g_bkg_mean_mass->Draw("PL SAME");
        g_combined_mean_mass->Draw("PL SAME");

        // Add custom x-axis labels
        TAxis* xaxis = g_mc_mean_mass->GetXaxis();
        for (int m = 0; m < nMassCuts; m++) {
            xaxis->SetBinLabel(xaxis->FindBin(m), massCutLabels[m]);
        }
        xaxis->SetLabelSize(0.05);

        // Legend
        TLegend* mass_mean_legend = new TLegend(0.6, 0.75, 0.89, 0.89);
        mass_mean_legend->SetBorderSize(0);
        mass_mean_legend->AddEntry(g_mc_mean_mass, SigLab(), "pl");
        mass_mean_legend->AddEntry(g_data_mean_mass, "Data", "pl");
        mass_mean_legend->AddEntry(g_bkg_mean_mass, BkgLab(), "pl");
        mass_mean_legend->AddEntry(g_combined_mean_mass, "Combined MC (Signal+Bkg)", "pl");
        mass_mean_legend->Draw();

        AtlasUtils::DrawAtlasLabel(0.18, 0.85);
        TLatex mass_mean_text;
        mass_mean_text.SetTextSize(0.035);
        mass_mean_text.DrawLatexNDC(0.2, 0.76, CurrentPtBinLabel(bin));

        mass_mean_canvas->SetGridx();
        mass_mean_canvas->SaveAs(Form(((RunConfig::PlotDir() + "/%s/Response_Fit_pTbins_Mcuts_")
                                       + CurrentChannel() + "/%s/R_DB_mean_vs_masscut_ptbin%d.pdf")
                                          .c_str(),
                                      cfg.plot_subdir.c_str(), wp_label, bin));

        delete mass_mean_canvas;
        delete g_mc_mean_mass;
        delete g_data_mean_mass;
        delete g_bkg_mean_mass;
        delete g_combined_mean_mass;
        delete mass_mean_legend;

        // Also create comparison of sigmas across mass cuts for this pT bin
        TCanvas* mass_sigma_canvas = new TCanvas(
            Form("mass_sigma_canvas_ptbin%d_%s", bin, wp_label),
            Form("R_response Sigma vs Mass Cut (pT bin: %s, WP: %s)", CurrentPtBinLabel(bin), wp_label), 800,
            600);
        mass_sigma_canvas->cd();

        // Extract sigmas for this pT bin across different mass cuts
        double mc_sigma_values[MAX_MASS_CUTS], mc_sigma_errors[MAX_MASS_CUTS];
        double data_sigma_values[MAX_MASS_CUTS], data_sigma_errors[MAX_MASS_CUTS];
        double bkg_sigma_values[MAX_MASS_CUTS], bkg_sigma_errors[MAX_MASS_CUTS];
        double combined_sigma_values[MAX_MASS_CUTS], combined_sigma_errors[MAX_MASS_CUTS];

        for (int m = 0; m < nMassCuts; m++) {
            mc_sigma_values[m] = mc_sigmas[m][bin];
            mc_sigma_errors[m] = mc_sigma_errs[m][bin];
            data_sigma_values[m] = data_sigmas[m][bin];
            data_sigma_errors[m] = data_sigma_errs[m][bin];
            bkg_sigma_values[m] = bkg_sigmas[m][bin];
            bkg_sigma_errors[m] = bkg_sigma_errs[m][bin];
            combined_sigma_values[m] = combined_sigmas[m][bin];
            combined_sigma_errors[m] = combined_sigma_errs[m][bin];
        }

        // Create graphs
        TGraphErrors* g_mc_sigma_mass =
            new TGraphErrors(nMassCuts, massCutX, mc_sigma_values, nullptr, mc_sigma_errors);
        TGraphErrors* g_data_sigma_mass =
            new TGraphErrors(nMassCuts, massCutX, data_sigma_values, nullptr, data_sigma_errors);
        TGraphErrors* g_bkg_sigma_mass =
            new TGraphErrors(nMassCuts, massCutX, bkg_sigma_values, nullptr, bkg_sigma_errors);
        TGraphErrors* g_combined_sigma_mass =
            new TGraphErrors(nMassCuts, massCutX, combined_sigma_values, nullptr, combined_sigma_errors);

        // Style graphs
        g_mc_sigma_mass->SetMarkerStyle(20);
        g_mc_sigma_mass->SetMarkerColor(kRed);
        g_mc_sigma_mass->SetLineColor(kRed);
        g_mc_sigma_mass->SetMarkerSize(1.2);

        g_data_sigma_mass->SetMarkerStyle(21);
        g_data_sigma_mass->SetMarkerColor(kBlack);
        g_data_sigma_mass->SetLineColor(kBlack);
        g_data_sigma_mass->SetMarkerSize(1.2);

        g_bkg_sigma_mass->SetMarkerStyle(22);
        g_bkg_sigma_mass->SetMarkerColor(kBlue);
        g_bkg_sigma_mass->SetLineColor(kBlue);
        g_bkg_sigma_mass->SetMarkerSize(1.2);

        g_combined_sigma_mass->SetMarkerStyle(23);
        g_combined_sigma_mass->SetMarkerColor(kGreen + 2);
        g_combined_sigma_mass->SetLineColor(kGreen + 2);
        g_combined_sigma_mass->SetMarkerSize(1.2);

        // Draw graphs
        g_mc_sigma_mass->SetTitle(
            Form("%s Width vs Mass Cut for %s;Mass Cut;#sigma (%s)", RLab(), CurrentPtBinLabel(bin), RLab()));
        g_mc_sigma_mass->GetYaxis()->SetRangeUser(0.0, 0.2); // Adjust range as needed
        g_mc_sigma_mass->GetXaxis()->SetLimits(-0.5, nMassCuts - 0.5);
        g_mc_sigma_mass->Draw("APL");
        g_data_sigma_mass->Draw("PL SAME");
        g_bkg_sigma_mass->Draw("PL SAME");
        g_combined_sigma_mass->Draw("PL SAME");

        // Add custom x-axis labels
        TAxis* xaxis_sigma = g_mc_sigma_mass->GetXaxis();
        for (int m = 0; m < nMassCuts; m++) {
            xaxis_sigma->SetBinLabel(xaxis_sigma->FindBin(m), massCutLabels[m]);
        }
        xaxis_sigma->SetLabelSize(0.05);

        // Legend
        TLegend* mass_sigma_legend = new TLegend(0.6, 0.7, 0.89, 0.89);
        mass_sigma_legend->SetBorderSize(0);
        mass_sigma_legend->AddEntry(g_mc_sigma_mass, SigLab(), "pl");
        mass_sigma_legend->AddEntry(g_data_sigma_mass, "Data", "pl");
        mass_sigma_legend->AddEntry(g_bkg_sigma_mass, BkgLab(), "pl");
        mass_sigma_legend->AddEntry(g_combined_sigma_mass, "Combined MC (Signal+Bkg)", "pl");
        mass_sigma_legend->Draw();

        AtlasUtils::DrawAtlasLabel(0.18, 0.85);
        TLatex mass_sigma_text;
        mass_sigma_text.SetTextSize(0.035);
        mass_sigma_text.DrawLatexNDC(0.2, 0.76, CurrentPtBinLabel(bin));

        mass_sigma_canvas->SetGridx();
        mass_sigma_canvas->SaveAs(Form(((RunConfig::PlotDir() + "/%s/Response_Fit_pTbins_Mcuts_")
                                        + CurrentChannel() + "/%s/R_DB_sigma_vs_masscut_ptbin%d.pdf")
                                           .c_str(),
                                       cfg.plot_subdir.c_str(), wp_label, bin));

        delete mass_sigma_canvas;
        delete g_mc_sigma_mass;
        delete g_data_sigma_mass;
        delete g_bkg_sigma_mass;
        delete g_combined_sigma_mass;
        delete mass_sigma_legend;

        // Create data/MC ratio comparison across mass cuts
        TCanvas* mass_ratio_canvas = new TCanvas(
            Form("mass_ratio_canvas_ptbin%d_%s", bin, wp_label),
            Form("Data/MC Ratio vs Mass Cut (pT bin: %s, WP: %s)", CurrentPtBinLabel(bin), wp_label), 800,
            600);
        mass_ratio_canvas->cd();

        // Calculate ratios for each mass cut
        double data_mc_ratio_values[MAX_MASS_CUTS], data_mc_ratio_errors[MAX_MASS_CUTS];
        double data_combined_ratio_values[MAX_MASS_CUTS], data_combined_ratio_errors[MAX_MASS_CUTS];
        double bkg_mc_ratio_values[MAX_MASS_CUTS], bkg_mc_ratio_errors[MAX_MASS_CUTS];

        for (int m = 0; m < nMassCuts; m++) {
            // Data/MC ratio
            if (mc_mean_values[m] != 0) {
                data_mc_ratio_values[m] = data_mean_values[m] / mc_mean_values[m];
                data_mc_ratio_errors[m] = data_mc_ratio_values[m]
                                          * std::sqrt(std::pow(data_mean_errors[m] / data_mean_values[m], 2)
                                                      + std::pow(mc_mean_errors[m] / mc_mean_values[m], 2));
            } else {
                data_mc_ratio_values[m] = 0.0;
                data_mc_ratio_errors[m] = 0.0;
            }

            // Data/Combined ratio
            if (combined_mean_values[m] != 0) {
                data_combined_ratio_values[m] = data_mean_values[m] / combined_mean_values[m];
                data_combined_ratio_errors[m] =
                    data_combined_ratio_values[m]
                    * std::sqrt(std::pow(data_mean_errors[m] / data_mean_values[m], 2)
                                + std::pow(combined_mean_errors[m] / combined_mean_values[m], 2));
            } else {
                data_combined_ratio_values[m] = 0.0;
                data_combined_ratio_errors[m] = 0.0;
            }

            // Background/MC ratio
            if (mc_mean_values[m] != 0) {
                bkg_mc_ratio_values[m] = bkg_mean_values[m] / mc_mean_values[m];
                bkg_mc_ratio_errors[m] = bkg_mc_ratio_values[m]
                                         * std::sqrt(std::pow(bkg_mean_errors[m] / bkg_mean_values[m], 2)
                                                     + std::pow(mc_mean_errors[m] / mc_mean_values[m], 2));
            } else {
                bkg_mc_ratio_values[m] = 0.0;
                bkg_mc_ratio_errors[m] = 0.0;
            }
        }

        // Create graphs
        TGraphErrors* g_data_mc_ratio_mass =
            new TGraphErrors(nMassCuts, massCutX, data_mc_ratio_values, nullptr, data_mc_ratio_errors);
        TGraphErrors* g_data_combined_ratio_mass = new TGraphErrors(
            nMassCuts, massCutX, data_combined_ratio_values, nullptr, data_combined_ratio_errors);
        TGraphErrors* g_bkg_mc_ratio_mass =
            new TGraphErrors(nMassCuts, massCutX, bkg_mc_ratio_values, nullptr, bkg_mc_ratio_errors);

        // Style graphs
        g_data_mc_ratio_mass->SetMarkerStyle(20);
        g_data_mc_ratio_mass->SetMarkerColor(kRed);
        g_data_mc_ratio_mass->SetLineColor(kRed);
        g_data_mc_ratio_mass->SetMarkerSize(1.2);

        g_data_combined_ratio_mass->SetMarkerStyle(21);
        g_data_combined_ratio_mass->SetMarkerColor(kBlack);
        g_data_combined_ratio_mass->SetLineColor(kBlack);
        g_data_combined_ratio_mass->SetMarkerSize(1.2);

        g_bkg_mc_ratio_mass->SetMarkerStyle(22);
        g_bkg_mc_ratio_mass->SetMarkerColor(kBlue);
        g_bkg_mc_ratio_mass->SetLineColor(kBlue);
        g_bkg_mc_ratio_mass->SetMarkerSize(1.2);

        // Draw graphs
        g_data_mc_ratio_mass->SetTitle(
            Form("Ratio vs Mass Cut for %s;Mass Cut;Ratio", CurrentPtBinLabel(bin)));
        g_data_mc_ratio_mass->GetYaxis()->SetRangeUser(0.95, 1.05); // Adjust range as needed
        g_data_mc_ratio_mass->GetXaxis()->SetLimits(-0.5, nMassCuts - 0.5);
        g_data_mc_ratio_mass->Draw("APL");
        g_data_combined_ratio_mass->Draw("PL SAME");
        g_bkg_mc_ratio_mass->Draw("PL SAME");

        // Add reference line at 1
        TLine* ratio_line = new TLine(-0.5, 1.0, nMassCuts - 0.5, 1.0);
        ratio_line->SetLineStyle(2);
        ratio_line->SetLineColor(kRed);
        ratio_line->Draw();

        // Add custom x-axis labels
        TAxis* xaxis_ratio = g_data_mc_ratio_mass->GetXaxis();
        for (int m = 0; m < nMassCuts; m++) {
            xaxis_ratio->SetBinLabel(xaxis_ratio->FindBin(m), massCutLabels[m]);
        }
        xaxis_ratio->SetLabelSize(0.05);

        // Legend
        TLegend* mass_ratio_legend = new TLegend(0.6, 0.7, 0.89, 0.89);
        mass_ratio_legend->SetBorderSize(0);
        mass_ratio_legend->AddEntry(g_data_mc_ratio_mass, "Data/Signal MC", "pl");
        mass_ratio_legend->AddEntry(g_data_combined_ratio_mass, "Data/Combined MC", "pl");
        mass_ratio_legend->AddEntry(g_bkg_mc_ratio_mass, "Multijets/Signal MC", "pl");
        mass_ratio_legend->AddEntry(ratio_line, "Ratio = 1", "l");
        mass_ratio_legend->Draw();

        AtlasUtils::DrawAtlasLabel(0.18, 0.85);
        TLatex mass_ratio_text;
        mass_ratio_text.SetTextSize(0.035);
        mass_ratio_text.DrawLatexNDC(0.2, 0.76, CurrentPtBinLabel(bin));

        mass_ratio_canvas->SetGridx();
        mass_ratio_canvas->SaveAs(Form(((RunConfig::PlotDir() + "/%s/Response_Fit_pTbins_Mcuts_")
                                        + CurrentChannel() + "/%s/R_DB_ratio_vs_masscut_ptbin%d.pdf")
                                           .c_str(),
                                       cfg.plot_subdir.c_str(), wp_label, bin));

        delete mass_ratio_canvas;
        delete g_data_mc_ratio_mass;
        delete g_data_combined_ratio_mass;
        delete g_bkg_mc_ratio_mass;
        delete ratio_line;
        delete mass_ratio_legend;
    }

    // Create a summary plot comparing all mass cuts in one figure (for data/combined MC ratio)
    TCanvas* all_ratio_canvas =
        new TCanvas("all_ratio_canvas", "Data/Combined MC Ratio for All pT and Mass Cuts", 1000, 700);
    all_ratio_canvas->cd();

    // We'll create a 2D plot to visualize all the ratios
    TH2F* ratio_map =
        new TH2F("ratio_map", Form("Data/Combined MC Ratio (WP: %s);Mass Cut;p_{T} Bin", wp_label), nMassCuts,
                 -0.5, nMassCuts - 0.5, CurrentNPtBins(), -0.5, CurrentNPtBins() - 0.5);

    // Fill the 2D plot
    for (int m = 0; m < nMassCuts; m++) {
        for (int bin = 0; bin < CurrentNPtBins(); bin++) {
            double ratio = data_means[m][bin] / combined_means[m][bin];
            ratio_map->SetBinContent(m + 1, bin + 1, ratio);

            // Calculate error
            double ratio_err =
                (combined_means[m][bin] != 0)
                    ? ratio
                          * std::sqrt(std::pow(data_mean_errs[m][bin] / data_means[m][bin], 2)
                                      + std::pow(combined_mean_errs[m][bin] / combined_means[m][bin], 2))
                    : 0.0;
            ratio_map->SetBinError(m + 1, bin + 1, ratio_err);
        }
    }

    // Style the plot
    ratio_map->SetStats(0);
    ratio_map->SetMarkerSize(1.5);

    // Set custom bin labels
    for (int m = 0; m < nMassCuts; m++) {
        ratio_map->GetXaxis()->SetBinLabel(m + 1, massCutLabels[m]);
    }

    for (int bin = 0; bin < CurrentNPtBins(); bin++) {
        // Simplify pT bin labels for readability
        TString simplifiedLabel =
            Form("%d-%d GeV", (int)CurrentPtBinsDisplay()[bin], (int)CurrentPtBinsDisplay()[bin + 1]);
        ratio_map->GetYaxis()->SetBinLabel(bin + 1, simplifiedLabel.Data());
    }

    ratio_map->GetXaxis()->SetLabelSize(0.04);
    ratio_map->GetYaxis()->SetLabelSize(0.04);

    gStyle->SetNumberContours(100);

    // Set a nice range centered at 1
    ratio_map->SetMinimum(0.95);
    ratio_map->SetMaximum(1.05);

    // Draw the plot
    ratio_map->Draw("COLZ TEXT");

    // Add ATLAS label
    AtlasUtils::DrawAtlasLabel(0.18, 0.85);

    all_ratio_canvas->SaveAs(Form(((RunConfig::PlotDir() + "/%s/Response_Fit_pTbins_Mcuts_")
                                   + CurrentChannel() + "/%s/R_DB_all_ratio_summary.pdf")
                                      .c_str(),
                                  cfg.plot_subdir.c_str(), wp_label));

    delete all_ratio_canvas;
    delete ratio_map;

    // Same map for background / signal MC
    TCanvas* bkg_ratio_canvas =
        new TCanvas("bkg_ratio_canvas", "Multijets/Signal MC Ratio for All pT and Mass Cuts", 1000, 700);
    bkg_ratio_canvas->cd();

    // Create a 2D plot for background to MC ratio
    TH2F* bkg_ratio_map =
        new TH2F("bkg_ratio_map", Form("Multijets/Signal MC Ratio (WP: %s);Mass Cut;p_{T} Bin", wp_label),
                 nMassCuts, -0.5, nMassCuts - 0.5, CurrentNPtBins(), -0.5, CurrentNPtBins() - 0.5);

    // Fill the 2D plot
    for (int m = 0; m < nMassCuts; m++) {
        for (int bin = 0; bin < CurrentNPtBins(); bin++) {
            double ratio = (mc_means[m][bin] != 0) ? bkg_means[m][bin] / mc_means[m][bin] : 0.0;
            bkg_ratio_map->SetBinContent(m + 1, bin + 1, ratio);

            // Calculate error
            double ratio_err = (mc_means[m][bin] != 0)
                                   ? ratio
                                         * std::sqrt(std::pow(bkg_mean_errs[m][bin] / bkg_means[m][bin], 2)
                                                     + std::pow(mc_mean_errs[m][bin] / mc_means[m][bin], 2))
                                   : 0.0;
            bkg_ratio_map->SetBinError(m + 1, bin + 1, ratio_err);
        }
    }

    // Style the plot
    bkg_ratio_map->SetStats(0);
    bkg_ratio_map->SetMarkerSize(1.5);

    // Set custom bin labels (same as previous)
    for (int m = 0; m < nMassCuts; m++) {
        bkg_ratio_map->GetXaxis()->SetBinLabel(m + 1, massCutLabels[m]);
    }

    for (int bin = 0; bin < CurrentNPtBins(); bin++) {
        TString simplifiedLabel =
            Form("%d-%d GeV", (int)CurrentPtBinsDisplay()[bin], (int)CurrentPtBinsDisplay()[bin + 1]);
        bkg_ratio_map->GetYaxis()->SetBinLabel(bin + 1, simplifiedLabel.Data());
    }

    bkg_ratio_map->GetXaxis()->SetLabelSize(0.04);
    bkg_ratio_map->GetYaxis()->SetLabelSize(0.04);

    // Same range as the data / combined MC map
    bkg_ratio_map->SetMinimum(0.95);
    bkg_ratio_map->SetMaximum(1.05);

    // Draw the plot
    bkg_ratio_map->Draw("COLZ TEXT");

    // Add ATLAS label
    AtlasUtils::DrawAtlasLabel(0.18, 0.85);

    bkg_ratio_canvas->SaveAs(Form(((RunConfig::PlotDir() + "/%s/Response_Fit_pTbins_Mcuts_")
                                   + CurrentChannel() + "/%s/R_DB_bkg_ratio_summary.pdf")
                                      .c_str(),
                                  cfg.plot_subdir.c_str(), wp_label));

    delete bkg_ratio_canvas;
    delete bkg_ratio_map;

    // JER summary in Summary/: sigma, R_JES and R_JER vs pT, heat maps, text and
    // console tables.
    {
        std::string jerDir = Form(
            ((RunConfig::PlotDir() + "/%s/Response_Fit_pTbins_Mcuts_") + CurrentChannel() + "/%s/Summary/")
                .c_str(),
            cfg.plot_subdir.c_str(), wp_label);
        PathBuilder::MakeDir(jerDir);

        Int_t colors_mass[MAX_MASS_CUTS] = {kRed, kBlue, kGreen + 2};
        Int_t markers_mass[MAX_MASS_CUTS] = {20, 21, 22};
        double offsets[MAX_MASS_CUTS] = {-8.0, 0.0, 8.0};

        // One quantity vs pT, one marker set per mass bin. Cells without a
        // valid truth fit are drawn at 0.
        auto drawSigmaVsPt = [&](const std::string& canvasName, const std::string& yTitle,
                                 const std::string& saveName, double vals[MAX_MASS_CUTS][MAX_PT_BINS],
                                 double errs[MAX_MASS_CUTS][MAX_PT_BINS], double ylo = 0.0,
                                 double yhi = 0.35) {
            TCanvas* c = new TCanvas(canvasName.c_str(), canvasName.c_str(), 800, 600);
            c->SetLeftMargin(0.14);
            c->SetBottomMargin(0.14);
            c->SetGridy();

            TH1F* frame =
                new TH1F((canvasName + "_fr").c_str(), "", CurrentNPtBins(), CurrentPtBinsDisplay());
            frame->SetMinimum(ylo);
            frame->SetMaximum(yhi);
            frame->GetXaxis()->SetTitle("p_{T} [GeV]");
            frame->GetYaxis()->SetTitle(yTitle.c_str());
            frame->GetYaxis()->SetTitleOffset(1.4);
            frame->SetStats(0);
            frame->Draw("AXIS");

            TLegend* leg = new TLegend(0.60, 0.70, 0.92, 0.88);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            leg->SetTextSize(0.030);

            for (int m = 0; m < nMassCuts; m++) {
                double x[MAX_PT_BINS], y[MAX_PT_BINS], xe[MAX_PT_BINS], ye[MAX_PT_BINS];
                for (int b = 0; b < CurrentNPtBins(); b++) {
                    x[b] = ptBinCenters[b] + offsets[m];
                    y[b] = vals[m][b];
                    xe[b] = 0;
                    ye[b] = std::abs(errs[m][b]);
                }
                TGraphErrors* gr = new TGraphErrors(CurrentNPtBins(), x, y, xe, ye);
                gr->SetMarkerColor(colors_mass[m]);
                gr->SetLineColor(colors_mass[m]);
                gr->SetMarkerStyle(markers_mass[m]);
                gr->SetMarkerSize(1.2);
                gr->SetLineWidth(2);
                gr->Draw("P SAME");
                leg->AddEntry(gr, massCutLabels[m], "p");
            }
            leg->Draw();

            TLegend* procLeg = new TLegend(0.15, 0.82, 0.45, 0.88);
            procLeg->SetBorderSize(0);
            procLeg->SetFillStyle(0);
            procLeg->SetTextSize(0.030);
            procLeg->AddEntry((TObject*)0, SigLab(), "");
            procLeg->Draw();

            AtlasUtils::DrawAtlasLabel(0.18, 0.85);
            c->SaveAs((jerDir + saveName).c_str());
            delete frame;
            delete c;
            delete leg;
            delete procLeg;
        };

        // One quantity as a (mass bin, pT bin) heat map
        auto drawHeatmap = [&](const std::string& name, const std::string& zTitle,
                               const std::string& saveName, double vals[MAX_MASS_CUTS][MAX_PT_BINS],
                               double errs[MAX_MASS_CUTS][MAX_PT_BINS]) {
            TCanvas* c = new TCanvas(name.c_str(), name.c_str(), 800, 600);
            c->SetTopMargin(0.10);
            c->SetRightMargin(0.16);
            c->SetLeftMargin(0.14);
            c->SetBottomMargin(0.12);

            TH2F* h2 = new TH2F(("h2_" + name).c_str(), Form(";Mass Bin;p_{T} Bin"), nMassCuts, -0.5,
                                nMassCuts - 0.5, CurrentNPtBins(), -0.5, CurrentNPtBins() - 0.5);

            for (int m = 0; m < nMassCuts; m++)
                for (int b = 0; b < CurrentNPtBins(); b++) {
                    h2->SetBinContent(m + 1, b + 1, vals[m][b]);
                    h2->SetBinError(m + 1, b + 1, errs[m][b]);
                }

            for (int m = 0; m < nMassCuts; m++) h2->GetXaxis()->SetBinLabel(m + 1, massCutLabels[m]);
            for (int b = 0; b < CurrentNPtBins(); b++)
                h2->GetYaxis()->SetBinLabel(b + 1, Form("%d-%d GeV", (int)CurrentPtBinsDisplay()[b],
                                                        (int)CurrentPtBinsDisplay()[b + 1]));

            h2->SetStats(0);
            h2->SetMarkerSize(1.4);
            h2->GetXaxis()->SetLabelSize(0.045);
            h2->GetYaxis()->SetLabelSize(0.045);
            gStyle->SetNumberContours(100);
            gStyle->SetPaintTextFormat("4.4f");

            // Colour range from the non-zero cells, plus 10% padding
            double minV = 1e9, maxV = -1e9;
            for (int m = 0; m < nMassCuts; m++)
                for (int b = 0; b < CurrentNPtBins(); b++) {
                    if (vals[m][b] != 0) {
                        minV = std::min(minV, vals[m][b]);
                        maxV = std::max(maxV, vals[m][b]);
                    }
                }
            if (minV < 1e9) {
                double pad = std::max(0.01, (maxV - minV) * 0.1);
                h2->SetMinimum(minV - pad);
                h2->SetMaximum(maxV + pad);
            }

            h2->Draw("COLZ TEXT");
            gPad->Update();

            TPaletteAxis* pal = (TPaletteAxis*)h2->GetListOfFunctions()->FindObject("palette");
            if (pal) {
                pal->SetX1NDC(0.85);
                pal->SetX2NDC(0.88);
                TLatex plab;
                plab.SetNDC();
                plab.SetTextSize(0.030);
                plab.SetTextAlign(22);
                plab.SetTextAngle(90);
                plab.DrawLatex(0.96, 0.5, zTitle.c_str());
            }
            AtlasUtils::DrawAtlasLabel(0.2, 0.92);
            c->SaveAs((jerDir + saveName).c_str());
            delete h2;
            delete c;
        };

        // sigma, R_JES and R_JER vs pT
        drawSigmaVsPt("sig_RDB_data", std::string("#sigma(") + RLab() + ")_{data}",
                      "sigma_RDB_data_vs_pT_" + std::string(wp_label) + ".pdf", jer_sigma_RDB_data,
                      jer_sigma_RDB_data_err);

        drawSigmaVsPt("sig_RDB_MC", std::string("#sigma(") + RLab() + ")_{MC}",
                      "sigma_RDB_MC_vs_pT_" + std::string(wp_label) + ".pdf", jer_sigma_RDB_MC,
                      jer_sigma_RDB_MC_err);

        drawSigmaVsPt("sig_truth", "#sigma(R_{truth,DB})_{MC}",
                      "sigma_truth_vs_pT_" + std::string(wp_label) + ".pdf", jer_sigma_truth,
                      jer_sigma_truth_err);

        drawSigmaVsPt("sig_reco_truth", "#sigma(p_{T}^{reco}/p_{T}^{truth})_{MC}",
                      "sigma_reco_truth_vs_pT_" + std::string(wp_label) + ".pdf", jer_sigma_reco_truth,
                      jer_sigma_reco_truth_err);

        drawSigmaVsPt("sig_insitu", "#sigma_{in situ}",
                      "sigma_insitu_vs_pT_" + std::string(wp_label) + ".pdf", jer_sigma_insitu,
                      jer_sigma_insitu_err);

        drawSigmaVsPt("sig_MC_insitu", "#sigma^{MC}_{in situ}",
                      "sigma_MC_insitu_vs_pT_" + std::string(wp_label) + ".pdf", jer_sigma_MC_insitu,
                      jer_sigma_MC_insitu_err);

        drawSigmaVsPt("sig_MC", "#sigma^{MC}", "sigma_MC_vs_pT_" + std::string(wp_label) + ".pdf",
                      jer_sigma_MC, jer_sigma_MC_err);
        drawSigmaVsPt("R_JES",
                      std::string("R_{JES} = #LT ") + RLab() + " #GT_{data} / #LT " + RLab() + " #GT_{MC}",
                      "R_JES_vs_pT_" + std::string(wp_label) + ".pdf", jes_ratio, jes_ratio_err, 0.95, 1.05);

        drawHeatmap("hm_R_JES",
                    std::string("R_{JES} = #LT ") + RLab() + " #GT_{data} / #LT " + RLab() + " #GT_{MC}",
                    "heatmap_R_JES_" + std::string(wp_label) + ".pdf", jes_ratio, jes_ratio_err);

        drawSigmaVsPt("R_JER", "R_{JER} = #sigma_{insitu}/#sigma^{MC}_{insitu}",
                      "R_JER_vs_pT_" + std::string(wp_label) + ".pdf", jer_ratio, jer_ratio_err, 0.5, 1.5);

        drawHeatmap("hm_R_JER", "R_{JER} = #sigma_{insitu}/#sigma^{MC}_{insitu}",
                    "heatmap_R_JER_" + std::string(wp_label) + ".pdf", jer_ratio, jer_ratio_err);

        // R_JES (top) and R_JER (bottom) vs pT, one figure per WP. Only cells
        // with a valid truth fit are drawn.
        {
            TCanvas* c = new TCanvas(Form("c_calib_%s", wp_label), "JES + JER Calibration Summary", 800, 800);

            // Top pad: R_JES vs pT
            TPad* p1 = new TPad("p1_jes", "", 0, 0.50, 1, 1.0);
            p1->SetTopMargin(0.08);
            p1->SetBottomMargin(0.02);
            p1->SetLeftMargin(0.14);
            p1->SetRightMargin(0.05);
            p1->SetGridy();
            p1->Draw();

            // Bottom pad: R_JER vs pT
            TPad* p2 = new TPad("p2_jer", "", 0, 0.0, 1, 0.50);
            p2->SetTopMargin(0.02);
            p2->SetBottomMargin(0.18);
            p2->SetLeftMargin(0.14);
            p2->SetRightMargin(0.05);
            p2->SetGridy();
            p2->Draw();

            // One ratio panel, used for both pads. The array parameters must
            // have the compile-time size of the caller's arrays; only the first
            // CurrentNPtBins() entries are used.
            auto drawRatioPanel = [&](TPad* pad, double vals[MAX_MASS_CUTS][MAX_PT_BINS],
                                      double errs[MAX_MASS_CUTS][MAX_PT_BINS], const char* yTitle, double ylo,
                                      double yhi, bool drawXLabels) {
                pad->cd();
                TH1F* fr = new TH1F(Form("fr_calib_%s_%p", wp_label, (void*)pad), "", CurrentNPtBins(),
                                    CurrentPtBinsDisplay());
                fr->SetMinimum(ylo);
                fr->SetMaximum(yhi);
                fr->SetStats(0);
                fr->GetYaxis()->SetTitle(yTitle);
                fr->GetYaxis()->SetTitleSize(0.06);
                fr->GetYaxis()->SetTitleOffset(1.1);
                fr->GetYaxis()->SetLabelSize(0.05);
                fr->GetYaxis()->SetNdivisions(505);
                if (drawXLabels) {
                    fr->GetXaxis()->SetTitle("p_{T} [GeV]");
                    fr->GetXaxis()->SetTitleSize(0.06);
                    fr->GetXaxis()->SetLabelSize(0.05);
                } else {
                    fr->GetXaxis()->SetLabelSize(0);
                    fr->GetXaxis()->SetTitleSize(0);
                }
                fr->Draw("AXIS");

                // Reference line at 1
                double xMin = CurrentPtBinsDisplay()[0], xMax = CurrentPtBinsDisplay()[MAX_PT_BINS];
                TLine* ref = new TLine(xMin, 1.0, xMax, 1.0);
                ref->SetLineStyle(2);
                ref->SetLineColor(kGray + 2);
                ref->SetLineWidth(2);
                ref->Draw("SAME");

                // One TGraphErrors per mass cut, with horizontal offsets so they don't overlap
                for (int m = 0; m < nMassCuts; m++) {
                    double x[MAX_PT_BINS], y[MAX_PT_BINS], xe[MAX_PT_BINS], ye[MAX_PT_BINS];
                    int nValid = 0;
                    for (int b = 0; b < CurrentNPtBins(); b++) {
                        if (!jer_valid[m][b]) continue;
                        x[nValid] = ptBinCenters[b] + offsets[m];
                        y[nValid] = vals[m][b];
                        xe[nValid] = 0;
                        ye[nValid] = std::abs(errs[m][b]);
                        ++nValid;
                    }
                    if (nValid == 0) continue;
                    TGraphErrors* gr = new TGraphErrors(nValid, x, y, xe, ye);
                    gr->SetMarkerColor(colors_mass[m]);
                    gr->SetLineColor(colors_mass[m]);
                    gr->SetMarkerStyle(markers_mass[m]);
                    gr->SetMarkerSize(1.3);
                    gr->SetLineWidth(2);
                    gr->Draw("P SAME");
                }
            };

            // Top: R_JES (expected ~1, narrow range)
            drawRatioPanel(p1, jes_ratio, jes_ratio_err,
                           Form("R_{JES} = #LT %s #GT_{data} / #LT %s #GT_{MC}", RLab(), RLab()), 0.95, 1.05,
                           /*drawXLabels=*/false);

            // ATLAS label + WP info on top pad only
            p1->cd();
            AtlasUtils::DrawAtlasLabel(0.2, 0.85);
            TLatex lt;
            lt.SetNDC();
            lt.SetTextSize(0.05);
            lt.DrawLatex(0.18, 0.78, Form("%s  WP: %s", cfg.label.c_str(), wp_label));

            // Mass-cut legend on top pad
            TLegend* leg = new TLegend(0.60, 0.70, 0.92, 0.90);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            leg->SetTextSize(0.045);
            for (int m = 0; m < nMassCuts; m++) {
                TGraph* dummy = new TGraph();
                dummy->SetMarkerColor(colors_mass[m]);
                dummy->SetMarkerStyle(markers_mass[m]);
                dummy->SetMarkerSize(1.3);
                dummy->SetLineColor(colors_mass[m]);
                dummy->SetLineWidth(2);
                leg->AddEntry(dummy, massCutLabels[m], "lp");
            }
            leg->Draw();

            // Bottom: R_JER, on a wider range than R_JES
            drawRatioPanel(p2, jer_ratio, jer_ratio_err, "R_{JER} = #sigma_{insitu} / #sigma^{MC}_{insitu}",
                           0.5, 1.5, /*drawXLabels=*/true);

            c->SaveAs((jerDir + "Calibration_summary_" + std::string(wp_label) + ".pdf").c_str());
            delete c;
        }

        drawSigmaVsPt("delta_closure", "#Delta#sigma_{closure}",
                      "delta_closure_vs_pT_" + std::string(wp_label) + ".pdf", jer_delta_closure,
                      jer_delta_closure_err, -0.5, 0.5);

        // All sigma quantities on one canvas per mass bin
        for (int m = 0; m < nMassCuts; m++) {
            std::string cname = Form("jer_overlay_m%d_%s", m, wp_label);
            TCanvas* cov = new TCanvas(cname.c_str(), cname.c_str(), 900, 650);
            cov->SetLeftMargin(0.14);
            cov->SetBottomMargin(0.14);
            cov->SetGridy();

            TH1F* fr = new TH1F((cname + "_fr").c_str(), "", CurrentNPtBins(), CurrentPtBinsDisplay());
            fr->SetMinimum(0);
            fr->SetMaximum(0.40);
            fr->GetXaxis()->SetTitle("p_{T} [GeV]");
            fr->GetYaxis()->SetTitle("#sigma");
            fr->GetYaxis()->SetTitleOffset(1.4);
            fr->SetStats(0);
            fr->Draw("AXIS");

            struct PlotDef {
                double* v;
                double* e;
                int col;
                int mkr;
                const char* label;
                double offset;
            };
            std::string lblRDBdata = std::string("#sigma(") + RLab() + ")_{data}",
                        lblRDBmc = std::string("#sigma(") + RLab() + ")_{MC}";
            PlotDef plots[] = {
                {jer_sigma_RDB_data[m], jer_sigma_RDB_data_err[m], kBlack, 20, lblRDBdata.c_str(), -12.5},
                {jer_sigma_RDB_MC[m], jer_sigma_RDB_MC_err[m], kRed, 21, lblRDBmc.c_str(), -7.5},
                {jer_sigma_truth[m], jer_sigma_truth_err[m], kAzure + 2, 22, "#sigma(R_{truth})_{MC}", -2.5},
                {jer_sigma_reco_truth[m], jer_sigma_reco_truth_err[m], kOrange + 7, 23,
                 "#sigma^{MC} (p_{T}^{reco}/p_{T}^{truth})", 2.5},
                {jer_sigma_insitu[m], jer_sigma_insitu_err[m], kViolet - 1, 29, "#sigma_{in situ}", 7.5},
                {jer_sigma_MC_insitu[m], jer_sigma_MC_insitu_err[m], kGreen + 2, 33, "#sigma^{MC}_{in situ}",
                 12.5},
            };

            TLegend* leg = new TLegend(0.50, 0.55, 0.92, 0.88);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            leg->SetTextSize(0.028);

            for (auto& p : plots) {
                double x[MAX_PT_BINS], xe[MAX_PT_BINS];
                for (int b = 0; b < CurrentNPtBins(); b++) {
                    x[b] = ptBinCenters[b] + p.offset;
                    xe[b] = 0;
                }
                TGraphErrors* gr = new TGraphErrors(CurrentNPtBins(), x, p.v, xe, p.e);
                gr->SetMarkerColor(p.col);
                gr->SetLineColor(p.col);
                gr->SetMarkerStyle(p.mkr);
                gr->SetMarkerSize(1.2);
                gr->SetLineWidth(2);
                gr->Draw("P SAME");
                leg->AddEntry(gr, p.label, "p");
            }
            leg->Draw();
            AtlasUtils::DrawAtlasLabel(0.18, 0.85);

            TLatex lt;
            lt.SetNDC();
            lt.SetTextSize(0.035);
            lt.DrawLatex(0.16, 0.85, Form("Mass: %s", massCutLabels[m]));

            cov->SaveAs((jerDir + Form("JER_all_sigmas_m%d_%s.pdf", m, wp_label)).c_str());
            delete fr;
            delete cov;
            delete leg;
        }

        // Heat maps
        drawHeatmap("hm_sig_data", std::string("#sigma(") + RLab() + ")_{data}",
                    "heatmap_sigma_RDB_data_" + std::string(wp_label) + ".pdf", jer_sigma_RDB_data,
                    jer_sigma_RDB_data_err);
        drawHeatmap("hm_sig_MC", std::string("#sigma(") + RLab() + ")_{MC}",
                    "heatmap_sigma_RDB_MC_" + std::string(wp_label) + ".pdf", jer_sigma_RDB_MC,
                    jer_sigma_RDB_MC_err);
        drawHeatmap("hm_sig_truth", "#sigma(R_{truth})_{MC}",
                    "heatmap_sigma_truth_" + std::string(wp_label) + ".pdf", jer_sigma_truth,
                    jer_sigma_truth_err);
        drawHeatmap("hm_sig_recotruth", "#sigma(p_{T}^{reco}/p_{T}^{truth})",
                    "heatmap_sigma_reco_truth_" + std::string(wp_label) + ".pdf", jer_sigma_reco_truth,
                    jer_sigma_reco_truth_err);
        drawHeatmap("hm_sig_insitu", "#sigma_{in situ}",
                    "heatmap_sigma_insitu_" + std::string(wp_label) + ".pdf", jer_sigma_insitu,
                    jer_sigma_insitu_err);
        drawHeatmap("hm_sig_MCinsitu", "#sigma^{MC}_{in situ}",
                    "heatmap_sigma_MC_insitu_" + std::string(wp_label) + ".pdf", jer_sigma_MC_insitu,
                    jer_sigma_MC_insitu_err);
        drawHeatmap("hm_sig_MC93", "#sigma^{MC}", "heatmap_sigma_MC_" + std::string(wp_label) + ".pdf",
                    jer_sigma_MC, jer_sigma_MC_err);
        drawHeatmap("hm_delta_clos", "#Delta#sigma_{closure}",
                    "heatmap_delta_closure_" + std::string(wp_label) + ".pdf", jer_delta_closure,
                    jer_delta_closure_err);

        // Text summary
        {
            std::string txtPath = jerDir + "JER_full_summary_" + std::string(wp_label) + ".txt";
            std::ofstream txt(txtPath);
            txt << "# JER Full Summary - " << cfg.key << "  WP: " << wp_label << "\n";
            txt << "# Truth match cuts: none (R_truth_DB > 0 only)\n";
            txt << "#\n";
            txt << std::left << std::setw(14) << "# MassBin" << std::setw(16) << "pTBin" << std::right
                << std::setw(14) << "sigma_RDB_data" << std::setw(10) << "±err" << std::setw(14) << "sigma_RDB_MC"
                << std::setw(10) << "±err" << std::setw(14) << "sigma_truth" << std::setw(10) << "±err"
                << std::setw(14) << "sigma_reco/tr" << std::setw(10) << "±err" << std::setw(14) << "sigma_insitu"
                << std::setw(10) << "±err" << std::setw(14) << "sigma^MC_insitu" << std::setw(10) << "±err"
                << std::setw(14) << "sigma^MC" << std::setw(10) << "±err" << std::setw(14) << "Deltasigma_clos"
                << std::setw(10) << "±err" << std::setw(14) << "R_JES" << std::setw(10) << "±err"
                << std::setw(14) << "R_JER" << std::setw(10) << "±err"
                << "\n";
            txt << std::string(234, '-') << "\n";

            for (int m = 0; m < nMassCuts; m++) {
                for (int b = 0; b < CurrentNPtBins(); b++) {
                    txt << std::left << std::setw(14) << fileSystemMassCutLabels[m] << std::setw(16)
                        << Form("%d-%dGeV", (int)CurrentPtBinsDisplay()[b],
                                (int)CurrentPtBinsDisplay()[b + 1]);
                    if (jer_valid[m][b]) {
                        txt << std::right << std::fixed << std::setprecision(5) << std::setw(14)
                            << jer_sigma_RDB_data[m][b] << std::setw(10) << jer_sigma_RDB_data_err[m][b]
                            << std::setw(14) << jer_sigma_RDB_MC[m][b] << std::setw(10)
                            << jer_sigma_RDB_MC_err[m][b] << std::setw(14) << jer_sigma_truth[m][b]
                            << std::setw(10) << jer_sigma_truth_err[m][b] << std::setw(14)
                            << jer_sigma_reco_truth[m][b] << std::setw(10) << jer_sigma_reco_truth_err[m][b]
                            << std::setw(14) << jer_sigma_insitu[m][b] << std::setw(10)
                            << jer_sigma_insitu_err[m][b] << std::setw(14) << jer_sigma_MC_insitu[m][b]
                            << std::setw(10) << jer_sigma_MC_insitu_err[m][b] << std::setw(14)
                            << jer_sigma_MC[m][b] << std::setw(10) << jer_sigma_MC_err[m][b] << std::setw(14)
                            << jer_delta_closure[m][b] << std::setw(10) << jer_delta_closure_err[m][b]
                            << std::setw(14) << jes_ratio[m][b] << std::setw(10) << jes_ratio_err[m][b]
                            << std::setw(14) << jer_ratio[m][b] << std::setw(10) << jer_ratio_err[m][b];
                    } else {
                        txt << "  [NO TRUTH FIT]";
                    }
                    txt << "\n";
                }
            }
            txt.close();
            Log::info(("Rdb" + CurrentChannel()).c_str(), "JER summary: " + txtPath);
        }

        // Console table
        std::cout << "\n";
        std::cout
            << "||===========================================================================================================||\n";
        std::cout << "||  JER Summary - " << std::left << std::setw(15) << cfg.key << " WP: " << std::setw(10)
                  << wp_label << std::right << std::setw(63) << "||\n";
        std::cout
            << "||===========================================================================================================||\n";
        std::cout << "|| " << std::left << std::setw(12) << "Mass" << std::setw(12) << "pT" << std::setw(12)
                  << "sigma_data" << std::setw(12) << "sigma_truth" << std::setw(12) << "sigma_insitu" << std::setw(12)
                  << "sigma^MC(9.3)" << std::setw(12) << "sigma^MC_in" << std::setw(10) << "Deltaclos" << std::setw(12)
                  << "R_JES" << std::setw(12) << "R_JER" << std::right << "   ||\n";
        std::cout
            << "||===========================================================================================================||\n";
        for (int m = 0; m < nMassCuts; m++) {
            for (int b = 0; b < CurrentNPtBins(); b++) {
                std::cout << "|| " << std::left << std::setw(12) << fileSystemMassCutLabels[m] << std::setw(12)
                          << Form("%d-%d", (int)CurrentPtBinsDisplay()[b],
                                  (int)CurrentPtBinsDisplay()[b + 1]);
                if (jer_valid[m][b]) {
                    std::cout << std::fixed << std::setprecision(4) << std::setw(12)
                              << jer_sigma_RDB_data[m][b] << std::setw(12) << jer_sigma_truth[m][b]
                              << std::setw(12) << jer_sigma_insitu[m][b] << std::setw(12)
                              << jer_sigma_MC[m][b] << std::setw(12) << jer_sigma_MC_insitu[m][b]
                              << std::setprecision(3) << std::setw(10) << jer_delta_closure[m][b]
                              << std::setprecision(4) << std::setw(12) << jes_ratio[m][b]
                              << std::setprecision(4) << std::setw(12) << jer_ratio[m][b];
                } else {
                    std::cout << std::setw(70) << "[no truth fit]";
                }
                std::cout << std::right << "   ║\n";
            }
        }
        std::cout
            << "||===========================================================================================================||\n";
    }

    // Binning recommendation, printed under "[BINNING-REC]" and written to
    // binning_recommendation.txt, ready to paste into VariantConfig.h.
    // For each cell, from the filled combined-MC histograms:
    //   fit-friendly  ~EASYBJETS_BINNING_CORE_BINS bins across the +/-2 sigma core (recommended)
    //   Scott         3.49 sigma N_eff^(-1/3), density-optimal (reference only)
    // reco comes from the Response histogram (nBinsR); truth takes the larger count of
    // truth/ref and reco/truth (nBinsTruth). The pasted value per pT bin is the
    // maximum over the mass bins.
    {
        const int rMin = 10, rMax = 200, tMin = 10, tMax = 200;
        auto scott = [&](TH1F* h, int mn, int mx) -> int {
            return (h && h->GetEntries() > 20) ? BinningConfig::ComputeScott(h, 0.0, 2.0, mn, mx) : 0;
        };
        auto fitfr = [&](TH1F* h, int mn, int mx) -> int {
            return (h && h->GetEntries() > 20) ? BinningConfig::ComputeFitFriendly(h, 0.0, 2.0, mn, mx) : 0;
        };
        int scR[MAX_MASS_CUTS][MAX_PT_BINS] = {}, scT[MAX_MASS_CUTS][MAX_PT_BINS] = {};
        int ffR[MAX_MASS_CUTS][MAX_PT_BINS] = {}, ffT[MAX_MASS_CUTS][MAX_PT_BINS] = {};
        int maxR[MAX_PT_BINS] = {0}, maxT[MAX_PT_BINS] = {0};
        for (int m = 0; m < nMassCuts; m++)
            for (int b = 0; b < CurrentNPtBins(); b++) {
                scR[m][b] = scott(h_combined_MC_R_balance_ptbin[m][b], rMin, rMax);
                scT[m][b] = std::max(scott(h_comb_truth_R_balance_ptbin[m][b], tMin, tMax),
                                     scott(h_comb_reco_truth_ptbin[m][b], tMin, tMax));
                ffR[m][b] = fitfr(h_combined_MC_R_balance_ptbin[m][b], rMin, rMax);
                ffT[m][b] = std::max(fitfr(h_comb_truth_R_balance_ptbin[m][b], tMin, tMax),
                                     fitfr(h_comb_reco_truth_ptbin[m][b], tMin, tMax));
                if (ffR[m][b] > maxR[b]) maxR[b] = ffR[m][b];
                if (ffT[m][b] > maxT[b]) maxT[b] = ffT[m][b];
            }
        auto emit = [&](std::ostream& os) {
            os << "==== BINNING RECOMMENDATION -- " << cfg.key << "  WP " << wp_label << " ====\n";
            os << "  fit-friendly = ~N bins in +/-2sigma core (set EASYBJETS_BINNING_CORE_BINS, default 10);"
                  " Scott = density-optimal (reference).\n";
            for (int m = 0; m < nMassCuts; m++)
                for (int b = 0; b < CurrentNPtBins(); b++)
                    os << "  [" << massCutLabels[m] << "] " << CurrentPtBinLabel(b)
                       << " : reco(fit=" << ffR[m][b] << ", scott=" << scR[m][b] << ")"
                       << "  truth(fit=" << ffT[m][b] << ", scott=" << scT[m][b] << ")\n";
            os << "  paste into common/VariantConfig.h (fit-friendly, per-pT max over mass cuts; "
                  "truth often capped ~100):\n";
            os << "  configs[\"" << cfg.key << "\"].nBinsR_perPtBin     = {";
            for (int b = 0; b < CurrentNPtBins(); b++)
                os << maxR[b] << (b + 1 < CurrentNPtBins() ? ", " : "");
            os << "};\n";
            os << "  configs[\"" << cfg.key << "\"].nBinsTruth_perPtBin = {";
            for (int b = 0; b < CurrentNPtBins(); b++)
                os << maxT[b] << (b + 1 < CurrentNPtBins() ? ", " : "");
            os << "};\n";
        };
        std::cout << "\n+----------------------------- [BINNING-REC] -----------------------------+\n";
        emit(std::cout);
        std::cout << "+-------------------------------------------------------------------------+\n\n";
        std::string recPath = Form(((RunConfig::PlotDir() + "/%s/Response_Fit_pTbins_Mcuts_")
                                    + CurrentChannel() + "/%s/binning_recommendation.txt")
                                       .c_str(),
                                   cfg.plot_subdir.c_str(), wp_label);
        std::ofstream rec(recPath);
        if (rec.is_open()) {
            emit(rec);
            rec.close();
            Log::info(("Rdb" + CurrentChannel()).c_str(), "Binning recommendation written: " + recPath);
        }
    }

    // Clean up all dynamically allocated histograms for this working point
    for (int m = 0; m < nMassCuts; m++) {
        for (int i = 0; i < CurrentNPtBins(); i++) {
            delete h_mc_R_balance_ptbin[m][i];
            delete h_data_R_balance_ptbin[m][i];
            delete h_multijets_R_balance_ptbin[m][i];
            delete h_combined_MC_R_balance_ptbin[m][i];
            delete h_mc_truth_R_balance_ptbin[m][i];
            delete h_mc_reco_truth_ptbin[m][i];
            delete h_bkg_truth_R_balance_ptbin[m][i];
            delete h_bkg_reco_truth_ptbin[m][i];
            delete h_comb_truth_R_balance_ptbin[m][i];
            delete h_comb_reco_truth_ptbin[m][i];
        }
    }

    // Clean up extra files (SafeTFiles fdata/fmc/fmultijets auto-close on scope exit)
    for (int s = 0; s < SampleConfig::kNSamples; ++s)
        if (fextra[s]) {
            fextra[s]->Close();
            delete fextra[s];
        }

    std::cout << "All plots and summary for WP " << wp_label << " generated successfully!" << std::endl;
}

// variant: a key such as "DB_bJER10v00", or "ALL" for the 4 variants of the channel;
// channel: "DB" or "MJB"; massMode: "cuts" (3 mass bins) or "inclusive" (one bin, no mass cut).
inline void Rdb_pT_bins_MCut_Run(const std::string& variant, const std::string& channel,
                                 const std::string& massMode = "cuts") {
    CurrentChannel() = channel;
    MassInclusiveMode() = (massMode == "inclusive");
    ConfigureMassBins();
    AtlasUtils::LoadAtlasStyle();
    MassCutConfig::PrintBanner("Rdb_pT_bins_MCut");
    Log::info("Rdb_pT_bins_MCut",
              std::string("Mass mode: ")
                  + (MassInclusiveMode() ? "inclusive (1 bin, no mass cut)" : "cuts (3 mass bins)"));
    // Keep histograms out of the input files, so closing a file does not delete them
    TH1::AddDirectory(kFALSE);

    if (variant == "ALL") {
        for (const auto& key : GetPrimaryVariantKeys(CurrentChannel())) {
            Log::info("Rdb_pT_bins_MCut", "Running for: " + key);
            Rdb_pT_bins_MCut_Run(key, channel, massMode);
        }
        return;
    }

    VariantConfig cfg = GetVariant(variant);

    auto workingPoints = AnalysisWPConfig::StandardWPs();
    for (const auto& wp : workingPoints) {
        ProcessWorkingPoint(cfg, wp.c_str());
    }
    std::cout << "\n--- All working points processed ---" << std::endl;
}

#endif // RDB_PT_BINS_MCUT_H
