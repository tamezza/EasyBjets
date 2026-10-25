/**
 * ============================================================================
 *  Response_Fit_Inclusive.h - inclusive fit of the balance R, R_JES and R_JER
 * ============================================================================
 * ===========================================
 * @author:    Thiziri AMEZZA
 * @email:     thiziri.amezza@cern.ch
 * @date:      08 April 2026
 * ===========================================
 *
 * Shared by both channels: ZbbGamma/macros/Response_Fit_Inclusive.cxx
 * and ZbbJets/macros/Response_Fit_Inclusive.cxx call FitResponse_Run(variant,
 * "DB" or "MJB"). For each variant and working point, inclusive in pT:
 *
 *   1. Fill the jet mass and the balance R = pT(jet)/pT(ref) for data, the
 *      Z(bb) signal MC, the main background (gamma+jets for DB, multijets for
 *      MJB) and the optional minor backgrounds. With EASYBJETS_MASS_CUT set,
 *      only jets inside the mass window are kept, in every loop.
 *
 *   2. Take the mass-dependent background scale factor SF(m) from the cache
 *      written by Mass_Distributions (sf2_params_<wp>.txt), or refit it from
 *      the mass sidebands if the cache is missing or was fitted with other
 *      sidebands, and apply it event by event to the main background.
 *
 *   3. Fit the R peak with the iterative Gaussian fit: data first, then the
 *      signal MC, the background and the combined MC (signal + backgrounds). 
 *      The MC fits start from the final data window, but each then iterates 
 *      its own window (with EASYBJETS_SEED_MODE=raw the data window is not 
 *      used at all).
 *
 *   4. Fit the MC truth widths: sigma_truth, the width of the truth balance
 *      pT(truth jet)/pT(ref), and sigma^MC, the width of pT(reco)/pT(truth).
 *      From these and the R widths of data and combined MC:
 *        sigma_insitu    = sqrt(sigma_data^2  - sigma_truth^2)
 *        sigma^MC_insitu = sqrt(sigma_combMC^2 - sigma_truth^2)
 *        Delta_closure   = (sigma^MC - sigma^MC_insitu) / sigma_insitu
 *        R_JES = mu_data / mu_combMC
 *        R_JER = sigma_insitu / sigma^MC_insitu
 *      These are computed only when the MC has truth branches; otherwise the
 *      summary shows "[not computed]". If a difference of squares is
 *      negative, R_JER is not computed and a warning is printed.
 *
 *   5. Save the plots, the fit parameter files read by
 *      Response_Ratios_Inclusive.cxx, the JER summary tables, and a bin-count
 *      recommendation to paste into VariantConfig.h.
 *
 * The closure relation sigma^2(R) = sigma^2(R_truth) + sigma^2(reco/truth) only
 * holds if the three histograms see the same events, so the MC R histograms
 * used in the fits require a truth match with dR(reco, truth) < cfg.maxTruthDR.
 * The "_all" histograms have no truth requirement and are used for the data/MC
 * stack, which must describe the whole population that data sees.
 *
 * How to run (from ZbbGamma/macros for DB, from ZbbJets/macros for MJB):
 * 
 *   export EASYBJETS_RUN=run3 (run2 for 13 TeV)
 *   export EASYBJETS_DATA_DIR=$HOME/lpnatlas_data/ZbbGamma    (ZbbJets for MJB)
 *   export EASYBJETS_MASS_CUT=on       keep only 50 < m_J < 150 GeV in every loop
 *   root 'Response_Fit_Inclusive.cxx()'                 nominal DB (MJB in ZbbJets)
 *   root 'Response_Fit_Inclusive.cxx("DB_bJER10v00")'   one variant, full key
 *   root 'Response_Fit_Inclusive.cxx("MJB_MuonCorr")'   from ZbbJets/macros
 *   root 'Response_Fit_Inclusive.cxx("ALL")'            the four variants of the channel
 * 
 * Run Mass_Distributions.cxx first with the same options so the SF cache
 * exists (otherwise the SF is refitted here).
 * Then run Response_Ratios_Inclusive.cxx with the same options.
 *
 * Output: ../Run3Plots/<variant>/FitResponse_<path_tag>/<wp>/  (Run2Plots for run2,
 *   e.g. ../Run3Plots/DB/FitResponse_nominal/0p46/)
 *   {data,mc,bkg,combined,truth_ref,reco_truth}_fit_parameters.txt
 *   Fit_R_balance_*.pdf, JER_results.txt, JER_per_sample.txt,
 *   R_DB_summary_table.txt, binning_recommendation.txt,
 *   jet_mass_comparison.pdf, R_DB_stack_ratio_combined.pdf
 *
 * Options (environment variables, read once per ROOT session):
 *   EASYBJETS_RUN                    run2 | run3 (default)
 *   EASYBJETS_DATA_DIR               folder with NtupleSlim_Run2/3 (default ../../DataSets/<ZbbGamma|ZbbJets>)
 *   EASYBJETS_MASS_CUT               unset/off = no cut (default) | on = 50-150 GeV | lo:hi
 *   EASYBJETS_RESPONSE_FIT_LO        lowest R allowed in the fit window (default 0.3, negative = no floor)
 *   EASYBJETS_FIT_MODE               scan (default) | scan_dense | scan_log | legacy | legacy7grid
 *   EASYBJETS_FIT_NSIGMA             window half-width in legacy mode (default 2.0, capped at 1.5 here)
 *   EASYBJETS_FIT_NSIGMA_LIST        scan windows (default "1.0,1.5"; entries above 1.5 are dropped here)
 *   EASYBJETS_FIT_MODEL              gauss (default) | all | bukin | cb | dscb | dgauss; non-default
 *                                    models write to Run{2,3}Plots_<model>
 *   EASYBJETS_FIT_MODEL_NSIGMA_LO/HI window of the alternative models in sigma (default 3 / 2)
 *   EASYBJETS_FIT_ERR_SCALE          pdg = scale shown/saved errors by sqrt(chi2/ndf) when > 1 (default off)
 *   EASYBJETS_SEED_MODE              truncated (default) | raw = seed from raw mean +/- 1.5 RMS,
 *                                    ignoring the windows passed here
 *   EASYBJETS_RESPONSE_BINNING_MODE  static (default, cfg.nBinsR / nBinsTruth) | dynamic
 *   EASYBJETS_BINNING                bin algorithm in dynamic mode: fixed (default) | fd | scott |
 *                                    scottjes | sturges | peakaware | variable | fitfriendly
 *   EASYBJETS_BINNING_TRUTH          truth binning in dynamic mode: match (default) | auto | fixed | <algorithm>
 *   EASYBJETS_BINNING_MAx_BINS(_TRUTH) lower the bin-count cap in dynamic mode (default unset)
 *   EASYBJETS_BINNING_CORE_BINS      bins in the +/-2 sigma core for the recommendation (default 10)
 *   EASYBJETS_SIDEBAND1_LOW/HIGH, EASYBJETS_SIDEBAND2_LOW/MAx
 *                                    override the SF sidebands of every variant (default unset)
 *
 * ==========================================================================
 * @version:  1.0
 * ==========================================================================
 */
#ifndef FIT_RESPONSE_H
#define FIT_RESPONSE_H

#include "VariantConfig.h"
#include "CommonUtils.h"
#include "SFUtils.h"
#include "BinningConfig.h"
#include "MassCutConfig.h"
#include "Response_Fit_Common.h" // RDBFitResults + channel/pT-bin accessors (shared)
#include "Response_Fit_Core.h"   // shared response-fit core functions

#include <fstream>
#include <limits>

#include "TBox.h"

using namespace std;

// Output folders
inline std::string FitRdbDir(const VariantConfig& cfg) {
    return cfg.path_tag.empty() ? "FitResponse" : "FitResponse_" + cfg.path_tag;
}

inline std::string BasePlotPath(const VariantConfig& cfg, const char* wp) {
    return (RunConfig::PlotDir() + "/") + cfg.plot_subdir + "/" + FitRdbDir(cfg) + "/" + wp;
}

// Load the SF(m) cached by Mass_Distributions; refit it from the mass
// histograms only if the cache is missing or was made with other sidebands.
static SFUtils::LinearSF GetOrFitSF(TH1F* hData, TH1F* hBkg, const VariantConfig& cfg,
                                    const std::string& wp_raw, const std::string& plotPath) {
    std::string sfPath = PathBuilder::OutputDirWP(cfg, "Mass_Distributions", wp_raw) + "sf2_params_"
                         + PathBuilder::WPtoFileWP(wp_raw) + ".txt";

    // The sidebands are passed so a cache fitted with other windows (or an
    // old cache without the sideband stamp) is rejected and refitted below.
    SFUtils::LinearSF cached =
        SFUtils::LoadLinearSF(sfPath, cfg.sideband1_low, cfg.sideband1_high, cfg.sideband2_low,
                              static_cast<double>(cfg.sideband2_max));
    if (cached.valid) {
        Log::info("FitResponse", "Loaded cached SF from: " + sfPath);
        return cached;
    }

    SFUtils::LinearSF sf = SFUtils::FitMassDependentSF(
        hData, hBkg, cfg.sideband1_low, cfg.sideband1_high, cfg.sideband2_low, cfg.sideband2_max,
        plotPath + "/scale_factor_fit.pdf", cfg.key + " " + PathBuilder::WPtoFileWP(wp_raw));

    if (!sf.valid) Log::warning("FitResponse", "SF fit failed — using SF=1");
    return sf;
}

// Full analysis of one variant at one working point
inline void RunRdbAnalysis(const VariantConfig& cfg, const char* wp_raw) {
    TString fileWP_ts = PathBuilder::WPtoFileWP(std::string(wp_raw)).c_str();
    const char* wp_label = fileWP_ts.Data();
    TString basePath = BasePlotPath(cfg, wp_label).c_str();
    PathBuilder::MakeDir(std::string(basePath.Data()));

    Log::info("FitResponse", Form("--- Processing: %s WP %s ---", cfg.key.c_str(), wp_label));
    MassCutConfig::PrintBanner("FitResponse");

    const int nBins = cfg.nBinsR;

    // Open input files
    std::string dataPath = PathBuilder::SamplePath(cfg, "data", std::string(wp_raw));
    std::string zbbPath =
        PathBuilder::SamplePath(cfg, SampleConfig::FileNames(cfg)[SampleConfig::kZbb], std::string(wp_raw));
    std::string mjPath = PathBuilder::SamplePath(cfg, SampleConfig::FileNames(cfg)[SampleConfig::kMultijets],
                                                 std::string(wp_raw));

    auto scFileNames_rdb = SampleConfig::FileNames(cfg);
    auto scDisplayNames_rdb = SampleConfig::DisplayNames(cfg);
    auto scColors_rdb = SampleConfig::Colors();
    // Channel-aware signal/background legend labels (DB: Z(bb)+gamma / gamma+jets,
    // MJB: Z(bb)+jets / Multijets). Point into scDisplayNames_rdb, which lives
    // for the whole function, so these const char* are safe in arrays too.
    const char* sigLab = scDisplayNames_rdb[SampleConfig::kZbb].c_str();
    const char* bkgLab = scDisplayNames_rdb[SampleConfig::kMultijets].c_str();
    auto allPaths_rdb = SampleConfig::AllSamplePaths(cfg, std::string(wp_raw));
    const int nSC_rdb = SampleConfig::kNSamples;

    TFile* fmc = TFile::Open(zbbPath.c_str(), "READ");
    TFile* fdata = TFile::Open(dataPath.c_str(), "READ");
    TFile* fmultijets = TFile::Open(mjPath.c_str(), "READ");

    if (!fmc || fmc->IsZombie() || !fdata || fdata->IsZombie() || !fmultijets || fmultijets->IsZombie()) {
        Log::error("FitResponse", "Cannot open input files for WP " + std::string(wp_label));
        if (fmc) {
            fmc->Close();
            delete fmc;
        }
        if (fdata) {
            fdata->Close();
            delete fdata;
        }
        if (fmultijets) {
            fmultijets->Close();
            delete fmultijets;
        }
        return;
    }

    TFile* fextra_rdb[SampleConfig::kNSamples] = {nullptr};
    bool extraLoaded_rdb[SampleConfig::kNSamples] = {false};
    for (int s = 0; s < nSC_rdb; ++s) {
        if (s == SampleConfig::kData || s == SampleConfig::kMultijets || s == SampleConfig::kZbb) continue;
        // Skip padding sentinels silently (from SampleConfig padding when the
        // other channel's minor-bkg flag is on but this one's is off).
        if (scFileNames_rdb[s].rfind("__unused_", 0) == 0) continue;
        fextra_rdb[s] = TFile::Open(allPaths_rdb[s].c_str(), "READ");
        if (!fextra_rdb[s] || fextra_rdb[s]->IsZombie()) {
            if (fextra_rdb[s]) {
                delete fextra_rdb[s];
                fextra_rdb[s] = nullptr;
            }
        } else {
            extraLoaded_rdb[s] = true;
        }
    }

    TTree* tree_mc = Validate::GetTree(fmc, "nominal", "FitResponse");
    TTree* tree_data = Validate::GetTree(fdata, "nominal", "FitResponse");
    TTree* tree_multijets = Validate::GetTree(fmultijets, "nominal", "FitResponse");
    if (!tree_mc || !tree_data || !tree_multijets) {
        Log::error("FitResponse", "Tree 'nominal' not found for WP " + std::string(wp_label));
        fmc->Close();
        delete fmc;
        fdata->Close();
        delete fdata;
        fmultijets->Close();
        delete fmultijets;
        for (int s = 0; s < nSC_rdb; ++s)
            if (fextra_rdb[s]) {
                fextra_rdb[s]->Close();
                delete fextra_rdb[s];
            }
        return;
    }

    TTree* tree_extra_rdb[SampleConfig::kNSamples] = {nullptr};
    for (int s = 0; s < nSC_rdb; ++s) {
        if (extraLoaded_rdb[s]) {
            tree_extra_rdb[s] = (TTree*)fextra_rdb[s]->Get("nominal");
            if (!tree_extra_rdb[s]) extraLoaded_rdb[s] = false;
        }
    }

    // Create histograms
    auto MkHS = [&](const char* name, const char* title, int nb, double lo, double hi) {
        TH1F* h = new TH1F(Form("%s_%s_%s", name, cfg.key.c_str(), wp_label), title, nb, lo, hi);
        h->Sumw2();
        return h;
    };

    TH1F* h_mc_mass = MkHS("h_mc_mass", "ljet_mass", 30, 50, 150);
    TH1F* h_data_mass = MkHS("h_data_mass", "ljet_mass", 30, 50, 150);
    TH1F* h_multijets_mass = MkHS("h_multijets_mass", "ljet_mass", 30, 50, 150);

    // Response-histogram binning: static or dynamic
    // Static (default): cfg.nBinsR and cfg.nBinsTruth.
    // Dynamic (EASYBJETS_RESPONSE_BINNING_MODE=dynamic): fill a 100-bin probe
    // histogram from MC and let BinningConfig::ComputeResponse choose the
    // binning with the EASYBJETS_BINNING algorithm.
    BinningConfig::Result resp_bin;
    resp_bin.nbins = nBins; // static default
    BinningConfig::Result truth_bin;
    truth_bin.nbins = cfg.nBinsTruth; // static default for truth
    bool use_dynamic_resp = BinningConfig::ResponseDynamicEnabled();

    // Truth binning policy (EASYBJETS_BINNING_TRUTH), only used in dynamic mode.
    // The default "match" gives the truth histograms the reco binning.
    BinningConfig::Mode truth_alg = BinningConfig::Mode::Fixed; // unused unless policy=Algorithm
    BinningConfig::TruthPolicy truth_policy = BinningConfig::CurrentTruthPolicy(truth_alg);

    if (use_dynamic_resp) {
        std::cout << "[FitResponse] dynamic response binning enabled "
                  << "(reco_mode=" << BinningConfig::ModeTag()
                  << ", truth_policy=" << BinningConfig::TruthPolicyName(truth_policy);
        if (truth_policy == BinningConfig::TruthPolicy::Algorithm)
            std::cout << "(" << BinningConfig::ModeTag(truth_alg) << ")";
        std::cout << ")" << std::endl;

        // A separate truth probe is only needed when the truth binning is
        // computed from the truth distribution itself (auto or an explicit
        // algorithm); match and fixed do not look at it.
        const bool need_truth_probe = (truth_policy == BinningConfig::TruthPolicy::Auto)
                                      || (truth_policy == BinningConfig::TruthPolicy::Algorithm);

        // Both probes are filled in one pass over the MC tree.
        TH1F* h_probe_R = new TH1F(Form("h_probe_R_%s_%s", cfg.key.c_str(), wp_label),
                                   "R_response probe (reco/ref)", 100, 0.0, 2.0);
        h_probe_R->Sumw2();

        TH1F* h_probe_T = nullptr;
        if (need_truth_probe) {
            h_probe_T = new TH1F(Form("h_probe_T_%s_%s", cfg.key.c_str(), wp_label),
                                 "R_response probe (reco/truth)", 100, 0.0, 2.0);
            h_probe_T->Sumw2();
        }

        bool have_truth_branch =
            FillResponseBinningHists(tree_mc, cfg, need_truth_probe, h_probe_R, h_probe_T, "[FitResponse]");
        // Without the truth branch there is nothing to bin on: fall back to match.
        if (need_truth_probe && !have_truth_branch) truth_policy = BinningConfig::TruthPolicy::Match;

        // Wide bounds on the bin count; the static values are 30-50.
        const int min_r = 10, max_r = 100;
        resp_bin = BinningConfig::ComputeResponse(h_probe_R, cfg.channel, 0.0, 2.0, min_r, max_r);

        switch (truth_policy) {
        case BinningConfig::TruthPolicy::Match:
            truth_bin = resp_bin;
            break;
        case BinningConfig::TruthPolicy::Fixed:
            // nbins was already set to cfg.nBinsTruth above.
            truth_bin.tag = "fixed-truth";
            truth_bin.caption = "binning: fixed-truth, " + std::to_string(truth_bin.nbins) + " bins";
            break;
        case BinningConfig::TruthPolicy::Auto:
            // Same algorithm as reco, run on the reco/truth probe.
            truth_bin = BinningConfig::ComputeResponseWithMode(h_probe_T, cfg.channel, 0.0, 2.0, min_r, max_r,
                                                               BinningConfig::CurrentMode(), "-truth");
            break;
        case BinningConfig::TruthPolicy::Algorithm:
            // Explicit algorithm (e.g. EASYBJETS_BINNING_TRUTH=scott).
            truth_bin = BinningConfig::ComputeResponseWithMode(h_probe_T, cfg.channel, 0.0, 2.0, min_r, max_r,
                                                               truth_alg, "-truth");
            break;
        }

        delete h_probe_R;
        if (h_probe_T) delete h_probe_T;
    }

    auto MkR = [&](const std::string& stem) -> TH1F* {
        // R histogram: BinningConfig::MakeHist in dynamic mode (it handles
        // variable bin edges), otherwise nBins uniform bins on [0, 2].
        std::string fullName = stem + "_" + cfg.key + "_" + std::string(wp_label);
        if (use_dynamic_resp) {
            TH1F* h = BinningConfig::MakeHist(fullName, cfg.branch_R, resp_bin, 0.0, 2.0);
            h->Sumw2();
            return h;
        }
        return MkHS(stem.c_str(), cfg.branch_R.c_str(), nBins, 0.0, 2.0);
    };

    // Truth histogram: like MkR but with its own bin count, because the truth
    // peaks are narrower than the reco one. Static mode uses cfg.nBinsTruth;
    // dynamic mode uses truth_bin, already resolved from the policy above
    // (match = reco binning, auto/algorithm = from the reco/truth probe,
    // fixed = cfg.nBinsTruth). The title is passed explicitly because these
    // histograms are not named after a branch.
    auto MkRTruth = [&](const std::string& stem, const std::string& customTitle) -> TH1F* {
        std::string fullName = stem + "_" + cfg.key + "_" + std::string(wp_label);
        if (use_dynamic_resp) {
            TH1F* h = BinningConfig::MakeHist(fullName, customTitle, truth_bin, 0.0, 2.0);
            h->Sumw2();
            return h;
        }
        return MkHS(stem.c_str(), customTitle.c_str(), cfg.nBinsTruth, 0.0, 2.0);
    };

    // R histograms used in the fits. The MC ones are truth-matched (see the
    // file header); data has no truth requirement.
    TH1F* h_mc_R_balance = MkR("h_mc_" + cfg.branch_R);
    TH1F* h_data_R_balance = MkR("h_data_" + cfg.branch_R);
    TH1F* h_multijets_R_balance = MkR("h_multijets_" + cfg.branch_R);

    // Same without the truth requirement, for the data/MC stack plot.
    TH1F* h_mc_R_balance_all = MkR("h_mc_" + cfg.branch_R + "_all");
    TH1F* h_multijets_R_balance_all = MkR("h_multijets_" + cfg.branch_R + "_all");

    // Truth/ref and reco/truth for signal, main background and their sum
    // (cfg.nBinsTruth bins in static mode: 30 for DB, 90 for MJB).
    TH1F* h_mc_truth_R_balance = MkRTruth("h_mc_truth_" + cfg.branch_R, "p_{T}^{truth}/p_{T}^{ref,db}");
    TH1F* h_mc_reco_truth = MkRTruth("h_mc_reco_truth", "p_{T}^{j1}/p_{T}^{truth}");
    TH1F* h_bkg_truth_R_balance = MkRTruth("h_bkg_truth_" + cfg.branch_R, "Bkg truth/ref");
    TH1F* h_bkg_reco_truth = MkRTruth("h_bkg_reco_truth", "Bkg reco/truth");
    TH1F* h_comb_truth_R_balance = MkRTruth("h_comb_truth_" + cfg.branch_R, "Comb truth/ref");
    TH1F* h_comb_reco_truth = MkRTruth("h_comb_reco_truth", "Comb reco/truth");
    h_mc_truth_R_balance->SetLineColor(kAzure + 2);
    h_mc_reco_truth->SetLineColor(kOrange + 1);

    // Kinematic control histograms: signal, main background, data, sum of the two MC
    static const int kNS = 4;
    const char* sNames[kNS] = {"sig", "bkg", "data", "comb"};
    TH1F *h_jet_pt[kNS], *h_ref_pt[kNS], *h_dR_samp[kNS];
    for (int s = 0; s < kNS; ++s) {
        h_jet_pt[s] = MkHS(Form("h_jet_pt_%s", sNames[s]), ";p_{T}^{jet} [GeV];Events", 50, 0, 800);
        h_ref_pt[s] = MkHS(Form("h_ref_pt_%s", sNames[s]), ";p_{T}^{ref} [GeV];Events", 50, 0, 800);
        h_dR_samp[s] = MkHS(Form("h_dR_%s", sNames[s]), ";#DeltaR(reco, truth);Events", 50, 0, 1.2);
    }

    // Colors
    h_data_mass->SetLineColor(kBlack);
    h_data_R_balance->SetLineColor(kBlack);
    h_multijets_mass->SetLineColor(kBlue);
    h_multijets_R_balance->SetLineColor(kBlue);
    h_mc_mass->SetLineColor(kRed);
    h_mc_R_balance->SetLineColor(kRed);
    h_data_mass->SetFillColor(kBlack);
    h_data_R_balance->SetFillColor(kBlack);
    h_multijets_mass->SetFillColor(kBlue);
    h_multijets_R_balance->SetFillColor(kBlue);
    h_mc_mass->SetFillColor(kRed);
    h_mc_R_balance->SetFillColor(kRed);
    h_mc_R_balance_all->SetLineColor(kRed);
    h_mc_R_balance_all->SetFillColor(kRed);
    h_multijets_R_balance_all->SetLineColor(kBlue);
    h_multijets_R_balance_all->SetFillColor(kBlue);

    TH1F* h_extra_mass[SampleConfig::kNSamples];
    TH1F* h_extra_R_balance[SampleConfig::kNSamples];
    for (int s = 0; s < nSC_rdb; ++s) {
        h_extra_mass[s] = MkHS(Form("h_%s_mass", scFileNames_rdb[s].c_str()), "ljet_mass", 30, 50, 150);
        h_extra_R_balance[s] = MkR(Form(("h_%s_" + cfg.branch_R).c_str(), scFileNames_rdb[s].c_str()));
        h_extra_mass[s]->SetFillColor(scColors_rdb[s]);
        h_extra_mass[s]->SetLineColor(scColors_rdb[s]);
        h_extra_R_balance[s]->SetFillColor(scColors_rdb[s]);
        h_extra_R_balance[s]->SetLineColor(scColors_rdb[s]);
    }

    // Branch setup
    std::vector<float>* ljet_m = nullptr;
    std::vector<float>* ljet_pt_vec = nullptr;
    double weight, weight_multijets, weight_extra_rdb;
    double R_response;
    int selected_ljet_index = -1;

    tree_mc->SetBranchAddress(cfg.branch_mass.c_str(), &ljet_m);
    tree_multijets->SetBranchAddress(cfg.branch_mass.c_str(), &ljet_m);
    tree_mc->SetBranchAddress(cfg.branch_pt.c_str(), &ljet_pt_vec);
    tree_multijets->SetBranchAddress(cfg.branch_pt.c_str(), &ljet_pt_vec);

    // Variant-branch consistency check on data
    // Data and MC must be read from the same branches, otherwise every plot of
    // the variant compares nominal data with corrected MC. The slimmer writes
    // the corrected branches for data too, so a missing branch means a broken
    // ntuple: stop here instead of falling back to the nominal branches.
    {
        const bool dataHasVariantMass = Validate::HasBranch(tree_data, cfg.branch_mass);
        const bool dataHasVariantPt = Validate::HasBranch(tree_data, cfg.branch_pt);
        const bool dataHasRDB = Validate::HasBranch(tree_data, cfg.branch_R);

        if (!dataHasVariantMass || !dataHasVariantPt || !dataHasRDB) {
            Log::error("FitResponse",
                       Form("Variant '%s' branches missing on DATA tree for WP %s:\n"
                            "    %s : %s\n"
                            "    %s : %s\n"
                            "    %s : %s\n"
                            "  Data file: %s\n"
                            "  Refusing to fall back to nominal branches — that would compare\n"
                            "  nominal-mass data against %s-mass MC. Re-run the slimmer so the\n"
                            "  data ntuple carries the variant branches, or remove this variant\n"
                            "  from GetVariantConfigs().",
                            cfg.key.c_str(), wp_label, cfg.branch_mass.c_str(),
                            dataHasVariantMass ? "OK" : "MISSING", cfg.branch_pt.c_str(),
                            dataHasVariantPt ? "OK" : "MISSING", cfg.branch_R.c_str(),
                            dataHasRDB ? "OK" : "MISSING", dataPath.c_str(), cfg.path_tag.c_str()));

            fmc->Close();
            delete fmc;
            fdata->Close();
            delete fdata;
            fmultijets->Close();
            delete fmultijets;
            for (int s = 0; s < nSC_rdb; ++s)
                if (fextra_rdb[s]) {
                    fextra_rdb[s]->Close();
                    delete fextra_rdb[s];
                }
            return;
        }

        tree_data->SetBranchAddress(cfg.branch_mass.c_str(), &ljet_m);
        tree_data->SetBranchAddress(cfg.branch_pt.c_str(), &ljet_pt_vec);
    }

    tree_mc->SetBranchAddress("selected_ljet_index", &selected_ljet_index);
    tree_data->SetBranchAddress("selected_ljet_index", &selected_ljet_index);
    tree_multijets->SetBranchAddress("selected_ljet_index", &selected_ljet_index);
    tree_mc->SetBranchAddress(cfg.branch_R.c_str(), &R_response);
    tree_data->SetBranchAddress(cfg.branch_R.c_str(), &R_response);
    tree_multijets->SetBranchAddress(cfg.branch_R.c_str(), &R_response);
    tree_mc->SetBranchAddress("weight", &weight);
    tree_multijets->SetBranchAddress("weight", &weight_multijets);

    double R_truth_response = -1.0, R_truth_response_bkg = -1.0;
    bool hasTruthBranches = false, hasTruthBkg = false;
    std::vector<float>* matched_tjet_dR_sig = nullptr;
    std::vector<float>* matched_tjet_dR_bkg = nullptr;

    if (Validate::HasBranch(tree_mc, ("R_truth_" + cfg.channel).c_str())) {
        tree_mc->SetBranchAddress(("R_truth_" + cfg.channel).c_str(), &R_truth_response);
        hasTruthBranches = true;
    }
    if (Validate::HasBranch(tree_mc, "matched_tjet_dR"))
        tree_mc->SetBranchAddress("matched_tjet_dR", &matched_tjet_dR_sig);
    if (Validate::HasBranch(tree_multijets, ("R_truth_" + cfg.channel).c_str())) {
        tree_multijets->SetBranchAddress(("R_truth_" + cfg.channel).c_str(), &R_truth_response_bkg);
        hasTruthBkg = true;
    }
    if (Validate::HasBranch(tree_multijets, "matched_tjet_dR"))
        tree_multijets->SetBranchAddress("matched_tjet_dR", &matched_tjet_dR_bkg);

    for (int s = 0; s < nSC_rdb; ++s) {
        if (tree_extra_rdb[s]) {
            tree_extra_rdb[s]->SetBranchAddress(cfg.branch_mass.c_str(), &ljet_m);
            tree_extra_rdb[s]->SetBranchAddress("selected_ljet_index", &selected_ljet_index);
            tree_extra_rdb[s]->SetBranchAddress(cfg.branch_R.c_str(), &R_response);
            tree_extra_rdb[s]->SetBranchAddress("weight", &weight_extra_rdb);
        }
    }

    // Fill loops (inclusive in pT; mass window only with EASYBJETS_MASS_CUT)
    // Truth-match acceptance counters, for the printout only (the fits use
    // Sumw2 histograms). Raw counts alone are misleading for weighted MC, so
    // the weighted yield and N_eff = (sum w)^2 / sum w^2 are kept as well.
    Long64_t n_mc_total = 0, n_mc_truth_pass = 0;
    Long64_t n_bkg_total = 0, n_bkg_truth_pass = 0;
    double sumw_mc_total = 0, sumw2_mc_total = 0;
    double sumw_mc_pass = 0, sumw2_mc_pass = 0;
    double sumw_bkg_total = 0, sumw2_bkg_total = 0;
    double sumw_bkg_pass = 0, sumw2_bkg_pass = 0;

    // MC signal
    for (Long64_t i = 0; i < tree_mc->GetEntries(); i++) {
        tree_mc->GetEntry(i);
        if (selected_ljet_index < 0 || !ljet_m || selected_ljet_index >= (int)ljet_m->size()) continue;
        if (MassCutConfig::Reject(ljet_m, selected_ljet_index))
            continue; // optional mass window (EASYBJETS_MASS_CUT)
        double jetMass = (*ljet_m)[selected_ljet_index];
        h_mc_mass->Fill(jetMass, weight);

        h_mc_R_balance_all->Fill(R_response, weight);
        ++n_mc_total;
        sumw_mc_total += weight;
        sumw2_mc_total += weight * weight;

        // dR(reco, truth) of the selected jet; -1 if there is no match
        double dR_truth = -1.0;
        if (matched_tjet_dR_sig && selected_ljet_index < (int)matched_tjet_dR_sig->size())
            dR_truth = matched_tjet_dR_sig->at(selected_ljet_index);

        bool truthOK = hasTruthBranches && R_truth_response > 0 && dR_truth >= 0 && dR_truth < cfg.maxTruthDR;

        if (truthOK) {
            // The three closure histograms get exactly the same events
            h_mc_R_balance->Fill(R_response, weight);
            h_mc_truth_R_balance->Fill(R_truth_response, weight);
            h_mc_reco_truth->Fill(R_response / R_truth_response, weight);
            ++n_mc_truth_pass;
            sumw_mc_pass += weight;
            sumw2_mc_pass += weight * weight;
        }

        if (ljet_pt_vec && selected_ljet_index < (int)ljet_pt_vec->size()) {
            double jpt = (*ljet_pt_vec)[selected_ljet_index];
            h_jet_pt[0]->Fill(jpt, weight);
            if (R_response > 0) h_ref_pt[0]->Fill(jpt / R_response, weight);
        }
        if (dR_truth >= 0) h_dR_samp[0]->Fill(dR_truth, weight);
    }

    // Data
    for (Long64_t i = 0; i < tree_data->GetEntries(); i++) {
        tree_data->GetEntry(i);
        if (selected_ljet_index < 0 || !ljet_m || selected_ljet_index >= (int)ljet_m->size()) continue;
        if (MassCutConfig::Reject(ljet_m, selected_ljet_index))
            continue; // optional mass window (EASYBJETS_MASS_CUT)
        double jetMass = (*ljet_m)[selected_ljet_index];
        h_data_mass->Fill(jetMass);
        h_data_R_balance->Fill(R_response);
        if (ljet_pt_vec && selected_ljet_index < (int)ljet_pt_vec->size()) {
            double jpt = (*ljet_pt_vec)[selected_ljet_index];
            h_jet_pt[2]->Fill(jpt);
            if (R_response > 0) h_ref_pt[2]->Fill(jpt / R_response);
        }
    }

    // Main background, first pass: mass only. R needs SF(m), which is fitted
    // on this mass histogram, so it is filled in a second pass below.
    for (Long64_t i = 0; i < tree_multijets->GetEntries(); i++) {
        tree_multijets->GetEntry(i);
        if (selected_ljet_index < 0 || !ljet_m || selected_ljet_index >= (int)ljet_m->size()) continue;
        if (MassCutConfig::Reject(ljet_m, selected_ljet_index))
            continue; // optional mass window (EASYBJETS_MASS_CUT)
        h_multijets_mass->Fill((*ljet_m)[selected_ljet_index], weight_multijets);
    }

    // Minor backgrounds: no SF and no truth requirement
    for (int s = 0; s < nSC_rdb; ++s) {
        if (!tree_extra_rdb[s]) continue;
        for (Long64_t i = 0; i < tree_extra_rdb[s]->GetEntries(); i++) {
            tree_extra_rdb[s]->GetEntry(i);
            if (selected_ljet_index < 0 || !ljet_m || selected_ljet_index >= (int)ljet_m->size()) continue;
            if (MassCutConfig::Reject(ljet_m, selected_ljet_index))
                continue; // optional mass window (EASYBJETS_MASS_CUT)
            double jetMass = (*ljet_m)[selected_ljet_index];
            h_extra_mass[s]->Fill(jetMass, weight_extra_rdb);
            h_extra_R_balance[s]->Fill(R_response, weight_extra_rdb);
        }
    }

    // Background scale factor SF(m), one per variant and working point
    SFUtils::LinearSF sf =
        GetOrFitSF(h_data_mass, h_multijets_mass, cfg, std::string(wp_raw), basePath.Data());
    Log::info("FitResponse",
              Form("SF(m) = %.4f (+/-%.4f) + %.8f (+/-%.8f) * m", sf.p0, sf.p0_err, sf.p1, sf.p1_err));

    // Scaled background mass histogram, for the stack plot
    TH1F* h_multijets_mass_scaled =
        (TH1F*)h_multijets_mass->Clone(Form("h_multijets_mass_scaled_%s_%s", cfg.key.c_str(), wp_label));
    SFUtils::ApplySFBinByBin(h_multijets_mass_scaled, sf);

    // Main background, second pass: R weighted by SF(m) event by event
    for (Long64_t i = 0; i < tree_multijets->GetEntries(); ++i) {
        tree_multijets->GetEntry(i);
        if (selected_ljet_index < 0 || !ljet_m || selected_ljet_index >= (int)ljet_m->size()) continue;
        if (MassCutConfig::Reject(ljet_m, selected_ljet_index))
            continue; // optional mass window (EASYBJETS_MASS_CUT)
        double jetMass = (*ljet_m)[selected_ljet_index];
        double w_corr = SFUtils::GetEventWeight(weight_multijets, jetMass, sf);

        h_multijets_R_balance_all->Fill(R_response, w_corr);
        ++n_bkg_total;
        sumw_bkg_total += w_corr;
        sumw2_bkg_total += w_corr * w_corr;

        double dR_truth = -1.0;
        if (matched_tjet_dR_bkg && selected_ljet_index < (int)matched_tjet_dR_bkg->size())
            dR_truth = matched_tjet_dR_bkg->at(selected_ljet_index);

        bool truthOK = hasTruthBkg && R_truth_response_bkg > 0 && dR_truth >= 0 && dR_truth < cfg.maxTruthDR;

        if (truthOK) {
            h_multijets_R_balance->Fill(R_response, w_corr);
            h_bkg_truth_R_balance->Fill(R_truth_response_bkg, w_corr);
            h_bkg_reco_truth->Fill(R_response / R_truth_response_bkg, w_corr);
            ++n_bkg_truth_pass;
            sumw_bkg_pass += w_corr;
            sumw2_bkg_pass += w_corr * w_corr;
        }

        if (ljet_pt_vec && selected_ljet_index < (int)ljet_pt_vec->size()) {
            double jpt = (*ljet_pt_vec)[selected_ljet_index];
            h_jet_pt[1]->Fill(jpt, w_corr);
            if (R_response > 0) h_ref_pt[1]->Fill(jpt / R_response, w_corr);
        }
        if (dR_truth >= 0) h_dR_samp[1]->Fill(dR_truth, w_corr);
    }
    Log::info("FitResponse", "Pass 2: gamma+jets R_response filled with per-event SF");

    // Diagnostic: truth-match acceptance
    // The weighted fraction is the one to compare with the cutflow tables;
    // N_eff is what sets the fit uncertainties.
    auto Neff = [](double sw, double sw2) {
        return (sw2 > 0) ? (sw * sw) / sw2 : 0.0;
    };
    Log::info("FitResponse", Form("Truth-match selection (DR < %.2f):", cfg.maxTruthDR));
    if (n_mc_total > 0) {
        Log::info("FitResponse",
                  Form("  MC signal: raw %lld/%lld (%.1f%%) | "
                       "Sum_w %.2f/%.2f (%.1f%%) | N_eff %.0f/%.0f",
                       n_mc_truth_pass, n_mc_total, 100.0 * n_mc_truth_pass / n_mc_total, sumw_mc_pass,
                       sumw_mc_total, sumw_mc_total > 0 ? 100.0 * sumw_mc_pass / sumw_mc_total : 0.0,
                       Neff(sumw_mc_pass, sumw2_mc_pass), Neff(sumw_mc_total, sumw2_mc_total)));
    }
    if (n_bkg_total > 0) {
        Log::info("FitResponse",
                  Form("  gamma+jets: raw %lld/%lld (%.1f%%) | "
                       "Sum_w %.2f/%.2f (%.1f%%) | N_eff %.0f/%.0f",
                       n_bkg_truth_pass, n_bkg_total, 100.0 * n_bkg_truth_pass / n_bkg_total, sumw_bkg_pass,
                       sumw_bkg_total, sumw_bkg_total > 0 ? 100.0 * sumw_bkg_pass / sumw_bkg_total : 0.0,
                       Neff(sumw_bkg_pass, sumw2_bkg_pass), Neff(sumw_bkg_total, sumw2_bkg_total)));
    }

    // Mass stack plot
    {
        CanvasGuard c(Form("mass_%s_%s", cfg.key.c_str(), wp_label), "Jet Mass", 800, 600);
        THStack* stk = new THStack("mc_stack", "Large-R Jet Mass;M [GeV];Events");
        for (int s = 0; s < nSC_rdb; ++s) {
            if (s == SampleConfig::kData || s == SampleConfig::kMultijets || s == SampleConfig::kZbb)
                continue;
            if (extraLoaded_rdb[s]) {
                h_extra_mass[s]->SetFillStyle(1001);
                h_extra_mass[s]->SetBit(kCanDelete, kFALSE);
                stk->Add(h_extra_mass[s]);
            }
        }
        h_multijets_mass_scaled->SetFillColor(kBlue);
        h_multijets_mass_scaled->SetLineColor(kBlue);
        h_multijets_mass_scaled->SetBit(kCanDelete, kFALSE);
        stk->Add(h_multijets_mass_scaled);
        h_mc_mass->SetFillColor(kRed);
        h_mc_mass->SetLineColor(kRed);
        h_mc_mass->SetBit(kCanDelete, kFALSE);
        stk->Add(h_mc_mass);
        h_data_mass->Draw("E");
        stk->Draw("HIST SAME");
        h_data_mass->Draw("E SAME");
        TLegend* leg = new TLegend(0.55, 0.55, 0.89, 0.89);
        leg->SetBorderSize(0);
        leg->SetTextSize(0.025);
        leg->AddEntry(h_data_mass, "Data", "lp");
        leg->AddEntry(h_mc_mass, sigLab, "f");
        leg->AddEntry(h_multijets_mass_scaled, bkgLab, "f");
        for (int s = 0; s < nSC_rdb; ++s) {
            if (s == SampleConfig::kData || s == SampleConfig::kMultijets || s == SampleConfig::kZbb)
                continue;
            if (extraLoaded_rdb[s]) leg->AddEntry(h_extra_mass[s], scDisplayNames_rdb[s].c_str(), "f");
        }
        leg->Draw();
        AtlasUtils::DrawAtlasLabel(0.18, 0.85);
        c->SaveAs(Form("%s/jet_mass_comparison.pdf", basePath.Data()));
    }

    // Combined MC = signal + main background + minor backgrounds
    // Truth-matched sum, used in the fits
    TH1F* h_combined_R_balance =
        (TH1F*)h_mc_R_balance->Clone(Form("h_combined_R_balance_%s_%s", cfg.key.c_str(), wp_label));
    h_combined_R_balance->Add(h_multijets_R_balance);
    for (int s = 0; s < nSC_rdb; ++s) {
        if (s == SampleConfig::kData || s == SampleConfig::kMultijets || s == SampleConfig::kZbb) continue;
        if (extraLoaded_rdb[s]) h_combined_R_balance->Add(h_extra_R_balance[s]);
    }
    h_combined_R_balance->SetFillColor(kGreen - 2);
    h_combined_R_balance->SetLineColor(kGreen - 2);

    // Sum without truth requirement, for the data/MC ratio of the stack plot
    TH1F* h_combined_R_balance_all =
        (TH1F*)h_mc_R_balance_all->Clone(Form("h_combined_R_balance_all_%s_%s", cfg.key.c_str(), wp_label));
    h_combined_R_balance_all->Add(h_multijets_R_balance_all);
    for (int s = 0; s < nSC_rdb; ++s) {
        if (s == SampleConfig::kData || s == SampleConfig::kMultijets || s == SampleConfig::kZbb) continue;
        if (extraLoaded_rdb[s]) h_combined_R_balance_all->Add(h_extra_R_balance[s]);
    }
    h_combined_R_balance_all->SetFillColor(kGreen - 2);
    h_combined_R_balance_all->SetLineColor(kGreen - 2);

    // The truth sums only contain signal and main background (the minor
    // backgrounds have no truth histograms for now - produce new 
    // ntuples that include the truth if needed).
    h_comb_truth_R_balance->Add(h_mc_truth_R_balance);
    if (hasTruthBkg) h_comb_truth_R_balance->Add(h_bkg_truth_R_balance);
    h_comb_reco_truth->Add(h_mc_reco_truth);
    if (hasTruthBkg) h_comb_reco_truth->Add(h_bkg_reco_truth);
    h_jet_pt[3]->Add(h_jet_pt[0]);
    h_jet_pt[3]->Add(h_jet_pt[1]);
    h_ref_pt[3]->Add(h_ref_pt[0]);
    h_ref_pt[3]->Add(h_ref_pt[1]);
    h_dR_samp[3]->Add(h_dR_samp[0]);
    h_dR_samp[3]->Add(h_dR_samp[1]);

    // Fit of one R histogram
    // Fills mean, sigma and chi2/ndf into `out` and returns the full result
    // (amplitude, final window) for the callers that need it. Without an
    // explicit seed window, the window is derived from the histogram itself.
    auto DoFit = [&](TH1F* h, const char* tag_prefix, const char* label, RDBFitResults& out,
                     double seedLo = -1, double seedHi = -1) -> FitUtils::GausFitResult {
        // Data-driven seed window: mu +/- 2 RMS, with mean and RMS recomputed
        // once inside mu +/- 3 RMS so the long tails do not inflate them. A wider
        // distribution gets a wider seed without any per-variant tuning; +/-2
        // sigma holds about 95% of a Gaussian, enough for the scan to converge inward.
        if (seedLo < 0 || seedHi < 0) {
            double m1 = h->GetMean();
            double s1 = h->GetRMS();
            if (s1 > 0) {
                int saveFirst = h->GetXaxis()->GetFirst();
                int saveLast = h->GetXaxis()->GetLast();
                h->GetXaxis()->SetRangeUser(m1 - 3 * s1, m1 + 3 * s1);
                double m2 = h->GetMean();
                double s2 = h->GetRMS();
                h->GetXaxis()->SetRange(saveFirst, saveLast); // later fits must see the full range
                if (s2 > 0) {
                    m1 = m2;
                    s1 = s2;
                }
            }
            const double seed_n_sigma = 2.;
            double dataLo = m1 - seed_n_sigma * s1;
            double dataHi = m1 + seed_n_sigma * s1;
            dataLo = std::max(dataLo, h->GetXaxis()->GetXmin());
            dataHi = std::min(dataHi, h->GetXaxis()->GetXmax());
            if (seedLo < 0) seedLo = dataLo;
            if (seedHi < 0) seedHi = dataHi;
            Log::info("FitResponse",
                      Form("Data-driven seed [%s]: mu=%.4f rms_trunc=%.4f -> window [%.3f, %.3f]", tag_prefix,
                           m1, s1, seedLo, seedHi));
        }
        const char* tag = Form("%s_%s_%s", tag_prefix, cfg.key.c_str(), wp_label);

        // Raw mean and RMS, to compare with the fitted sigma: an RMS much
        // smaller than sigma means the fit locked onto the wrong window.
        // N_eff is printed next to the raw entries because it, not the raw
        // count, sets the fit error for weighted MC.
        Log::info(
            "FitResponse",
            Form("Pre-fit [%s]: raw=%lld N_eff=%.0f integral=%.2f mean=%.4f rms=%.4f peakBin=%.1f nbins=%d",
                 tag_prefix, (Long64_t)h->GetEntries(), h->GetEffectiveEntries(), h->Integral(), h->GetMean(),
                 h->GetRMS(), h->GetBinContent(h->GetMaximumBin()), h->GetNbinsX()));

        auto fit = FitUtils::IterativeGausFit(
            h, tag, Form("%s  %s", label, wp_label),
            Form("%s/Fit_R_balance_%s_scan.pdf", basePath.Data(), tag_prefix), seedLo, seedHi,
            /*maxNSigma*/ 1.5,
            /*physicalLo*/ []() {
                // Lowest R the fit window may reach. Far below the peak the
                // histogram is shaped by the jet pT threshold, not by the
                // resolution; 0.3 is well below the core (sigma ~ 0.1) and
                // keeps the scan out of that region. A negative value
                // disables the floor.
                const char* e = std::getenv("EASYBJETS_RESPONSE_FIT_LO");
                return (e && *e) ? std::atof(e) : 0.3;
            }());
        out.mean = fit.mean;
        out.meanErr = fit.meanErr;
        out.sigma = fit.sigma;
        out.sigmaErr = fit.sigmaErr;
        out.chi2ndf = fit.chi2ndf;

        CanvasGuard c(Form("c_%s", tag), label, 800, 600);
        h->Draw(std::string(tag_prefix) == "data" ? "EP" : "HIST");

        // gf is the fitted Gaussian extended to +/-4 sigma; it is set up but
        // not drawn. Only the curve over the actual fit window is shown.
        double drawLo = std::max(h->GetXaxis()->GetXmin(), fit.mean - 4.0 * fit.sigma);
        double drawHi = std::min(h->GetXaxis()->GetXmax(), fit.mean + 4.0 * fit.sigma);
        TF1* gf = new TF1(Form("gf_%s", tag), "gaus", drawLo, drawHi);
        gf->SetParameters(fit.amp, fit.mean, fit.sigma);
        gf->SetLineColor(kViolet - 4);
        gf->SetLineStyle(2);
        gf->SetLineWidth(1);
        // Fitted curve over the fit window, drawn with the fitted model
        // (Gaussian by default, or the EASYBJETS_FIT_MODEL shape).
        TF1* gfFitRange = FitUtils::MakeModelTF1(fit, Form("gfR_%s", tag));
        gfFitRange->SetLineColor(kViolet - 4);
        gfFitRange->SetLineStyle(1);
        gfFitRange->SetLineWidth(3);
        gfFitRange->Draw("same");

        TLegend* leg = new TLegend(0.6, 0.70, 0.9, 0.9);
        leg->SetBorderSize(0);
        leg->SetTextSize(0.025);
        leg->AddEntry(h, label, std::string(tag_prefix) == "data" ? "lep" : "l");
        leg->AddEntry(gfFitRange, FitUtils::ModelDisplayName(fit.modelName).c_str(), "l");
        TLatex txt;
        txt.SetTextSize(0.03);
        // Shown and saved errors are the scaled ones; they equal the raw
        // errors unless EASYBJETS_FIT_ERR_SCALE=pdg and chi2/ndf > 1.
        txt.DrawLatexNDC(0.65, 0.65, Form("#mu = %.3f #pm %.3f", fit.mean, fit.meanErrScaled()));
        txt.DrawLatexNDC(0.65, 0.60, Form("#sigma = %.3f #pm %.3f", fit.sigma, fit.sigmaErrScaled()));
        txt.DrawLatexNDC(0.65, 0.55, Form("#chi^{2}/NDF = %.3f", fit.chi2ndf));
        txt.DrawLatexNDC(0.65, 0.50, Form("window: [%.3f, %.3f]", fit.fitLo, fit.fitHi));
        leg->Draw();
        AtlasUtils::DrawAtlasLabel(0.18, 0.85);
        c->SaveAs(Form("%s/Fit_R_balance_%s.pdf", basePath.Data(), tag_prefix));

        FitOutputUtils::SaveFitParameters(Form("%s/%s_fit_parameters.txt", basePath.Data(), tag_prefix),
                                          Form("%s (WP: %s)", label, wp_label), fit.mean, fit.meanErrScaled(),
                                          fit.sigma, fit.sigmaErrScaled(), fit.amp, fit.ampErrScaled(),
                                          fit.chi2, fit.ndf);
        return fit;
    };

    // Fixed-window fit for the truth distributions (currently not called)
    // The truth/ref distribution has a long low tail from imperfect matching,
    // and with it the chi2/ndf window choice of the scan can jump between
    // windows when the binning changes, moving sigma a lot. This alternative
    // has no scan: peak at the median, width from the RMS recomputed in
    // +/-2 RMS around it, then one Gaussian fit in +/-1 sigma.
    auto FitTruthCore = [&](TH1F* h, const char* tag_prefix, const char* label,
                            RDBFitResults& out) -> FitUtils::GausFitResult {
        // The median is robust both to the asymmetric tail (which pulls the
        // mean) and to a flat top (which makes the maximum bin jump).
        double mu;
        {
            const int nq = 1;
            double xq[nq] = {0.5};
            double yq[nq];
            h->GetQuantiles(nq, yq, xq);
            mu = yq[0];
        }
        double sigma = h->GetRMS();
        // Recompute the RMS in mu +/- 2 RMS; mu stays at the median.
        if (sigma > 0) {
            int save_first = h->GetXaxis()->GetFirst();
            int save_last = h->GetXaxis()->GetLast();
            h->GetXaxis()->SetRangeUser(mu - 2 * sigma, mu + 2 * sigma);
            double sigma2 = h->GetRMS();
            h->GetXaxis()->SetRange(save_first, save_last);
            if (sigma2 > 0) sigma = sigma2;
        }

        // Core only: [mu - 1 sigma, mu + 1 sigma]
        double lo = std::max(mu - 1.0 * sigma, h->GetXaxis()->GetXmin());
        double hi = std::min(mu + 1.0 * sigma, h->GetXaxis()->GetXmax());

        const char* tag = Form("%s_%s_%s", tag_prefix, cfg.key.c_str(), wp_label);

        Log::info("FitResponse",
                  Form("FitTruthCore [%s]: mu_trunc=%.4f rms_trunc=%.4f -> window [%.3f, %.3f] nbins=%d",
                       tag_prefix, mu, sigma, lo, hi, h->GetNbinsX()));

        TF1* gf = new TF1(Form("gf_%s", tag), "gaus", lo, hi);
        gf->SetParameters(h->GetMaximum(), mu, sigma);
        h->Fit(gf, "RQN0");

        FitUtils::GausFitResult r;
        r.mean = gf->GetParameter(1);
        r.meanErr = gf->GetParError(1);
        r.sigma = gf->GetParameter(2);
        r.sigmaErr = gf->GetParError(2);
        r.chi2 = gf->GetChisquare();
        r.ndf = gf->GetNDF();
        r.chi2ndf = (r.ndf > 0) ? (r.chi2 / r.ndf) : 0.0;
        r.fitLo = lo;
        r.fitHi = hi;
        r.amp = gf->GetParameter(0);
        r.ampErr = gf->GetParError(0);

        out.mean = r.mean;
        out.meanErr = r.meanErr;
        out.sigma = r.sigma;
        out.sigmaErr = r.sigmaErr;
        out.chi2ndf = r.chi2ndf;

        // Same plot and parameter file as DoFit
        {
            CanvasGuard c(Form("c_%s", tag), label, 800, 600);
            h->Draw("HIST");
            gf->SetLineColor(kViolet - 4);
            gf->SetLineWidth(2);
            gf->Draw("same");
            TLegend* leg = new TLegend(0.6, 0.75, 0.9, 0.9);
            leg->SetBorderSize(0);
            leg->SetTextSize(0.025);
            leg->AddEntry(h, label, "l");
            leg->AddEntry(gf, "Gaussian Fit (core #pm 1#sigma)", "l");
            TLatex txt;
            txt.SetTextSize(0.03);
            txt.DrawLatexNDC(0.65, 0.7, Form("#mu = %.3f #pm %.3f", r.mean, r.meanErr));
            txt.DrawLatexNDC(0.65, 0.65, Form("#sigma = %.3f #pm %.3f", r.sigma, r.sigmaErr));
            txt.DrawLatexNDC(0.65, 0.6, Form("#chi^{2}/NDF = %.3f", r.chi2ndf));
            leg->Draw();
            AtlasUtils::DrawAtlasLabel(0.18, 0.85);
            c->SaveAs(Form("%s/Fit_R_balance_%s.pdf", basePath.Data(), tag_prefix));
        }

        FitOutputUtils::SaveFitParameters(Form("%s/%s_fit_parameters.txt", basePath.Data(), tag_prefix),
                                          Form("%s (WP: %s)", label, wp_label), r.mean, r.meanErr, r.sigma,
                                          r.sigmaErr, r.amp, r.ampErr, r.chi2, r.ndf);

        delete gf;
        return r;
    };


    RDBFitResults mc_fit, bkg_fit, data_fit, comb_fit;
    FitUtils::JERRatio jer;
    FitUtils::JESRatio jes;

    // Data first, with the data-driven seed.
    auto fit_data = DoFit(h_data_R_balance, "data", "Data", data_fit);

    // The MC fits are seeded with the final data window so that sigma_data
    // and sigma_MC, which enter sigma_insitu and sigma^MC_insitu, are measured
    // on comparable windows. This is only the seed: the scan can still move
    // the MC window (and EASYBJETS_SEED_MODE=raw ignores it).
    auto fit_mc = DoFit(h_mc_R_balance, "mc", sigLab, mc_fit, fit_data.fitLo, fit_data.fitHi);
    auto fit_bkg = DoFit(h_multijets_R_balance, "bkg", bkgLab, bkg_fit, fit_data.fitLo, fit_data.fitHi);
    auto fit_comb =
        DoFit(h_combined_R_balance, "combined", "Combined MC", comb_fit, fit_data.fitLo, fit_data.fitHi);

    // b-JER from the combined MC (needs the truth branches)
    if (hasTruthBranches) {
        // The three closure inputs should have the same entries and integral
        Log::info("FitResponse",
                  Form("JER inputs (entries):  R_response^MC=%lld  R_truth=%lld  reco/truth=%lld",
                       (Long64_t)h_combined_R_balance->GetEntries(),
                       (Long64_t)h_comb_truth_R_balance->GetEntries(),
                       (Long64_t)h_comb_reco_truth->GetEntries()));
        Log::info("FitResponse",
                  Form("JER inputs (integral): R_response^MC=%.1f  R_truth=%.1f  reco/truth=%.1f",
                       h_combined_R_balance->Integral(), h_comb_truth_R_balance->Integral(),
                       h_comb_reco_truth->Integral()));

        // Unfitted truth distributions. Drawn from clones so the title, range
        // and colour changes do not reach the histograms fitted below.
        {
            TH1F* h_raw_truth = (TH1F*)h_comb_truth_R_balance->Clone(
                Form("h_comb_truth_R_balance_raw_%s_%s", cfg.key.c_str(), wp_label));
            TH1F* h_raw_reco = (TH1F*)h_comb_reco_truth->Clone(
                Form("h_comb_reco_truth_raw_%s_%s", cfg.key.c_str(), wp_label));
            h_raw_truth->SetDirectory(nullptr);
            h_raw_reco->SetDirectory(nullptr);

            CanvasGuard c(Form("c_raw_truth_%s_%s", cfg.key.c_str(), wp_label), "raw truth distributions",
                          1200, 500);
            c->Divide(2, 1);

            c->cd(1);
            gPad->SetLeftMargin(0.15);
            h_raw_truth->SetLineColor(kAzure + 2);
            h_raw_truth->SetLineWidth(2);
            h_raw_truth->GetXaxis()->SetRangeUser(0.0, 2.0);
            h_raw_truth->GetXaxis()->SetTitle("p_{T}^{truth}/p_{T}^{ref,db}");
            h_raw_truth->GetYaxis()->SetTitle("Events");
            h_raw_truth->SetTitle("Raw truth/ref (combined MC)");
            h_raw_truth->Draw("HIST E");
            TLatex lt;
            lt.SetNDC();
            lt.SetTextSize(0.035);
            lt.DrawLatex(0.18, 0.85, Form("mean = %.4f", h_raw_truth->GetMean()));
            lt.DrawLatex(0.18, 0.80, Form("rms  = %.4f", h_raw_truth->GetRMS()));
            lt.DrawLatex(0.18, 0.75, Form("nbins = %d", h_raw_truth->GetNbinsX()));

            c->cd(2);
            gPad->SetLeftMargin(0.15);
            h_raw_reco->SetLineColor(kOrange + 1);
            h_raw_reco->SetLineWidth(2);
            h_raw_reco->GetXaxis()->SetRangeUser(0.0, 2.0);
            h_raw_reco->GetXaxis()->SetTitle("p_{T}^{j1}/p_{T}^{truth}");
            h_raw_reco->GetYaxis()->SetTitle("Events");
            h_raw_reco->SetTitle("Raw reco/truth (combined MC)");
            h_raw_reco->Draw("HIST E");
            lt.DrawLatex(0.18, 0.85, Form("mean = %.4f", h_raw_reco->GetMean()));
            lt.DrawLatex(0.18, 0.80, Form("rms  = %.4f", h_raw_reco->GetRMS()));
            lt.DrawLatex(0.18, 0.75, Form("nbins = %d", h_raw_reco->GetNbinsX()));

            c->SaveAs(Form("%s/raw_truth_distributions.pdf", basePath.Data()));

            delete h_raw_truth;
            delete h_raw_reco;
        }

        double Sigma_truth_ref, SigmaErr_truth_ref, Sigma_reco_truth, SigmaErr_reco_truth;
        // The truth distributions are narrower than R (sigma ~ 0.07 against
        // ~ 0.11), so they get their own data-driven seed instead of the data
        // window.
        {
            RDBFitResults tref, rt;
            DoFit(h_comb_truth_R_balance, "truth_ref", "Comb MC truth/ref", tref);
            DoFit(h_comb_reco_truth, "reco_truth", "Comb MC reco/truth", rt);
            Sigma_truth_ref = tref.sigma;
            SigmaErr_truth_ref = tref.sigmaErr;
            Sigma_reco_truth = rt.sigma;
            SigmaErr_reco_truth = rt.sigmaErr;
        }

        double disc_d = data_fit.sigma * data_fit.sigma - Sigma_truth_ref * Sigma_truth_ref;
        double disc_m = comb_fit.sigma * comb_fit.sigma - Sigma_truth_ref * Sigma_truth_ref;

        // A negative difference means the truth reference is wider than the
        // measured R: the closure is meaningless, so warn and leave sigma at 0.
        double sig_insitu = 0.0;
        if (disc_d > 0) {
            sig_insitu = std::sqrt(disc_d);
        } else {
            Log::warning("FitResponse",
                         Form("sigma_data^2 (%.5f) < sigma_truth_ref^2 (%.5f) -- closure invalid: "
                              "the truth reference is wider than data. Check truth-match quality.",
                              data_fit.sigma * data_fit.sigma, Sigma_truth_ref * Sigma_truth_ref));
        }

        double sig_MC_insitu = 0.0;
        if (disc_m > 0) {
            sig_MC_insitu = std::sqrt(disc_m);
        } else {
            Log::warning(
                "FitResponse",
                Form(
                    "sigma_MC^2 (%.5f) < sigma_truth_ref^2 (%.5f) -- sigma_truth dominates over MC R_response width. "
                    "Check population consistency between R_response and R_truth histograms.",
                    comb_fit.sigma * comb_fit.sigma, Sigma_truth_ref * Sigma_truth_ref));
        }
        double sig_MC = Sigma_reco_truth;

        double err_insitu = 0, err_MC_insitu = 0;
        if (sig_insitu > 0)
            err_insitu = std::sqrt(std::pow(data_fit.sigma * data_fit.sigmaErr / sig_insitu, 2)
                                   + std::pow(Sigma_truth_ref * SigmaErr_truth_ref / sig_insitu, 2));
        if (sig_MC_insitu > 0)
            err_MC_insitu = std::sqrt(std::pow(comb_fit.sigma * comb_fit.sigmaErr / sig_MC_insitu, 2)
                                      + std::pow(Sigma_truth_ref * SigmaErr_truth_ref / sig_MC_insitu, 2));

        // R_JES = mu_data / mu_MC
        jes = FitUtils::ComputeJESRatio(data_fit.mean, data_fit.meanErr, comb_fit.mean, comb_fit.meanErr);

        // R_JER = sigma_insitu / sigma^MC_insitu
        jer = FitUtils::ComputeJERRatio(data_fit.sigma, data_fit.sigmaErr, comb_fit.sigma, comb_fit.sigmaErr,
                                        Sigma_truth_ref, SigmaErr_truth_ref);

        if (!jer.valid) {
            Log::warning("FitResponse",
                         "Cannot compute b-JER ratio: negative variance in sigma_insitu or sigma^MC_insitu");
        }

        double delta_clos = 0, err_clos = 0;
        if (sig_insitu > 1e-9 && sig_MC_insitu > 1e-9) {
            delta_clos = (sig_MC - sig_MC_insitu) / sig_insitu;

            // Error propagation with the independent fit results as inputs.
            // sigma_truth_ref (x3) enters both B and C, so it is propagated
            // once through both instead of being counted twice:
            //   A = sigma^MC (reco/truth fit), x1 = sigma_data (only in C),
            //   x2 = sigma_MC (only in B),     x3 = sigma_truth_ref (B and C)
            //   Delta = (A - B(x2,x3)) / C(x1,x3)
            //   dDelta/dA  = 1/C
            //   dDelta/dx1 = -(A-B) x1 / C^3
            //   dDelta/dx2 = -x2 / (B C)
            //   dDelta/dx3 = (x3/C) [1/B + (A-B)/C^2]
            const double A = sig_MC;           // sigma^MC
            const double B = sig_MC_insitu;    // sigma^MC_insitu
            const double C = sig_insitu;       // sigma_insitu
            const double x1 = data_fit.sigma;  // sigma_data
            const double x2 = comb_fit.sigma;  // sigma of R in combined MC
            const double x3 = Sigma_truth_ref; // sigma_truth_ref
            const double sA = SigmaErr_reco_truth;
            const double sx1 = data_fit.sigmaErr;
            const double sx2 = comb_fit.sigmaErr;
            const double sx3 = SigmaErr_truth_ref;
            const double AmB = A - B;

            double dA = 1.0 / C;
            double dx1 = -AmB * x1 / (C * C * C);
            double dx2 = -x2 / (B * C);
            double dx3 = (x3 / C) * (1.0 / B + AmB / (C * C));

            err_clos = std::sqrt(dA * dA * sA * sA + dx1 * dx1 * sx1 * sx1 + dx2 * dx2 * sx2 * sx2
                                 + dx3 * dx3 * sx3 * sx3);
        }

        Log::info(
            "FitResponse",
            Form("JER: sig_insitu=%.4f+/-%.4f  sig^MC_insitu=%.4f+/-%.4f  sig^MC=%.4f+/-%.4f  Dcl=%.4f+/-%.4f",
                 sig_insitu, err_insitu, sig_MC_insitu, err_MC_insitu, sig_MC, SigmaErr_reco_truth,
                 delta_clos, err_clos));

        // JER summary plot
        {
            CanvasGuard c(Form("jer_%s_%s", cfg.key.c_str(), wp_label), "b-JER", 800, 600);
            c->SetLeftMargin(0.15);
            double vals[3] = {sig_insitu, sig_MC_insitu, sig_MC},
                   errs[3] = {err_insitu, err_MC_insitu, SigmaErr_reco_truth};
            double xp[3] = {1, 2, 3}, xe[3] = {0, 0, 0};
            TGraphErrors* g = new TGraphErrors(3, xp, vals, xe, errs);
            g->SetTitle(Form("b-JER (%s %s);;#sigma", cfg.key.c_str(), wp_label));
            g->SetMarkerStyle(21);
            g->SetMarkerSize(1.5);
            g->SetLineWidth(2);
            g->GetXaxis()->SetBinLabel(g->GetXaxis()->FindBin(1), "#sigma_{in situ}");
            g->GetXaxis()->SetBinLabel(g->GetXaxis()->FindBin(2), "#sigma^{MC}_{in situ}");
            g->GetXaxis()->SetBinLabel(g->GetXaxis()->FindBin(3), "#sigma^{MC}");
            g->GetXaxis()->SetLabelSize(0.06);
            double maxY = *std::max_element(vals, vals + 3) + *std::max_element(errs, errs + 3) + 0.05;
            g->GetYaxis()->SetRangeUser(0, maxY > 0.05 ? maxY : 0.5);
            g->Draw("AP");
            TLatex lt;
            lt.SetNDC();
            lt.SetTextFont(42);
            lt.SetTextSize(0.04);
            lt.DrawLatex(0.18, 0.75, Form("#Delta#sigma_{cl}=%.4f#pm%.4f", delta_clos, err_clos));
            if (jes.valid) lt.DrawLatex(0.18, 0.70, Form("R_{JES}=%.4f#pm%.4f", jes.ratio, jes.ratioErr));
            if (jer.valid) lt.DrawLatex(0.18, 0.65, Form("R_{JER}=%.4f#pm%.4f", jer.ratio, jer.ratioErr));
            AtlasUtils::DrawAtlasLabel(0.18, 0.85);
            c->SaveAs(Form("%s/JER_summary.pdf", basePath.Data()));
        }

        // JER overlay
        {
            CanvasGuard c(Form("jer_ov_%s_%s", cfg.key.c_str(), wp_label), "R_response overlay", 800, 600);
            TH1F* hd = (TH1F*)h_combined_R_balance->Clone("ov_r");
            if (hd->Integral() > 0) hd->Scale(1. / hd->Integral());
            TH1F* ht = (TH1F*)h_comb_truth_R_balance->Clone("ov_t");
            if (ht->Integral() > 0) ht->Scale(1. / ht->Integral());
            TH1F* hr = (TH1F*)h_comb_reco_truth->Clone("ov_rt");
            if (hr->Integral() > 0) hr->Scale(1. / hr->Integral());
            hd->SetLineColor(kRed);
            ht->SetLineColor(kAzure + 2);
            hr->SetLineColor(kOrange + 1);
            for (auto* h : {hd, ht, hr}) {
                h->SetLineWidth(2);
                h->SetFillStyle(0);
            }
            double ym = std::max({hd->GetMaximum(), ht->GetMaximum(), hr->GetMaximum()}) * 1.4;
            hd->GetYaxis()->SetRangeUser(0, ym);
            hd->SetTitle((std::string("Combined MC;") + cfg.RAxisLabel() + ";Normalised").c_str());
            hd->Draw("HIST");
            ht->Draw("HIST SAME");
            hr->Draw("HIST SAME");
            TLegend* l = new TLegend(0.55, 0.65, 0.92, 0.88);
            l->SetBorderSize(0);
            l->SetTextSize(0.028);
            l->AddEntry(hd, Form("%s (#sigma=%.3f)", cfg.RAxisLabel().c_str(), comb_fit.sigma), "l");
            l->AddEntry(ht, Form("R_{truth} (#sigma=%.3f)", Sigma_truth_ref), "l");
            l->AddEntry(hr, Form("reco/truth (#sigma=%.3f)", Sigma_reco_truth), "l");
            l->Draw();
            AtlasUtils::DrawAtlasLabel(0.18, 0.85);
            c->SaveAs(Form("%s/JER_R_balance_overlay.pdf", basePath.Data()));
        }

        // JER text file
        {
            ofstream f(Form("%s/JER_results.txt", basePath.Data()));
            f << "b-JER -- " << cfg.key << " WP: " << wp_label << "\n"
              << std::string(60, '=') << "\n"
              << "Truth-match cut: dR < " << cfg.maxTruthDR << "\n\n"
              << "Sample statistics (raw / Sum_w / N_eff):\n"
              << Form(
                     "  MC signal:   raw=%lld/%lld  Sum_w=%.2f/%.2f  N_eff=%.0f/%.0f  (weighted pass = %.1f%%)\n",
                     n_mc_truth_pass, n_mc_total, sumw_mc_pass, sumw_mc_total,
                     Neff(sumw_mc_pass, sumw2_mc_pass), Neff(sumw_mc_total, sumw2_mc_total),
                     sumw_mc_total > 0 ? 100.0 * sumw_mc_pass / sumw_mc_total : 0.0)
              << Form(
                     "  gamma+jets:  raw=%lld/%lld  Sum_w=%.2f/%.2f  N_eff=%.0f/%.0f  (weighted pass = %.1f%%)\n",
                     n_bkg_truth_pass, n_bkg_total, sumw_bkg_pass, sumw_bkg_total,
                     Neff(sumw_bkg_pass, sumw2_bkg_pass), Neff(sumw_bkg_total, sumw2_bkg_total),
                     sumw_bkg_total > 0 ? 100.0 * sumw_bkg_pass / sumw_bkg_total : 0.0)
              << "\n"
              << Form("sigma_data         = %.6f +/- %.6f\n", data_fit.sigma, data_fit.sigmaErr)
              << Form("sigma_mc_combined  = %.6f +/- %.6f\n", comb_fit.sigma, comb_fit.sigmaErr)
              << Form("sigma_truth_ref    = %.6f +/- %.6f\n", Sigma_truth_ref, SigmaErr_truth_ref)
              << Form("sigma_reco_truth   = %.6f +/- %.6f\n", Sigma_reco_truth, SigmaErr_reco_truth)
              << "\nDerived:\n"
              << Form("sigma_insitu    = %.6f +/- %.6f\n", sig_insitu, err_insitu)
              << Form("sigma^MC_insitu = %.6f +/- %.6f\n", sig_MC_insitu, err_MC_insitu)
              << Form("sigma^MC        = %.6f +/- %.6f\n", sig_MC, SigmaErr_reco_truth)
              << Form("Delta_closure   = %.6f +/- %.6f (%.2f%%)\n", delta_clos, err_clos, delta_clos * 100)
              << Form("R_JES = data_mean/MC_mean         = %.4f +/- %.4f\n", jes.ratio, jes.ratioErr)
              << Form("R_JER = sigma_insitu/sigma^MC_insitu = %.4f +/- %.4f\n", jer.ratio, jer.ratioErr)
              << Form("JER SF closure test ~ %.4f\n", 1.0 + delta_clos);
            f.close();
        }

        // Same JER chain separately for signal, background and their sum.
        // Histograms with fewer than 20 entries are skipped (sigma = 0).
        auto QuickFit = [&](TH1F* h, const char* tag, const char* label, double seedLo = -1,
                            double seedHi = -1) -> std::pair<double, double> {
            if (!h || h->GetEntries() < 20) return {0, 0};
            RDBFitResults out;
            auto fit = DoFit(h, tag, label, out, seedLo, seedHi);
            return {fit.sigma, fit.sigmaErr};
        };

        // PerSampFit and ComputePerSampleJER are in Response_Fit_Core.h.
        // The truth fits use the data-driven seed, so the signal and background
        // cores each get a window sized to their own width.
        auto [st_sig, st_sig_e] = QuickFit(h_mc_truth_R_balance, "truth_sig", "Sig R_{truth,DB}");
        auto [sr_sig, sr_sig_e] = QuickFit(h_mc_reco_truth, "reco_truth_sig", "Sig reco/truth");
        auto [st_bkg, st_bkg_e] = QuickFit(h_bkg_truth_R_balance, "truth_bkg", "Bkg R_{truth,DB}");
        auto [sr_bkg, sr_bkg_e] = QuickFit(h_bkg_reco_truth, "reco_truth_bkg", "Bkg reco/truth");
        // The combined truth widths are the truth_ref / reco_truth fits above
        double st_comb = Sigma_truth_ref, st_comb_e = SigmaErr_truth_ref;
        double sr_comb = Sigma_reco_truth, sr_comb_e = SigmaErr_reco_truth;

        PerSampFit psf[3];
        psf[0] = ComputePerSampleJER(data_fit.sigma, data_fit.sigmaErr, mc_fit.sigma, mc_fit.sigmaErr, st_sig,
                                     st_sig_e, sr_sig, sr_sig_e);
        psf[1] = ComputePerSampleJER(data_fit.sigma, data_fit.sigmaErr, bkg_fit.sigma, bkg_fit.sigmaErr,
                                     st_bkg, st_bkg_e, sr_bkg, sr_bkg_e);
        psf[2] = ComputePerSampleJER(data_fit.sigma, data_fit.sigmaErr, comb_fit.sigma, comb_fit.sigmaErr,
                                     st_comb, st_comb_e, sr_comb, sr_comb_e);
        const char* psLabels[3] = {"Signal", "Bkg", "Combined"};

        // The six widths for each sample on one plot
        {
            CanvasGuard c(Form("csum_%s_%s", cfg.key.c_str(), wp_label), "sigma summary", 900, 600);
            c->SetLeftMargin(0.12);
            c->SetBottomMargin(0.14);
            c->SetGridy();
            double ymax = 0;
            for (int i = 0; i < 3; i++)
                ymax = std::max(
                    {ymax, psf[i].sig_RDB + psf[i].sig_RDB_err, psf[i].sig_RDB_MC + psf[i].sig_RDB_MC_err,
                     psf[i].sig_truth + psf[i].sig_truth_err, psf[i].sig_MC + psf[i].sig_MC_err,
                     psf[i].sig_insitu + psf[i].sig_insitu_err, psf[i].sig_MCi + psf[i].sig_MCi_err});
            TH1F* fr = new TH1F(Form("fr_%s_%s", cfg.key.c_str(), wp_label), ";;#sigma", 3, 0, 2);
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
            // The first two labels use the channel's R name (R_{DB} / R_{MJB}).
            std::string rLab = cfg.RAxisLabel();
            std::string qL[6] = {std::string("#sigma(") + rLab + ")_{data}",
                                 std::string("#sigma(") + rLab + ")_{MC}",
                                 "#sigma(R_{truth})_{MC}",
                                 "#sigma^{MC}(reco/truth)",
                                 "#sigma_{insitu}",
                                 "#sigma^{MC}_{insitu}"};
            double sOff[3] = {0.333, 1.0, 1.5};
            TLegend* leg = new TLegend(0.72, 0.62, 0.92, 0.92);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            leg->SetTextSize(0.024);
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
                leg->AddEntry(gr, qL[q].c_str(), "lp");
            }
            leg->Draw();
            AtlasUtils::DrawAtlasLabel(0.18, 0.85);
            TLatex lt;
            lt.SetNDC();
            lt.SetTextSize(0.030);
            lt.DrawLatex(0.18, 0.75, Form("%s  WP: %s", cfg.label.c_str(), wp_label));
            c->SaveAs(Form("%s/JER_sigma_summary_markers.pdf", basePath.Data()));
        }

        // Per-sample table on screen and in JER_per_sample.txt
        {
            std::cout << "\n+====================================================================+\n"
                      << "|  JER per sample -- " << cfg.key << " WP " << wp_label << "\n"
                      << "+====================================================================+\n";
            for (int i = 0; i < 3; i++)
                std::cout << Form(
                    "| %-10s sig_d=%.4f sig_MC=%.4f sig_tr=%.4f sig^MC=%.4f sig_in=%.4f sig^MCi=%.4f Dcl=%.3f",
                    psLabels[i], psf[i].sig_RDB, psf[i].sig_RDB_MC, psf[i].sig_truth, psf[i].sig_MC,
                    psf[i].sig_insitu, psf[i].sig_MCi, psf[i].dC)
                          << "\n";
            std::cout << "+====================================================================+\n\n";

            ofstream f(Form("%s/JER_per_sample.txt", basePath.Data()));
            f << "JER per sample -- " << cfg.key << " WP: " << wp_label << "\n"
              << std::string(95, '=') << "\n";
            f << Form("%-25s%22s%22s%22s\n", "", "Signal", "Bkg", "Combined") << std::string(95, '-') << "\n";
            auto WR = [&](const char* l, double a, double ae, double b, double be, double cc, double ce) {
                f << Form("%-25s%10.5f+/-%8.5f%10.5f+/-%8.5f%10.5f+/-%8.5f\n", l, a, ae, b, be, cc, ce);
            };
            WR("sig(R_response)_data", psf[0].sig_RDB, psf[0].sig_RDB_err, psf[1].sig_RDB, psf[1].sig_RDB_err,
               psf[2].sig_RDB, psf[2].sig_RDB_err);
            WR("sig(R_response)_MC", psf[0].sig_RDB_MC, psf[0].sig_RDB_MC_err, psf[1].sig_RDB_MC,
               psf[1].sig_RDB_MC_err, psf[2].sig_RDB_MC, psf[2].sig_RDB_MC_err);
            WR("sig(R_truth)_MC", psf[0].sig_truth, psf[0].sig_truth_err, psf[1].sig_truth,
               psf[1].sig_truth_err, psf[2].sig_truth, psf[2].sig_truth_err);
            WR("sig^MC", psf[0].sig_MC, psf[0].sig_MC_err, psf[1].sig_MC, psf[1].sig_MC_err, psf[2].sig_MC,
               psf[2].sig_MC_err);
            WR("sig_insitu", psf[0].sig_insitu, psf[0].sig_insitu_err, psf[1].sig_insitu,
               psf[1].sig_insitu_err, psf[2].sig_insitu, psf[2].sig_insitu_err);
            WR("sig^MC_insitu", psf[0].sig_MCi, psf[0].sig_MCi_err, psf[1].sig_MCi, psf[1].sig_MCi_err,
               psf[2].sig_MCi, psf[2].sig_MCi_err);
            WR("Dsig_clos", psf[0].dC, psf[0].dC_err, psf[1].dC, psf[1].dC_err, psf[2].dC, psf[2].dC_err);
            f.close();
        }
    } // end if hasTruthBranches

    // R stack with data/MC ratio
    // Uses the "_all" histograms (no truth requirement), so the MC describes
    // the same population as data. The fit curves are those of the
    // truth-matched fits.
    {
        CanvasGuard c(Form("cstack_%s_%s", cfg.key.c_str(), wp_label), "R_response Stack", 800, 800);
        TPad* p1 = new TPad("p1", "", 0, 0.3, 1, 1);
        p1->SetBottomMargin(0.02);
        p1->SetTopMargin(0.08);
        p1->SetLeftMargin(0.12);
        p1->SetRightMargin(0.05);
        p1->Draw();
        p1->cd();

        TH1F* hd = (TH1F*)h_data_R_balance->Clone();
        hd->SetMarkerStyle(21);
        hd->SetMarkerColor(kBlack);
        TH1F* hs = (TH1F*)h_mc_R_balance_all->Clone();
        hs->SetFillColor(kRed);
        hs->SetFillStyle(1001);
        TH1F* hb = (TH1F*)h_multijets_R_balance_all->Clone();
        hb->SetFillColor(kBlue);
        hb->SetFillStyle(1001);
        THStack* stk = new THStack(Form("stk_%s_%s", cfg.key.c_str(), wp_label), ";;Events");
        hb->SetBit(kCanDelete, kFALSE);
        hs->SetBit(kCanDelete, kFALSE);
        stk->Add(hb);
        stk->Add(hs);
        stk->Draw("HIST");
        stk->SetMaximum(stk->GetMaximum() * 1.5);
        stk->GetYaxis()->SetTitle("Events");
        stk->GetXaxis()->SetLabelSize(0);
        HistUtils::StyleDataPoints(hd);
        hd->Draw("EP SAME");
        // Data and combined MC fits, drawn with the fitted model over their fit windows.
        TF1* gfD = FitUtils::MakeModelTF1(fit_data, "gfD");
        gfD->SetLineColor(kViolet - 4);
        gfD->SetLineWidth(2);
        gfD->Draw("same");
        TF1* gfC = FitUtils::MakeModelTF1(fit_comb, "gfC");
        gfC->SetLineColor(kCyan - 7);
        gfC->SetLineWidth(2);
        gfC->Draw("same");
        TH1F* hMC = (TH1F*)stk->GetStack()->Last();
        TGraphAsymmErrors* eb = HistUtils::MakeUncertaintyBand(hMC);
        eb->Draw("2 SAME");
        TLegend* leg = new TLegend(0.6, 0.65, 0.89, 0.89);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.025);
        leg->AddEntry(hd, "Data", "lep");
        leg->AddEntry(hs, sigLab, "f");
        leg->AddEntry(hb, bkgLab, "f");
        leg->AddEntry(gfD, "Data Fit", "l");
        leg->AddEntry(gfC, "Total MC Fit", "l");
        leg->AddEntry(eb, "MC Stat. Unc.", "f");
        leg->Draw();
        AtlasUtils::DrawAtlasLabel(0.18, 0.85);
        c->cd();
        TPad* p2 = new TPad("p2", "", 0, 0, 1, 0.3);
        p2->SetTopMargin(0.05);
        p2->SetBottomMargin(0.35);
        p2->SetLeftMargin(0.12);
        p2->SetRightMargin(0.05);
        p2->Draw();
        p2->cd();
        TH1F* hR = (TH1F*)hd->Clone();
        hR->Divide(h_combined_R_balance_all);
        hR->SetTitle("");
        hR->GetYaxis()->SetTitle("Data/MC");
        hR->GetYaxis()->SetTitleSize(0.14);
        hR->GetYaxis()->SetTitleOffset(0.4);
        hR->GetYaxis()->SetLabelSize(0.10);
        hR->GetYaxis()->SetRangeUser(0.5, 1.5);
        hR->GetYaxis()->SetNdivisions(505);
        hR->GetXaxis()->SetTitle(cfg.RAxisLabel().c_str());
        hR->GetXaxis()->SetTitleSize(0.14);
        hR->GetXaxis()->SetLabelSize(0.10);
        hR->Draw("EP");
        TGraphAsymmErrors* rb = HistUtils::MakeRatioUncertaintyBand(h_combined_R_balance_all);
        rb->Draw("2 SAME");
        TLine* ln = new TLine(p2->GetUxmin(), 1, p2->GetUxmax(), 1);
        ln->SetLineStyle(2);
        ln->SetLineWidth(2);
        ln->Draw("SAME");
        hR->Draw("EP SAME");
        c->SaveAs(Form("%s/R_DB_stack_ratio_combined.pdf", basePath.Data()));
    }

    // Kinematic distributions
    {
        const char* sLabels[kNS] = {sigLab, bkgLab, "Data", "Combined MC"};
        const char* sTags[kNS] = {"sig", "bkg", "data", "comb"};
        int sColors[kNS] = {kRed, kBlue, kBlack, kGreen + 2};

        auto DrawKinematic = [&](TH1F* hists[kNS], const char* var, const char* xTitle) {
            for (int s = 0; s < kNS; s++) {
                if (hists[s]->GetEntries() < 10) continue;
                CanvasGuard c(Form("c_%s_%s_%s_%s", var, sTags[s], cfg.key.c_str(), wp_label), var, 800, 600);
                hists[s]->SetLineColor(sColors[s]);
                hists[s]->SetLineWidth(2);
                hists[s]->SetMinimum(0);
                hists[s]->Draw("HIST");
                TLegend* l = new TLegend(0.55, 0.75, 0.92, 0.88);
                l->SetBorderSize(0);
                l->SetTextSize(0.028);
                l->AddEntry(hists[s], sLabels[s], "l");
                l->Draw();
                AtlasUtils::DrawAtlasLabel(0.18, 0.85);
                c->SaveAs(Form("%s/%s_%s.pdf", basePath.Data(), var, sTags[s]));
            }
            CanvasGuard c(Form("c_%s_ov_%s_%s", var, cfg.key.c_str(), wp_label), "overlay", 800, 600);
            TH1F* hn[kNS];
            double ym = 0;
            for (int s = 0; s < kNS; s++) {
                hn[s] = (TH1F*)hists[s]->Clone(Form("n_%s_%d", var, s));
                if (hn[s]->Integral() > 0) hn[s]->Scale(1. / hn[s]->Integral());
                hn[s]->SetLineColor(sColors[s]);
                hn[s]->SetLineWidth(2);
                ym = std::max(ym, hn[s]->GetMaximum());
            }
            hn[0]->SetMaximum(ym * 1.5);
            hn[0]->SetMinimum(0);
            hn[0]->GetXaxis()->SetTitle(xTitle);
            hn[0]->GetYaxis()->SetTitle("Normalised");
            hn[0]->Draw("HIST");
            for (int s = 1; s < kNS; s++) hn[s]->Draw("HIST SAME");
            TLegend* l = new TLegend(0.55, 0.65, 0.92, 0.88);
            l->SetBorderSize(0);
            l->SetTextSize(0.028);
            for (int s = 0; s < kNS; s++) l->AddEntry(hn[s], sLabels[s], "l");
            l->Draw();
            AtlasUtils::DrawAtlasLabel(0.18, 0.85);
            c->SaveAs(Form("%s/%s_overlay.pdf", basePath.Data(), var));
        };

        TH1F* jpt[kNS] = {h_jet_pt[0], h_jet_pt[1], h_jet_pt[2], h_jet_pt[3]};
        TH1F* phpt[kNS] = {h_ref_pt[0], h_ref_pt[1], h_ref_pt[2], h_ref_pt[3]};
        DrawKinematic(jpt, "jet_pt", "p_{T}^{jet} [GeV]");
        DrawKinematic(phpt, "photon_pt", "p_{T}^{#gamma} [GeV]");

        // dR(reco, truth), MC only, with the truth-match cut as a vertical line
        for (int s : {0, 1, 3}) {
            if (h_dR_samp[s]->GetEntries() < 10) continue;
            CanvasGuard c(Form("c_dR_%s_%s_%s", sTags[s], cfg.key.c_str(), wp_label), "dR", 800, 600);
            h_dR_samp[s]->SetLineColor(sColors[s]);
            h_dR_samp[s]->SetLineWidth(2);
            h_dR_samp[s]->SetMinimum(0);
            h_dR_samp[s]->Draw("HIST");
            TLine* lcut = new TLine(cfg.maxTruthDR, 0, cfg.maxTruthDR, h_dR_samp[s]->GetMaximum() * 0.9);
            lcut->SetLineColor(kRed);
            lcut->SetLineStyle(2);
            lcut->SetLineWidth(2);
            lcut->Draw("SAME");
            TLegend* l = new TLegend(0.55, 0.75, 0.92, 0.88);
            l->SetBorderSize(0);
            l->SetTextSize(0.028);
            l->AddEntry(h_dR_samp[s], sLabels[s], "l");
            l->AddEntry(lcut, Form("#DeltaR < %.2f cut", cfg.maxTruthDR), "l");
            l->Draw();
            AtlasUtils::DrawAtlasLabel(0.18, 0.85);
            c->SaveAs(Form("%s/truth_match_dR_%s.pdf", basePath.Data(), sTags[s]));
        }
        {
            CanvasGuard c(Form("c_dR_ov_%s_%s", cfg.key.c_str(), wp_label), "dR ov", 800, 600);
            TH1F* hn[kNS] = {nullptr};
            double ym = 0;
            for (int s : {0, 1, 3}) {
                hn[s] = (TH1F*)h_dR_samp[s]->Clone(Form("n_dR_%d", s));
                if (hn[s]->Integral() > 0) hn[s]->Scale(1. / hn[s]->Integral());
                hn[s]->SetLineColor(sColors[s]);
                hn[s]->SetLineWidth(2);
                ym = std::max(ym, hn[s]->GetMaximum());
            }
            hn[0]->SetMaximum(ym * 1.5);
            hn[0]->SetMinimum(0);
            hn[0]->GetXaxis()->SetTitle("#DeltaR(reco,truth)");
            hn[0]->GetYaxis()->SetTitle("Normalised");
            hn[0]->Draw("HIST");
            hn[1]->Draw("HIST SAME");
            hn[3]->Draw("HIST SAME");
            TLine* lcut = new TLine(cfg.maxTruthDR, 0, cfg.maxTruthDR, ym * 1.4);
            lcut->SetLineColor(kRed);
            lcut->SetLineStyle(2);
            lcut->SetLineWidth(2);
            lcut->Draw("SAME");
            TLegend* l = new TLegend(0.55, 0.65, 0.92, 0.88);
            l->SetBorderSize(0);
            l->SetTextSize(0.028);
            for (int s : {0, 1, 3}) l->AddEntry(hn[s], sLabels[s], "l");
            l->AddEntry(lcut, Form("#DeltaR < %.2f cut", cfg.maxTruthDR), "l");
            l->Draw();
            AtlasUtils::DrawAtlasLabel(0.18, 0.85);
            c->SaveAs(Form("%s/truth_match_dR_overlay.pdf", basePath.Data()));
        }
    }

    // Summary table
    {
        ofstream tbl(Form("%s/R_DB_summary_table.txt", basePath.Data()));
        tbl << "R_response SUMMARY -- " << cfg.key << " WP: " << wp_label << " (inclusive, no mass cuts)\n"
            << "Truth-match cut for JER closure: dR < " << cfg.maxTruthDR << "\n"
            << std::string(80, '=') << "\n\n"
            << Form("%-15s %10s %10s %10s %10s %10s\n", "Sample", "Mean", "MeanErr", "Sigma", "SigmaErr",
                    "Chi2/NDF")
            << Form("%-15s %10.4f %10.4f %10.4f %10.4f %10.3f\n", "MC Signal", mc_fit.mean, mc_fit.meanErr,
                    mc_fit.sigma, mc_fit.sigmaErr, mc_fit.chi2ndf)
            << Form("%-15s %10.4f %10.4f %10.4f %10.4f %10.3f\n", "Data", data_fit.mean, data_fit.meanErr,
                    data_fit.sigma, data_fit.sigmaErr, data_fit.chi2ndf)
            << Form("%-15s %10.4f %10.4f %10.4f %10.4f %10.3f\n", "gamma+jets", bkg_fit.mean, bkg_fit.meanErr,
                    bkg_fit.sigma, bkg_fit.sigmaErr, bkg_fit.chi2ndf)
            << Form("%-15s %10.4f %10.4f %10.4f %10.4f %10.3f\n", "Combined MC", comb_fit.mean,
                    comb_fit.meanErr, comb_fit.sigma, comb_fit.sigmaErr, comb_fit.chi2ndf)
            << (jes.valid ? Form("R_JES (data/MC mean)  = %.4f +/- %.4f\n", jes.ratio, jes.ratioErr)
                          : std::string("R_JES (data/MC mean)  = [not computed]\n"))
            << (jer.valid
                    ? Form("R_JER (data/MC sigma) = %.4f +/- %.4f\n", jer.ratio, jer.ratioErr)
                    : std::string(
                          "R_JER (data/MC sigma) = [not computed — no truth branches or negative variance]\n"));
        tbl.close();
        Log::info("FitResponse", Form("Summary: %s/R_DB_summary_table.txt", basePath.Data()));
    }

    // Bin-count recommendation for VariantConfig.h
    // One reco count (nBinsR) and one truth count (nBinsTruth, from the
    // narrower of truth/ref and reco/truth, i.e. the larger count). The values
    // to paste are the fit-friendly ones; Scott is printed as a reference.
    // Both are computed from the width and N_eff of the filled combined-MC
    // histograms, so they hardly depend on the binning in use. Printed under
    // "[BINNING-REC]" and written to binning_recommendation.txt.
    {
        auto scott = [&](TH1F* h) -> int {
            return (h && h->GetEntries() > 20) ? BinningConfig::ComputeScott(h, 0.0, 2.0, 10, 200) : 0;
        };
        auto fitfr = [&](TH1F* h) -> int {
            return (h && h->GetEntries() > 20) ? BinningConfig::ComputeFitFriendly(h, 0.0, 2.0, 10, 200) : 0;
        };
        int scR = scott(h_combined_R_balance);
        int scT = std::max(scott(h_comb_truth_R_balance), scott(h_comb_reco_truth));
        int ffR = fitfr(h_combined_R_balance);
        int ffT = std::max(fitfr(h_comb_truth_R_balance), fitfr(h_comb_reco_truth));
        auto emit = [&](std::ostream& os) {
            os << "==== BINNING RECOMMENDATION -- " << cfg.key << "  WP " << wp_label
               << "  (inclusive) ====\n";
            os << "  fit-friendly = ~N bins in +/-2sigma core (set EASYBJETS_BINNING_CORE_BINS, default 10);"
                  " Scott = density-optimal (reference).\n";
            os << "  reco(fit=" << ffR << ", scott=" << scR << ")  "
               << "truth(fit=" << ffT << ", scott=" << scT << ")\n";
            os << "  paste into common/VariantConfig.h (fit-friendly; truth often capped ~100):\n";
            os << "  configs[\"" << cfg.key << "\"].nBinsR     = " << ffR << ";\n";
            os << "  configs[\"" << cfg.key << "\"].nBinsTruth = " << ffT << ";\n";
        };
        std::cout << "\n+----------------------------- [BINNING-REC] -----------------------------+\n";
        emit(std::cout);
        std::cout << "+-------------------------------------------------------------------------+\n\n";
        std::ofstream rec(Form("%s/binning_recommendation.txt", basePath.Data()));
        if (rec.is_open()) {
            emit(rec);
            rec.close();
            Log::info("FitResponse",
                      Form("Binning recommendation: %s/binning_recommendation.txt", basePath.Data()));
        }
    }

    // Cleanup
    fmc->Close();
    delete fmc;
    fdata->Close();
    delete fdata;
    fmultijets->Close();
    delete fmultijets;
    for (int s = 0; s < SampleConfig::kNSamples; s++)
        if (fextra_rdb[s]) {
            fextra_rdb[s]->Close();
            delete fextra_rdb[s];
        }
    Log::info("FitResponse", Form("Completed %s WP %s", cfg.key.c_str(), wp_label));
}

// Entry points
inline void RunVariant(const std::string& key) {
    // Full keys ("DB", "DB_*", "MJB", "MJB_*") are used as they are. A short
    // name such as "bJER10v00" is taken as a DB variant, so MJB variants must
    // always be given with their full key.
    std::string resolved = key;
    const bool is_db = (key == "DB") || (key.rfind("DB_", 0) == 0);
    const bool is_mjb = (key == "MJB") || (key.rfind("MJB_", 0) == 0);
    if (!is_db && !is_mjb) resolved = "DB_" + key;
    VariantConfig cfg = GetVariant(resolved);
    for (const auto& wp : AnalysisWPConfig::StandardWPs()) RunRdbAnalysis(cfg, wp.c_str());
    Log::info("FitResponse", "All WPs done for " + cfg.key);
}

inline void FitResponse_Run(const std::string& variant, const std::string& channel) {
    AtlasUtils::LoadAtlasStyle();
    TH1::AddDirectory(kFALSE);
    gROOT->SetBatch(kTRUE);
    if (variant == "ALL") {
        for (const auto& key : GetPrimaryVariantKeys(channel)) RunVariant(key);
        return; // "ALL" is not a variant key; GetVariant would throw on it
    }
    RunVariant(variant);
}

#endif // FIT_RESPONSE_H
