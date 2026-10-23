/**
 * ============================================================================
 *  Mass_Distributions.h - data/MC distributions of the selected large-R jet
 * ============================================================================
 * ===========================================
 * @author:    Thiziri AMEZZA
 * @email:     thiziri.amezza@cern.ch
 * @date:      10 March 2026
 * ===========================================
 *
 * For one variant and each WP, reads the slim ntuples of data, the Z->bb signal
 * MC and the main background MC (gamma+jets for DB, multijets for MJB), and:
 *   1. fills the jet mass (several ranges), the balance R (R_DB or R_MJB) and the
 *      asymmetry, for every sample
 *   2. fits the sideband scale factor SF(m) = p0 + p1*m of the background MC to
 *      data on the 50-150 GeV mass histogram (SFUtils), and applies it bin by bin
 *      to the background mass and event by event to its R and asymmetry
 *   3. draws data/MC stack plots with a ratio panel
 *   4. splits the background by the b and c hadrons in the jet and plots the
 *      composition: BB (>= 2 b), BC (1 b + c), BL (1 b, no c), CC (>= 2 c),
 *      CL (1 c), LL (light)
 * The fitted SF is also saved as sf2_params_<wp>.txt, which
 * Response_Fit_Inclusive.h reuses instead of refitting.
 * Mass_Distributions_pTbins.h reuses HistogramManager with a pT window.
 *
 * How to run (from ZbbGamma/macros for DB, from ZbbJets/macros for MJB):
 *   export EASYBJETS_RUN=run3
 *   export EASYBJETS_DATA_DIR=$HOME/lpnatlas_data/ZbbGamma    (ZbbJets for MJB)
 *   export EASYBJETS_MASS_CUT=on       for this analysis: jets with 50-150 GeV
 *   root -l -b -q 'Mass_Distributions.cxx()'                nominal variant
 *   root -l -b -q 'Mass_Distributions.cxx("DB_bJER10v00")'  one variant (MJB_... in ZbbJets)
 *   root -l -b -q 'Mass_Distributions.cxx("ALL")'           all four variants
 *
 * Output: <PlotDir>/<variant>/Mass_Distributions/<wp>/
 *   stack plots (Large_R_Jet_m*.pdf, Asymmetry*.pdf, R_DB_*.pdf),
 *   scale_factor_fit.pdf, sf2_params_<wp>.txt,
 *   background_composition.{pdf,root,txt}, R_DB_split.pdf, mass_split.pdf
 *
 * Options (environment variables, read once per ROOT session):
 *   EASYBJETS_RUN         run2 | run3 (default)
 *   EASYBJETS_DATA_DIR    folder holding NtupleSlim_Run2/3 of this channel
 *                         (default ../../DataSets/<ZbbGamma|ZbbJets>)
 *   EASYBJETS_MASS_CUT    on: keep jets with 50 <= m <= 150 GeV; lo:hi for another
 *                         window; not set = no cut (default). See MassCutConfig.h
 *   EASYBJETS_SIDEBAND1_LOW, _SIDEBAND1_HIGH, _SIDEBAND2_LOW, _SIDEBAND2_MAX
 *                         override the sidebands of the SF fit (GeV)
 *   EASYBJETS_FIT_MODEL   only changes the output folder (Run{2,3}Plots_<model>);
 *                         leave it unset here
 *
 * ==========================================================================
 * @version:  1.0
 * ==========================================================================
 */
#ifndef MASS_RESPONSE_H
#define MASS_RESPONSE_H

#include "VariantConfig.h"
#include "CommonUtils.h"
#include "SFUtils.h"
#include "MassCutConfig.h"
#include <fstream>

// One histogram per sample.
using SampleHists = HistUtils::HistArray<SampleConfig::kNSamples>;

// Balance title and symbol for the channel: DB or MJB.
inline std::string BalanceTitleFor(const std::string& ch) {
    return (ch == "MJB") ? "Multi-Jet Balance Ratio" : "Direct Balance Ratio";
}
inline std::string RTexLabelFor(const std::string& ch) { return (ch == "MJB") ? "R_{MJB}" : "R_{DB}"; }

// Fills, scales and draws all the histograms of one variant and WP.
class HistogramManager {
  private:
    const std::string ctx = "Mass_Distributions";

    // One histogram per sample. The suffixes
    // are the same variable with different ranges / binnings.
    SampleHists h_ljet_m;
    SampleHists h_ljet_m_50_350;
    SampleHists h_ljet_m_50_150;
    SampleHists h_ljet_m_50_150_2; // 20 bins, the one the SF is fitted on
    SampleHists h_Asym;
    SampleHists h_Asym_2;
    SampleHists h_R_DB_1;
    SampleHists h_R_DB_2;
    SampleHists h_R_DB_3;
    SampleHists h_R_DB;

    std::unique_ptr<TH1D> h_background_composition;
    std::vector<std::unique_ptr<TH1F>> h_R_DB_bins; // R per flavour category
    std::vector<std::unique_ptr<TH1F>> h_mass_bins; // mass per flavour category

    std::string workingPoint;
    std::string outputDir;
    VariantConfig cfg;

    bool sampleLoaded[SampleConfig::kNSamples];

    SFUtils::LinearSF sf;

    // Optional pT window, used by Mass_Distributions_pTbins. Off by default,
    // so the full pT range is used.
    double ptLo_ = -1e30, ptHi_ = 1e30;
    bool ptCutActive_ = false;
    std::string ptDispLabel_; // e.g. "200 < p_{T} < 300 GeV", written on the plots

    // Read the variant pT branch when a pT window is set. For data, use
    // "ljet_pt" if the variant branch does not exist (same as for the mass).
    void SetPtBranchIfNeeded(TTree* tree, std::vector<float>*& ljet_pt, bool isData) {
        if (!ptCutActive_) return;
        if (isData) {
            bool hasV = Validate::HasBranch(tree, cfg.branch_pt);
            tree->SetBranchAddress(hasV ? cfg.branch_pt.c_str() : "ljet_pt", &ljet_pt);
        } else {
            tree->SetBranchAddress(cfg.branch_pt.c_str(), &ljet_pt);
        }
    }
    bool PassPtWindow(std::vector<float>* ljet_pt, int idx) const {
        if (!ptCutActive_) return true;
        if (!ljet_pt || idx < 0 || idx >= (int)ljet_pt->size()) return false;
        double pt = (*ljet_pt)[idx];
        return (pt >= ptLo_ && pt < ptHi_);
    }
    // Write the pT-bin label on the plot (nothing for the full pT range).
    void DrawPtLabel() {
        if (ptDispLabel_.empty()) return;
        TLatex pl;
        pl.SetNDC();
        pl.SetTextFont(42);
        pl.SetTextSize(0.035);
        pl.DrawLatex(0.20, 0.80, ptDispLabel_.c_str());
    }

    // Create the histogram of one variable for every sample.
    void InitializeHistogramGroup(SampleHists& hists, const std::string& variable, const std::string& title,
                                  int nbins, float min, float max) {
        auto names = SampleConfig::HistNames(cfg);
        for (int i = 0; i < (int)names.size(); ++i) {
            hists.Set(i, HistUtils::CreateDetached("h_" + names[i] + "_" + variable, title + ";Entries",
                                                   nbins, min, max));
        }
    }

    void SetHistogramColors(SampleHists& hists) {
        auto colors = SampleConfig::Colors();
        auto fills = SampleConfig::FillStyles();
        for (int i = 0; i < SampleConfig::kNSamples; ++i) {
            if (hists[i]) {
                hists[i]->SetLineColor(colors[i]);
                hists[i]->SetFillColor(colors[i]);
                hists[i]->SetFillStyle(fills[i]);
            }
        }
    }

  public:
    HistogramManager(const std::string& wp, const std::string& dir, const VariantConfig& c)
        : workingPoint(wp), outputDir(dir), cfg(c), sf() {
        for (int i = 0; i < SampleConfig::kNSamples; ++i) sampleLoaded[i] = false;
        InitializeHistograms();
    }

    ~HistogramManager() = default;

    // Restrict to a pT window (Mass_Distributions_pTbins). dispLabel is written
    // on every plot.
    void SetPtWindow(double lo, double hi, const std::string& dispLabel = "") {
        ptLo_ = lo;
        ptHi_ = hi;
        ptCutActive_ = true;
        ptDispLabel_ = dispLabel;
    }

    void InitializeHistograms() {
        InitializeHistogramGroup(h_ljet_m, "ljet_m", "Large-R Jet Mass;M^{J} [GeV]", 30, 50,
                                 (int)cfg.mass_xmax);
        if (cfg.channel == "DB")
            InitializeHistogramGroup(h_ljet_m_50_350, "ljet_m_50_350", "Large-R Jet Mass;M^{J} [GeV]", 15, 50,
                                     150);
        InitializeHistogramGroup(h_ljet_m_50_150, "ljet_m_50_150", "Large-R Jet Mass;M^{J} [GeV]", 30, 50,
                                 150);
        InitializeHistogramGroup(h_ljet_m_50_150_2, "ljet_m_50_150_2", "Large-R Jet Mass;M^{J} [GeV]", 20, 50,
                                 150);
        InitializeHistogramGroup(h_Asym, "Asym", "Asymmetry;Asymmetry", 50, -1, 1);
        InitializeHistogramGroup(h_Asym_2, "Asym_2", "Asymmetry;Asymmetry", 20, -0.35, 0.35);
        // R_{DB} or R_{MJB} depending on the channel.
        std::string rLab = RTexLabelFor(cfg.channel);
        std::string rTitle = "Direct Balance Ratio;" + rLab;
        if (cfg.channel == "MJB") rTitle = "Multi-Jet Balance Ratio;" + rLab;
        InitializeHistogramGroup(h_R_DB_1, "R_DB_1", rTitle.c_str(), 35, 0, 2);
        InitializeHistogramGroup(h_R_DB_2, "R_DB_2", rTitle.c_str(), 50, 0, 2);
        InitializeHistogramGroup(h_R_DB_3, "R_DB_3", rTitle.c_str(), 55, 0, 2);
        InitializeHistogramGroup(h_R_DB, cfg.branch_R.c_str(), rTitle.c_str(), 60, 0, 2);

        SetHistogramColors(h_ljet_m);
        if (cfg.channel == "DB") SetHistogramColors(h_ljet_m_50_350);
        SetHistogramColors(h_ljet_m_50_150);
        SetHistogramColors(h_ljet_m_50_150_2);
        SetHistogramColors(h_Asym);
        SetHistogramColors(h_Asym_2);
        SetHistogramColors(h_R_DB_1);
        SetHistogramColors(h_R_DB_2);
        SetHistogramColors(h_R_DB_3);
        SetHistogramColors(h_R_DB);

        h_background_composition =
            std::make_unique<TH1D>("h_background_composition",
                                   (SampleConfig::DisplayNames(cfg)[SampleConfig::kMultijets]
                                    + " Background Composition;Category;Events")
                                       .c_str(),
                                   6, 0, 6);
        for (int i = 0; i < 6; ++i) {
            h_R_DB_bins.push_back(std::make_unique<TH1F>(("h_R_DB_bin_" + std::to_string(i)).c_str(),
                                                         "R_response;R_response;Events", 60, 0, 2));
            h_mass_bins.push_back(std::make_unique<TH1F>(("h_mass_bin_" + std::to_string(i)).c_str(),
                                                         "Mass;Mass [GeV];Events", 20, 50, 150));
        }
    }

    // Connect the branches. For data, use "ljet_m" if the variant mass branch
    // does not exist. The jet index is "index_tagged" or "selected_ljet_index".
    void SetBranchAddresses(TTree* tree, std::vector<float>*& ljet_m, double& weight, double& Asym,
                            double& R_response, int& selected_ljet_index, std::vector<int>*& ljet_n_bhadrons,
                            std::vector<int>*& ljet_n_chadrons, bool isData = false) {
        if (isData) {
            bool hasVariant = Validate::HasBranch(tree, cfg.branch_mass);
            tree->SetBranchAddress(hasVariant ? cfg.branch_mass.c_str() : "ljet_m", &ljet_m);
            if (!hasVariant) Log::info("Mass_Distributions", "Data mass branch: ljet_m [fallback]");
        } else {
            tree->SetBranchAddress(cfg.branch_mass.c_str(), &ljet_m);
        }
        if (Validate::HasBranch(tree, "index_tagged"))
            tree->SetBranchAddress("index_tagged", &selected_ljet_index);
        else if (Validate::HasBranch(tree, "selected_ljet_index"))
            tree->SetBranchAddress("selected_ljet_index", &selected_ljet_index);
        tree->SetBranchAddress("weight", &weight);
        tree->SetBranchAddress(cfg.branch_R.c_str(), &R_response);
        if (Validate::HasBranch(tree, "Asym")) {
            tree->SetBranchAddress("Asym", &Asym);
        }
        if (Validate::HasBranch(tree, "ljet_n_bhadrons") && Validate::HasBranch(tree, "ljet_n_chadrons")) {
            tree->SetBranchAddress("ljet_n_bhadrons", &ljet_n_bhadrons);
            tree->SetBranchAddress("ljet_n_chadrons", &ljet_n_chadrons);
        }
    }

    // Data and signal: fill all variables in one loop. "histograms" is the
    // list built in ProcessFiles (index 1 is null for MJB).
    void FillHistograms(TTree* tree, TH1F** histograms, bool isData) {
        std::vector<float>* ljet_m = nullptr;
        double weight = 0.0, R_response = 0.0, Asym = 0.0;
        int selected_ljet_index = -1;
        std::vector<int>* ljet_n_bhadrons = nullptr;
        std::vector<int>* ljet_n_chadrons = nullptr;

        bool hasAsymBranch = Validate::HasBranch(tree, "Asym");

        SetBranchAddresses(tree, ljet_m, weight, Asym, R_response, selected_ljet_index, ljet_n_bhadrons,
                           ljet_n_chadrons, isData);
        std::vector<float>* ljet_pt = nullptr;
        SetPtBranchIfNeeded(tree, ljet_pt, isData);

        Long64_t nentries = tree->GetEntries();
        for (Long64_t i = 0; i < nentries; ++i) {
            tree->GetEntry(i);
            double eventWeight = isData ? 1.0 : weight;

            if (!hasAsymBranch) {
                Asym = (1.0 + R_response > 1e-9) ? (1.0 - R_response) / (1.0 + R_response) : 0.0;
            }

            if (selected_ljet_index >= 0 && ljet_m && selected_ljet_index < (int)ljet_m->size()) {
                // Mass window preselection (EASYBJETS_MASS_CUT), off by default
                if (MassCutConfig::Reject(ljet_m, selected_ljet_index)) continue;
                if (!PassPtWindow(ljet_pt, selected_ljet_index)) continue;
                float mass = ljet_m->at(selected_ljet_index);
                histograms[0]->Fill(mass, eventWeight);
                if (histograms[1]) histograms[1]->Fill(mass, eventWeight);
                histograms[2]->Fill(mass, eventWeight);
                histograms[3]->Fill(mass, eventWeight);
                histograms[4]->Fill(Asym, eventWeight);
                histograms[5]->Fill(Asym, eventWeight);
                histograms[6]->Fill(R_response, eventWeight);
                histograms[7]->Fill(R_response, eventWeight);
                histograms[8]->Fill(R_response, eventWeight);
                histograms[9]->Fill(R_response, eventWeight);
            }
        }

        if (isData) Log::info(ctx, "Data entries: " + std::to_string(histograms[0]->GetEntries()));
    }

    // Background MC, first pass: only the mass, to fit the SF.
    void FillMassHistogramsOnly(TTree* tree, int sampleIdx) {
        std::vector<float>* ljet_m = nullptr;
        double weight = 0.0;
        int sel = -1;

        if (Validate::HasBranch(tree, "index_tagged"))
            tree->SetBranchAddress("index_tagged", &sel);
        else if (Validate::HasBranch(tree, "selected_ljet_index"))
            tree->SetBranchAddress("selected_ljet_index", &sel);
        bool isData = (sampleIdx == SampleConfig::kData);
        if (isData) {
            bool hasVariant = Validate::HasBranch(tree, cfg.branch_mass);
            tree->SetBranchAddress(hasVariant ? cfg.branch_mass.c_str() : "ljet_m", &ljet_m);
        } else {
            tree->SetBranchAddress(cfg.branch_mass.c_str(), &ljet_m);
        }
        tree->SetBranchAddress("weight", &weight);
        std::vector<float>* ljet_pt = nullptr;
        SetPtBranchIfNeeded(tree, ljet_pt, isData);
        for (Long64_t i = 0; i < tree->GetEntries(); ++i) {
            tree->GetEntry(i);
            double w = isData ? 1.0 : weight;
            if (sel < 0 || !ljet_m || sel >= (int)ljet_m->size()) continue;
            // Mass window preselection (EASYBJETS_MASS_CUT), off by default
            if (MassCutConfig::Reject(ljet_m, sel)) continue;
            if (!PassPtWindow(ljet_pt, sel)) continue;
            float mass = ljet_m->at(sel);
            h_ljet_m[sampleIdx]->Fill(mass, w);
            if (cfg.channel == "DB" && h_ljet_m_50_350[sampleIdx]) h_ljet_m_50_350[sampleIdx]->Fill(mass, w);
            h_ljet_m_50_150[sampleIdx]->Fill(mass, w);
            h_ljet_m_50_150_2[sampleIdx]->Fill(mass, w);
        }
        Log::info(ctx, Form("Mass-only fill done (%s, %lld entries)",
                            SampleConfig::FileNames(cfg)[sampleIdx].c_str(), tree->GetEntries()));
    }

    // Background MC, second pass: R and asymmetry with the per-event SF(m),
    // and the flavour composition.
    void FillResponseWithSF(TTree* tree, int sampleIdx) {
        std::vector<float>* ljet_m = nullptr;
        double weight = 0.0, R_response = 0.0, Asym = 0.0;
        int sel = -1;

        bool hasAsymBranch = Validate::HasBranch(tree, "Asym");
        bool isData = (sampleIdx == SampleConfig::kData);

        if (Validate::HasBranch(tree, "index_tagged"))
            tree->SetBranchAddress("index_tagged", &sel);
        else if (Validate::HasBranch(tree, "selected_ljet_index"))
            tree->SetBranchAddress("selected_ljet_index", &sel);
        if (isData) {
            bool hasVariant = Validate::HasBranch(tree, cfg.branch_mass);
            tree->SetBranchAddress(hasVariant ? cfg.branch_mass.c_str() : "ljet_m", &ljet_m);
        } else {
            tree->SetBranchAddress(cfg.branch_mass.c_str(), &ljet_m);
        }
        tree->SetBranchAddress("weight", &weight);
        tree->SetBranchAddress(cfg.branch_R.c_str(), &R_response);
        if (hasAsymBranch) tree->SetBranchAddress("Asym", &Asym);

        std::vector<int>* ljet_n_bhadrons = nullptr;
        std::vector<int>* ljet_n_chadrons = nullptr;
        bool hasComp = false;
        if (!isData && Validate::HasBranch(tree, "ljet_n_bhadrons")
            && Validate::HasBranch(tree, "ljet_n_chadrons")) {
            tree->SetBranchAddress("ljet_n_bhadrons", &ljet_n_bhadrons);
            tree->SetBranchAddress("ljet_n_chadrons", &ljet_n_chadrons);
            hasComp = true;
        } else if (!isData) {
            Log::warning(ctx, "ljet_n_bhadrons/ljet_n_chadrons not found in "
                                  + SampleConfig::FileNames(cfg)[sampleIdx]
                                  + " — background composition will be empty");
        }

        std::vector<float>* ljet_pt = nullptr;
        SetPtBranchIfNeeded(tree, ljet_pt, isData);

        for (Long64_t i = 0; i < tree->GetEntries(); ++i) {
            tree->GetEntry(i);
            if (sel < 0 || !ljet_m || sel >= (int)ljet_m->size()) continue;
            // Mass window preselection (EASYBJETS_MASS_CUT), off by default
            if (MassCutConfig::Reject(ljet_m, sel)) continue;
            if (!PassPtWindow(ljet_pt, sel)) continue;

            double mass = ljet_m->at(sel);
            double w = SFUtils::GetEventWeight(isData ? 1.0 : weight, mass, sf);

            if (!hasAsymBranch)
                Asym = (1.0 + R_response > 1e-9) ? (1.0 - R_response) / (1.0 + R_response) : 0.0;

            h_Asym[sampleIdx]->Fill(Asym, w);
            h_Asym_2[sampleIdx]->Fill(Asym, w);
            h_R_DB_1[sampleIdx]->Fill(R_response, w);
            h_R_DB_2[sampleIdx]->Fill(R_response, w);
            h_R_DB_3[sampleIdx]->Fill(R_response, w);
            h_R_DB[sampleIdx]->Fill(R_response, w);

            if (hasComp && ljet_n_bhadrons && ljet_n_chadrons && sel < (int)ljet_n_bhadrons->size()) {
                int nb = ljet_n_bhadrons->at(sel);
                int nc = ljet_n_chadrons->at(sel);
                // Flavour of the jet from its b / c hadron counts:
                // 0 BB, 1 BC, 2 BL, 3 CC, 4 CL, 5 LL
                int jet_bin = (nb >= 2)              ? 0
                              : (nb == 1 && nc >= 1) ? 1
                              : (nb == 1 && nc == 0) ? 2
                              : (nb == 0 && nc >= 2) ? 3
                              : (nb == 0 && nc == 1) ? 4
                                                     : 5;
                h_background_composition->Fill(jet_bin, w);
                if (jet_bin < (int)h_R_DB_bins.size()) {
                    h_R_DB_bins[jet_bin]->Fill(R_response, w);
                    h_mass_bins[jet_bin]->Fill(mass, w);
                }
            }
        }
    }

    // Fit SF(m) of the background MC to data in the mass sidebands (on the
    // 50-150 GeV, 20-bin mass histograms). Returns false (SF = 1) if data or
    // the background MC is missing.
    bool ComputeSFFromMass() {
        if (!sampleLoaded[SampleConfig::kData] || !sampleLoaded[SampleConfig::kMultijets]) {
            Log::warning(ctx, "Data or gamma+jets not loaded — SF=1 (no correction)");
            return false;
        }

        TH1F* hD = h_ljet_m_50_150_2[SampleConfig::kData];
        TH1F* hMJ = h_ljet_m_50_150_2[SampleConfig::kMultijets];

        sf = SFUtils::FitMassDependentSF(hD, hMJ, cfg.sideband1_low, cfg.sideband1_high, cfg.sideband2_low,
                                         cfg.sideband2_max, outputDir + "scale_factor_fit.pdf",
                                         cfg.key + " " + workingPoint);

        // Save the SF for Response_Fit_Inclusive. The sideband windows are saved
        // with it, so the reader can check they match before reusing the SF.
        if (sf.valid)
            SFUtils::SaveLinearSF(sf,
                                  outputDir + "sf2_params_" + PathBuilder::WPtoFileWP(workingPoint) + ".txt",
                                  cfg.key + " " + workingPoint, cfg.sideband1_low, cfg.sideband1_high,
                                  cfg.sideband2_low, static_cast<double>(cfg.sideband2_max));

        return sf.valid;
    }

    void ApplySFToMassHistograms() {
        const int mj = SampleConfig::kMultijets;
        SFUtils::ApplySFBinByBin(h_ljet_m[mj], sf);
        if (cfg.channel == "DB" && h_ljet_m_50_350[mj]) SFUtils::ApplySFBinByBin(h_ljet_m_50_350[mj], sf);
        SFUtils::ApplySFBinByBin(h_ljet_m_50_150[mj], sf);
        SFUtils::ApplySFBinByBin(h_ljet_m_50_150_2[mj], sf);
    }

    // Open the samples and run the two passes: fill data, signal and the
    // background mass; fit and apply the SF; then fill the background R.
    void ProcessFiles(const std::vector<std::string>& files) {
        TFile* rootFiles[SampleConfig::kNSamples] = {nullptr};
        TTree* trees[SampleConfig::kNSamples] = {nullptr};

        for (size_t i = 0; i < files.size(); ++i) {
            // Skip the "__unused_*" placeholder names that SampleConfig puts in
            // the minor-background slots this channel does not use.
            if (files[i].find("__unused_") != std::string::npos) continue;

            Log::info(ctx, "Opening: " + files[i]);
            TFile* f = TFile::Open(files[i].c_str(), "READ");
            if (!f || f->IsZombie()) {
                if (SampleConfig::IsOptional(i))
                    Log::debug(ctx, "  Optional — skipping");
                else
                    Log::error(ctx, "  Could not open: " + files[i]);
                if (f) delete f;
                continue;
            }
            rootFiles[i] = f;
            TTree* tree = Validate::GetTree(f, cfg.tree_name, ctx);
            if (tree) {
                trees[i] = tree;
                sampleLoaded[i] = true;
            }
        }

        // Pass 1: data and signal, all variables
        for (int i = 0; i < SampleConfig::kNSamples; ++i) {
            if (!trees[i] || i == SampleConfig::kMultijets) continue;
            bool isData = (i == SampleConfig::kData);
            TH1F* histograms[] = {h_ljet_m[i],        (cfg.channel == "DB" ? h_ljet_m_50_350[i] : nullptr),
                                  h_ljet_m_50_150[i], h_ljet_m_50_150_2[i],
                                  h_Asym[i],          h_Asym_2[i],
                                  h_R_DB_1[i],        h_R_DB_2[i],
                                  h_R_DB_3[i],        h_R_DB[i]};
            FillHistograms(trees[i], histograms, isData);
        }

        // Pass 1: background MC, mass only
        if (trees[SampleConfig::kMultijets])
            FillMassHistogramsOnly(trees[SampleConfig::kMultijets], SampleConfig::kMultijets);

        // Fit SF(m) and scale the background mass histograms
        ComputeSFFromMass();
        ApplySFToMassHistograms();

        // Pass 2: background MC, R and asymmetry with the per-event SF(m)
        if (trees[SampleConfig::kMultijets])
            FillResponseWithSF(trees[SampleConfig::kMultijets], SampleConfig::kMultijets);

        // Set negative bins of the MC (from negative weights) to 0
        for (int s = 1; s < SampleConfig::kNSamples; ++s) {
            ClampNegativeBins(h_ljet_m[s]);
            if (cfg.channel == "DB" && h_ljet_m_50_350[s]) ClampNegativeBins(h_ljet_m_50_350[s]);
            ClampNegativeBins(h_ljet_m_50_150[s]);
            ClampNegativeBins(h_ljet_m_50_150_2[s]);
            ClampNegativeBins(h_Asym[s]);
            ClampNegativeBins(h_Asym_2[s]);
            ClampNegativeBins(h_R_DB_1[s]);
            ClampNegativeBins(h_R_DB_2[s]);
            ClampNegativeBins(h_R_DB_3[s]);
            ClampNegativeBins(h_R_DB[s]);
        }
        for (auto& h : h_R_DB_bins) ClampNegativeBins(h.get());
        for (auto& h : h_mass_bins) ClampNegativeBins(h.get());

        // Close the files
        for (int i = 0; i < SampleConfig::kNSamples; ++i)
            if (rootFiles[i]) {
                rootFiles[i]->Close();
                delete rootFiles[i];
            }
    }


    // Draw and save all plots and the composition outputs.
    void SaveHistograms(bool createRatioPlots = true) {
        PathBuilder::MakeDir(outputDir);

        struct HistGroup {
            SampleHists& hists;
            std::string name;
        };
        // h_ljet_m_50_350 exists only for DB (never created for MJB).
        std::vector<HistGroup> groups;
        groups.push_back({h_ljet_m, "Large_R_Jet_m"});
        if (cfg.channel == "DB") groups.push_back({h_ljet_m_50_350, "Large_R_Jet_m_50_350"});
        groups.push_back({h_ljet_m_50_150, "Large_R_Jet_m_50_150"});
        groups.push_back({h_ljet_m_50_150_2, "Large_R_Jet_m_50_150_2"});
        groups.push_back({h_Asym, "Asymmetry"});
        groups.push_back({h_Asym_2, "Asymmetry_2"});
        groups.push_back({h_R_DB_1, "R_DB_1"});
        groups.push_back({h_R_DB_2, "R_DB_2"});
        groups.push_back({h_R_DB_3, "R_DB_3"});
        groups.push_back({h_R_DB, cfg.branch_R});

        for (auto& g : groups) SaveHistogramGroup(g.hists, g.name, createRatioPlots);

        // Flavour histograms in a ROOT file
        {
            SafeTFile outFile(outputDir + "background_composition.root", "RECREATE");
            if (outFile) {
                h_background_composition->Write();
                for (const auto& hist : h_R_DB_bins) hist->Write();
                for (const auto& hist : h_mass_bins) hist->Write();
            }
        }

        SaveBackgroundComposition();
        SaveSplitBackgroundHistograms();
    }

  private:
    void ClampNegativeBins(TH1* h) {
        if (!h) return;
        for (int bin = 1; bin <= h->GetNbinsX(); ++bin) {
            if (h->GetBinContent(bin) < 0) {
                h->SetBinContent(bin, 0);
                h->SetBinError(bin, 0);
            }
        }
    }

    // Data/MC stack plot of one variable, with a data/MC ratio panel if
    // createRatioPlot. Skipped if the data histogram does not exist or no
    // sample was loaded.
    void SaveHistogramGroup(SampleHists& hists, const std::string& name, bool createRatioPlot) {
        if (!hists[SampleConfig::kData]) return;

        bool anyMCLoaded = false;
        for (int i = 1; i < SampleConfig::kNSamples; ++i)
            if (sampleLoaded[i]) {
                anyMCLoaded = true;
                break;
            }
        if (!sampleLoaded[SampleConfig::kData] || !anyMCLoaded) return;

        CanvasGuard cg("canvas_" + name, name, 800, 800, outputDir + name + ".pdf");

        if (createRatioPlot) {
            auto pads = PadUtils::CreateRatioPads(cg.Canvas(), 0.30);
            pads.upper->cd();

            THStack* hs =
                new THStack(Form("hs_%s_%s_%s", name.c_str(), cfg.key.c_str(), workingPoint.c_str()), "");
            // BuildStack needs a plain TH1F* array
            TH1F* rawArr[SampleConfig::kNSamples];
            for (int s = 0; s < SampleConfig::kNSamples; ++s) rawArr[s] = hists[s];
            SampleConfig::BuildStack(hs, rawArr, sampleLoaded);
            hs->SetMaximum(hs->GetMaximum() * 1.5);
            hs->Draw("HIST");
            hs->GetXaxis()->SetTitle(hists[0]->GetXaxis()->GetTitle());
            hs->GetYaxis()->SetTitle(hists[0]->GetYaxis()->GetTitle());
            hs->GetXaxis()->SetTitleSize(0.05);
            hs->GetYaxis()->SetTitleSize(0.05);
            hs->GetXaxis()->SetLabelSize(0);

            if (!hs->GetStack() || hs->GetStack()->GetSize() == 0) {
                Log::warning(ctx, "Empty MC stack for " + name + " — skipping ratio panel");
                return;
            }

            TH1F* h_mc_total = static_cast<TH1F*>(hs->GetStack()->Last());
            if (!h_mc_total) {
                Log::warning(ctx, "Null summed MC histogram for " + name);
                return;
            }

            hs->SetMaximum(hs->GetMaximum() * 1.5);

            TGraphAsymmErrors* err_band = HistUtils::MakeUncertaintyBand(h_mc_total);
            err_band->Draw("2 same");

            HistUtils::StyleDataPoints(hists[0]);
            hists[0]->Draw("EP same");

            TLegend* legend = new TLegend(0.65, 0.65, 0.9, 0.93);
            legend->SetBorderSize(0);
            legend->SetTextSize(0.028);
            SampleConfig::AddLegendEntries(legend, rawArr, err_band, sampleLoaded, cfg);
            legend->Draw();

            AtlasUtils::DrawAtlasLabel();
            DrawPtLabel();

            pads.lower->cd();
            pads.lower->SetGridy(0);
            pads.lower->SetGridx(0);

            TH1F* h_ratio = HistUtils::CloneDetached(hists[0], "ratio_" + name);
            h_ratio->SetTitle("");
            h_ratio->Divide(h_mc_total);
            PadUtils::StyleRatioHist(h_ratio, std::string(hists[0]->GetXaxis()->GetTitle()), 0.5, 1.5);
            HistUtils::StyleDataPoints(h_ratio);
            h_ratio->Draw("EP");

            TGraphAsymmErrors* ratio_band = HistUtils::MakeRatioUncertaintyBand(h_mc_total);
            ratio_band->Draw("2 same");

            h_ratio->Draw("EP same");

            PadUtils::DrawRefLine(pads.lower, 1.0);
            pads.lower->Update();

        } else {
            THStack* hs =
                new THStack(Form("hs_%s_%s_%s", name.c_str(), cfg.key.c_str(), workingPoint.c_str()), "");
            TH1F* rawArr[SampleConfig::kNSamples];
            for (int s = 0; s < SampleConfig::kNSamples; ++s) rawArr[s] = hists[s];
            SampleConfig::BuildStack(hs, rawArr, sampleLoaded);
            hs->Draw("HIST");

            TH1F* h_mc_total_no_ratio = (TH1F*)hs->GetStack()->Last();
            TGraphAsymmErrors* err_band_no_ratio = HistUtils::MakeUncertaintyBand(h_mc_total_no_ratio);
            err_band_no_ratio->Draw("2 same");

            HistUtils::StyleDataPoints(hists[0]);
            hists[0]->Draw("EP same");
            AtlasUtils::DrawAtlasLabel();
            DrawPtLabel();
        }
    }

    // Flavour categories of the background composition: label and colour.
    static constexpr int kNFlavors = 6;

    const char* FlavorLabel(int i) const {
        static const char* labels[kNFlavors] = {"BB", "BC", "BL", "CC", "CL", "LL"};
        return (i >= 0 && i < kNFlavors) ? labels[i] : "?";
    }

    int FlavorColor(int i) const {
        static const int cols[kNFlavors] = {
            TColor::GetColor("#4e9a8a"), TColor::GetColor("#b3b3e6"), TColor::GetColor("#8fbc8f"),
            TColor::GetColor("#bdb85a"), TColor::GetColor("#f5c09a"), TColor::GetColor("#c79090"),
        };
        return (i >= 0 && i < kNFlavors) ? cols[i] : kBlack;
    }

    void CompositionStats(double yields[kNFlavors], double fracs[kNFlavors]) const {
        double total = 0;
        for (int i = 0; i < kNFlavors; ++i) {
            yields[i] = h_background_composition->GetBinContent(i + 1);
            total += yields[i];
        }
        for (int i = 0; i < kNFlavors; ++i) fracs[i] = (total > 0) ? 100.0 * yields[i] / total : 0.0;
    }

    std::string FlavorLegendLabel(int i, const double fracs[kNFlavors]) const {
        return Form("%s (%.2g%%)", FlavorLabel(i), fracs[i]);
    }

    // Bar chart of the background flavour composition, plus a text summary.
    void SaveBackgroundComposition() {
        double yields[kNFlavors], fracs[kNFlavors];
        CompositionStats(yields, fracs);

        {
            CanvasGuard cg("bkg_comp", "#gamma+jets Background Composition", 800, 800,
                           outputDir + "background_composition.pdf");
            cg->SetLeftMargin(0.15);
            cg->SetRightMargin(0.05);
            cg->SetTopMargin(0.07);
            cg->SetBottomMargin(0.14); // room for the x-axis title
            std::string suffix = "_" + cfg.key + "_" + workingPoint;

            THStack* hs =
                new THStack(("hs_comp" + suffix).c_str(), "#gamma+jets Background Composition;Flavor;Events");
            for (int i = 0; i < kNFlavors; ++i) {
                TH1D* slice =
                    new TH1D(Form("comp_slice_%d%s", i, suffix.c_str()), "", kNFlavors, 0, kNFlavors);
                slice->SetDirectory(nullptr);
                slice->SetBinContent(i + 1, yields[i]);
                slice->SetFillColor(FlavorColor(i));
                slice->SetLineColor(kBlack);
                slice->SetLineWidth(1);
                hs->Add(slice);
            }
            hs->Draw("BAR SAME");

            TH1D* hAxis =
                new TH1D(("hAxis_comp" + suffix).c_str(), "#gamma+jets Background Composition;Flavor;Events",
                         kNFlavors, 0, kNFlavors);
            hAxis->SetDirectory(nullptr);
            for (int i = 0; i < kNFlavors; ++i) hAxis->GetXaxis()->SetBinLabel(i + 1, FlavorLabel(i));
            hAxis->SetStats(0);
            hAxis->SetMaximum(hs->GetMaximum() * 1.5);
            hAxis->GetXaxis()->SetLabelSize(0.05);
            hAxis->GetXaxis()->SetTitleSize(0.045);
            hAxis->GetXaxis()->SetTitleOffset(1.1);
            hAxis->GetYaxis()->SetTitleOffset(1.5);
            hAxis->GetYaxis()->SetTitleSize(0.045);
            hAxis->GetYaxis()->SetLabelSize(0.04);
            hAxis->Draw("AXIS");
            hs->Draw("BAR SAME");

            TLegend* leg = new TLegend(0.65, 0.6, 0.9, 0.88);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            leg->SetTextSize(0.028);
            for (int i = 0; i < kNFlavors; ++i) {
                TH1D* dummy = new TH1D(Form("dummy_%d%s", i, suffix.c_str()), "", 1, 0, 1);
                dummy->SetDirectory(nullptr);
                dummy->SetFillColor(FlavorColor(i));
                dummy->SetLineColor(kBlack);
                leg->AddEntry(dummy, FlavorLegendLabel(i, fracs).c_str(), "f");
            }
            leg->Draw();
            // ATLAS label drawn by hand here: the plots have a wider left margin than
            // AtlasUtils::DrawAtlasLabel assumes.
            {
                TLatex l;
                l.SetNDC();
                l.SetTextFont(72);
                l.SetTextSize(0.04);
                l.DrawLatex(0.18, 0.88, "ATLAS");
                l.SetTextFont(42);
                l.DrawLatex(0.33, 0.88, "Internal");
                l.SetTextSize(0.035);
                l.DrawLatex(0.18, 0.83, RunConfig::LumiLabel().c_str());
                if (!ptDispLabel_.empty()) l.DrawLatex(0.18, 0.78, ptDispLabel_.c_str());
            }
        }

        SaveCompositionTxt(yields, fracs);
    }

    // Text table of the flavour yields and fractions, largest first.
    void SaveCompositionTxt(const double yields[kNFlavors], const double fracs[kNFlavors]) const {
        int order[kNFlavors] = {0, 1, 2, 3, 4, 5};
        std::sort(order, order + kNFlavors, [&](int a, int b) { return yields[a] > yields[b]; });
        double total = 0;
        for (int i = 0; i < kNFlavors; ++i) total += yields[i];

        std::ofstream txt(outputDir + "background_composition.txt");
        txt << "# gamma+jets background flavor composition (" << workingPoint << " WP)\n";
        txt << "# Generated by Mass_Distributions.cxx\n";
        txt << "# gamma+jets MC only (with sideband SF applied). Signal Z(bb)+gamma excluded.\n";
        txt << "# Per-variable fractions (mass_split, R_DB_split) differ and are logged to console.\n";
        txt << "#\n";
        txt << std::left << std::setw(6) << "# Bin" << std::left << std::setw(8) << "Flavor" << std::right
            << std::setw(18) << "Yield (weighted)" << std::right << std::setw(14) << "Fraction (%)"
            << "\n";
        txt << std::string(48, '-') << "\n";
        for (int k = 0; k < kNFlavors; ++k) {
            int i = order[k];
            txt << std::left << std::setw(6) << i << std::left << std::setw(8) << FlavorLabel(i) << std::right
                << std::fixed << std::setprecision(1) << std::setw(18) << yields[i] << std::right
                << std::fixed << std::setprecision(3) << std::setw(14) << fracs[i] << "\n";
        }
        txt << std::string(48, '-') << "\n";
        txt << std::left << std::setw(14) << "Total" << std::right << std::fixed << std::setprecision(1)
            << std::setw(18) << total << std::right << std::setw(14) << "100.000"
            << "\n";
        txt.close();
        Log::info(ctx, "Composition summary: " + outputDir + "background_composition.txt");
    }

    // R and mass of the background MC stacked by flavour, with data on top.
    void SaveSplitBackgroundHistograms() {
        auto drawStack = [&](const std::vector<std::unique_ptr<TH1F>>& hists, const std::string& fname,
                             const std::string& xTitle, TH1F* hData) {
            if (hists.size() != (size_t)kNFlavors) return;

            double integrals[kNFlavors], plotFracs[kNFlavors];
            double total = 0;
            for (int i = 0; i < kNFlavors; ++i) {
                integrals[i] = hists[i]->Integral();
                total += integrals[i];
            }
            for (int i = 0; i < kNFlavors; ++i)
                plotFracs[i] = (total > 0) ? 100.0 * integrals[i] / total : 0.0;

            int order[kNFlavors] = {0, 1, 2, 3, 4, 5};
            std::sort(order, order + kNFlavors, [&](int a, int b) { return integrals[a] > integrals[b]; });

            CanvasGuard cg(fname, fname, 800, 800, outputDir + fname + ".pdf");
            // Wide left margin so 5-digit y labels do not overlap the title.
            cg->SetLeftMargin(0.19);
            cg->SetRightMargin(0.05);
            cg->SetTopMargin(0.07);
            cg->SetBottomMargin(0.13);

            THStack* hs =
                new THStack(Form("hs_%s_%s_%s", fname.c_str(), cfg.key.c_str(), workingPoint.c_str()), "");
            for (int k = kNFlavors - 1; k >= 0; --k) {
                int i = order[k];
                hists[i]->SetFillColor(FlavorColor(i));
                hists[i]->SetLineColor(kBlack);
                hists[i]->SetLineWidth(1);
                hists[i]->SetFillStyle(1001);
                hs->Add(hists[i].get());
            }
            hs->Draw("HIST");
            hs->GetXaxis()->SetTitle(xTitle.c_str());
            hs->GetYaxis()->SetTitle("Events");
            hs->GetXaxis()->SetTitleSize(0.05);
            hs->GetYaxis()->SetTitleSize(0.05);
            hs->GetXaxis()->SetLabelSize(0.04);
            hs->GetYaxis()->SetLabelSize(0.04);
            hs->GetYaxis()->SetTitleOffset(1.6);
            hs->GetXaxis()->SetTitleOffset(1.1);

            double ymax = hs->GetMaximum();
            if (hData && hData->GetMaximum() > ymax) ymax = hData->GetMaximum();
            hs->SetMaximum(ymax * 1.5);

            // MC statistical uncertainty band
            TH1F* h_mc_total_flav = (TH1F*)hs->GetStack()->Last();
            TGraphAsymmErrors* err_band_flav = HistUtils::MakeUncertaintyBand(h_mc_total_flav);
            err_band_flav->Draw("2 SAME");

            if (hData && hData->GetEntries() > 0) {
                HistUtils::StyleDataPoints(hData);
                hData->Draw("EP SAME");
            }

            TLegend* leg = new TLegend(0.65, 0.6, 0.9, 0.88);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            leg->SetTextSize(0.028);
            if (hData && hData->GetEntries() > 0) leg->AddEntry(hData, "Data", "lep");
            for (int k = 0; k < kNFlavors; ++k) {
                int i = order[k];
                TH1F* dummy = new TH1F(
                    Form("dummy_%s_%d_%s_%s", fname.c_str(), i, cfg.key.c_str(), workingPoint.c_str()), "", 1,
                    0, 1);
                dummy->SetDirectory(nullptr);
                dummy->SetFillColor(FlavorColor(i));
                dummy->SetLineColor(kBlack);
                leg->AddEntry(dummy, FlavorLegendLabel(i, plotFracs).c_str(), "f");
            }
            leg->AddEntry(err_band_flav, "MC Stat. Unc.", "f");
            leg->Draw();
            // ATLAS label drawn by hand, shifted for the wider left margin
            {
                TLatex l;
                l.SetNDC();
                l.SetTextFont(72);
                l.SetTextSize(0.04);
                l.DrawLatex(0.22, 0.88, "ATLAS");
                l.SetTextFont(42);
                l.DrawLatex(0.37, 0.88, "Internal");
                l.SetTextSize(0.035);
                l.DrawLatex(0.22, 0.83, RunConfig::LumiLabel().c_str());
            }

            Log::info(ctx, fname + " flavor composition:");
            for (int k = 0; k < kNFlavors; ++k) {
                int i = order[k];
                Log::info(ctx, Form("  %s: %.2f%%", FlavorLabel(i), plotFracs[i]));
            }
        };

        TH1F* dataRDB = sampleLoaded[SampleConfig::kData] ? h_R_DB[SampleConfig::kData] : nullptr;
        TH1F* dataMass = sampleLoaded[SampleConfig::kData] ? h_ljet_m_50_150_2[SampleConfig::kData] : nullptr;

        drawStack(h_R_DB_bins, "R_DB_split", cfg.RAxisLabel(), dataRDB);
        drawStack(h_mass_bins, "mass_split", "M^{J} [GeV]", dataMass);
    }
};

// Called by Mass_Distributions.cxx in each channel's macros/folder
// (usage in the file header). "ALL" runs every variant of the channel.
inline void Mass_Response_Run(const std::string& variant, const std::string& channel) {
    AtlasUtils::LoadAtlasStyle();
    MassCutConfig::PrintBanner("Mass_Response");

    if (variant == "ALL") {
        for (const auto& key : GetPrimaryVariantKeys(channel)) {
            Log::info("Mass_Distributions", "Running for: " + key);
            Mass_Response_Run(key, channel);
        }
        return;
    }

    VariantConfig cfg = GetVariant(variant);
    auto workingPoints = AnalysisWPConfig::StandardWPs();

    for (const auto& wp : workingPoints) {
        Log::info("Mass_Distributions", "Processing WP " + wp);
        std::string outDir = PathBuilder::OutputDirWP(cfg, "Mass_Distributions", wp);
        PathBuilder::MakeDir(outDir);
        std::vector<std::string> inputFiles = SampleConfig::AllSamplePaths(cfg, wp);

        HistogramManager histManager(wp, outDir, cfg);
        histManager.ProcessFiles(inputFiles);
        histManager.SaveHistograms();

        Log::info("Mass_Distributions", "Completed WP " + wp);
    }
}

#endif // MASS_RESPONSE_H
