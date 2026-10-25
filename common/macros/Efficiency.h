/**
 * ============================================================================
 *  Efficiency.h - selection efficiency cut by cut (DB and MJB)
 * ============================================================================
 * ===========================================
 * @author:    Thiziri AMEZZA
 * @email:     thiziri.amezza@cern.ch
 * @date:      18 April 2026
 * ===========================================
 *
 * Reads the cutflow histogram "h_cutflow" written by the slimming producers
 * and turns it into an efficiency: events after cut i / events at "Initial".
 * Samples: data, the Z(bb) signal MC and the main background MC
 * (gamma+jets for DB, multijets for MJB).
 *
 * Plots, for each variant, in <PlotDir>/<variant>/Efficiency/:
 *   EfficiencyComparison_<wp>_{log,linear}.pdf
 *       one working point, the three samples overlaid
 *   WPComparison_<sample>_{log,linear}.pdf
 *       one sample, all working points overlaid (WPConfig::DefaultWPs)
 * and, in the console, a table of the efficiency after each cut.
 *
 * Error bars (see MakeEffHist):
 *   data  Clopper-Pearson 68.3% interval (unweighted counts), the larger side
 *         used as a symmetric error
 *   MC    binomial error with event weights, from the sum of weights squared
 *         stored as the cutflow bin errors
 *
 * The cuts and their number come from CutFlowUtils::DefaultCutNames
 * (DB: 12 cuts, MJB: 8 cuts). Nothing is refitted or reselected here, so
 * EASYBJETS_MASS_CUT has no effect on the numbers.
 *
 * How to run (from ZbbGamma/macros or ZbbJets/macros):
 *
 *   export EASYBJETS_RUN=run3                      run2 or run3
 *   export EASYBJETS_DATA_DIR=$HOME/lpnatlas_data/ZbbGamma    (ZbbJets for MJB)
 *
 *   rootsigma 'Efficiency.cxx()'                  nominal (default)
 *   rootsigma 'Efficiency.cxx("DB_bJER10v00")'    one variant
 *   rootsigma 'Efficiency.cxx("ALL")'             all variants
 *   (MJB: "MJB_bJER10v00", ...)
 *
 * ==========================================================================
 * @version:  2.1
 * ==========================================================================
 */

#ifndef EFFICIENCY_H
#define EFFICIENCY_H

#include "VariantConfig.h"            // found through rootlogon.C (common/ on the include path)
#include "CommonUtils.h"
#include "MassCutConfig.h"
#include "TEfficiency.h"              // ClopperPearson for efficiency bin errors

namespace EfficiencyImpl {

// ----------------------------------------------------------------------------
// Efficiency histogram from a cutflow:
//   bin i = events passing cuts 0..i / initial events (bin 1),
//   with the cut names as bin labels.
//
//   weighted = false (data):
//     Clopper-Pearson 68.3% interval on the (integer) counts. The interval
//     is asymmetric; its larger side is used as a symmetric bin error.
//
//   weighted = true (MC):
//     Efficiency from sums of event weights. The error takes into account
//     that the passing events are a subset of the initial ones:
//
//       eff      = W_pass / W_init
//       var(eff) = [(1 - 2*eff) * S_pass + eff^2 * S_init] / W_init^2
//
//     W = sum of event weights         (bin content)
//     S = sum of squared event weights (bin error squared)
//     First-order error propagation, with the weights taken as fixed.
// ----------------------------------------------------------------------------
inline TH1F* MakeEffHist(TH1F*                            hCF,
                          const std::vector<std::string>&  cutNames,
                          const std::string&               name,
                          const std::string&               title,
                          bool                             weighted) {
    if (!hCF) return nullptr;
    int nCuts = (int)cutNames.size();
    TH1F* hEff = HistUtils::CreateDetached(name.c_str(), title.c_str(),
                                            nCuts, 0, nCuts);
    hEff->SetDirectory(nullptr);
    const double initial = hCF->GetBinContent(1);
    const double S_init = hCF->GetBinError(1) * hCF->GetBinError(1);
    const double cl = 0.683;  // 1 sigma
    for (int b = 1; b <= nCuts && b <= hCF->GetNbinsX(); ++b) {
        double passed = hCF->GetBinContent(b);
        double eff = (initial > 0) ? passed / initial : 0.0;
        hEff->SetBinContent(b, eff);
        if (initial <= 0) {
            // nothing to divide by: bin left at 0
        } else if (weighted) {
            double S_pass = hCF->GetBinError(b) * hCF->GetBinError(b);
            double var = ((1.0 - 2.0 * eff) * S_pass + eff * eff * S_init) / (initial * initial);
            hEff->SetBinError(b, std::sqrt(std::max(0.0, var)));
        } else {
            int N = (int)std::round(initial);
            int k = std::max(0, std::min(N, (int)std::round(passed)));
            double lo = TEfficiency::ClopperPearson(N, k, cl, false);
            double hi = TEfficiency::ClopperPearson(N, k, cl, true);
            hEff->SetBinError(b, std::max(eff - lo, hi - eff));
        }
        hEff->GetXaxis()->SetBinLabel(b, cutNames[b - 1].c_str());
    }
    return hEff;
}

// ----------------------------------------------------------------------------
// Draws several efficiency histograms on one canvas and saves it to outPath.
// logY: log y axis (range from the smallest non-zero value) 
// ----------------------------------------------------------------------------
inline void DrawEffCanvas(const std::vector<TH1F*>&        hists,
                           const std::vector<std::string>&  labels,
                           const std::vector<int>&          colors,
                           const std::string&               yTitle,
                           const std::string&               outPath,
                           bool                             logY) {
    if (hists.empty()) return;

    std::string cname = "cEff_" + outPath.substr(outPath.rfind('/') + 1);
    CanvasGuard cg(cname, yTitle, 1000, 800, outPath);
    TCanvas* c = cg.Canvas();
    c->SetBottomMargin(0.22);
    if (logY) c->SetLogy();

    // y range from the non-zero contents of all histograms
    double yMin = 1.0, yMax = 0.0;
    for (auto* h : hists) {
        if (!h) continue;
        for (int b = 1; b <= h->GetNbinsX(); ++b) {
            double v = h->GetBinContent(b);
            if (v > 0) {
                yMin = std::min(yMin, v);
                yMax = std::max(yMax, v);
            }
        }
    }
    yMin = logY ? std::max(yMin * 0.3, 1e-7) : 0.0;
    yMax = logY ? std::min(yMax * 8.0, 10.0) : 1.30;

    TLegend* leg = new TLegend(0.62, 0.65, 0.92, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.030);

    bool first = true;
    for (size_t i = 0; i < hists.size(); ++i) {
        if (!hists[i]) continue;
        hists[i]->SetLineColor(colors[i]);
        hists[i]->SetMarkerColor(colors[i]);
        hists[i]->SetLineWidth(3);
        hists[i]->SetMarkerStyle(20 + (int)i);
        hists[i]->SetMarkerSize(1.4);
        hists[i]->GetYaxis()->SetRangeUser(yMin, yMax);
        hists[i]->GetYaxis()->SetTitle(yTitle.c_str());
        hists[i]->GetYaxis()->SetTitleOffset(1.3);
        hists[i]->GetXaxis()->SetLabelSize(0.03);
        hists[i]->GetXaxis()->SetLabelOffset(0.01);
        hists[i]->GetXaxis()->SetTitleSize(0.05);
        hists[i]->GetXaxis()->SetTitleOffset(1.8);
        hists[i]->GetXaxis()->LabelsOption("d");
        hists[i]->Draw(first ? "LP" : "LP SAME");
        leg->AddEntry(hists[i], labels[i].c_str(), "lp");
        first = false;
    }
    leg->Draw();
    AtlasUtils::DrawAtlasLabel();
}

// ----------------------------------------------------------------------------
// All efficiency plots and tables of one variant. The channel comes from
// cfg.channel, so the same code serves DB and MJB.
// ----------------------------------------------------------------------------
inline void RunEfficiency(const VariantConfig& cfg) {
    const std::string ctx = "Efficiency";
    Log::info(ctx, "=== " + cfg.label + " (" + cfg.key + ") ===");

    auto workingPoints = WPConfig::DefaultWPs();
    auto fileNames     = SampleConfig::FileNames(cfg);
    auto dispNames     = SampleConfig::DisplayNames(cfg);
    auto colors        = SampleConfig::Colors();

    std::string outDir = PathBuilder::OutputDir(cfg, "Efficiency");
    PathBuilder::MakeDir(outDir);

    // ==========================================================================
    // PART 1 - one canvas per WP: data, signal and background overlaid
    // ==========================================================================
    for (const auto& wp : workingPoints) {
        std::string fileWP  = PathBuilder::WPtoFileWP(wp);
        std::string labelWP = PathBuilder::WPtoLabel(wp);
        auto cutNames   = CutFlowUtils::DefaultCutNames(cfg.channel, labelWP);
        const int nCuts = (int)cutNames.size();

        // Cutflows of the three samples
        auto cfSet = CutFlowUtils::LoadAllCutflows(cfg, wp, ctx);

        if (!cfSet.IsLoaded(SampleConfig::kData)) {
            Log::warning(ctx, "Skipping WP " + fileWP + "-data file missing");
            continue;
        }

        TH1F* hEff[SampleConfig::kNSamples] = {};
        for (int s = 0; s < SampleConfig::kNSamples; ++s) {
            if (!cfSet.IsLoaded(s)) continue;
            hEff[s] = MakeEffHist(
                cfSet[s], cutNames,
                Form("heff_%s_%s_%s", cfg.key.c_str(), fileWP.c_str(), fileNames[s].c_str()),
                Form("%s;Cut;Efficiency", labelWP.c_str()),
                s != SampleConfig::kData);
        }

        std::vector<TH1F*>       toPlot;
        std::vector<std::string>  legLabels;
        std::vector<int>          legColors;
        for (int s = 0; s < SampleConfig::kNSamples; ++s) {
            if (!hEff[s]) continue;
            toPlot.push_back(hEff[s]);
            legLabels.push_back(dispNames[s]);
            legColors.push_back(colors[s]);
        }

        DrawEffCanvas(toPlot, legLabels, legColors,
                      "Selection Efficiency",
                      outDir + "EfficiencyComparison_" + fileWP + "_log.pdf",
                      true);
        DrawEffCanvas(toPlot, legLabels, legColors,
                      "Selection Efficiency",
                      outDir + "EfficiencyComparison_" + fileWP + "_linear.pdf",
                      false);

        // Console table: efficiency after each cut, one column per sample
        std::cout << "\n[Efficiency] " << cfg.key << "  WP=" << labelWP << "\n";
        std::cout << std::string(72, '-') << "\n";
        std::cout << std::left << std::setw(36) << "Cut";
        for (int s = 0; s < SampleConfig::kNSamples; ++s)
            if (hEff[s]) std::cout << std::right << std::setw(12) << dispNames[s];
        std::cout << "\n" << std::string(72, '-') << "\n";
        for (int b = 1; b <= nCuts; ++b) {
            std::cout << std::left << std::setw(36) << cutNames[b - 1];
            for (int s = 0; s < SampleConfig::kNSamples; ++s) {
                if (!hEff[s]) continue;
                std::cout << std::right << std::fixed << std::setprecision(4)
                          << std::setw(12) << hEff[s]->GetBinContent(b);
            }
            std::cout << "\n";
        }
        std::cout << std::string(72, '-') << "\n";

        // cfSet auto-cleans
        for (auto* h : hEff) delete h;
    }

    // ==========================================================================
    // PART 2 - one canvas per sample: all WPs overlaid
    // ==========================================================================
    const std::vector<int> wpColors = {kRed, kMagenta, kOrange+7,
                                       kBlue, kGreen+2, kCyan+2};

    // Efficiencies of every sample and WP, each file read once per WP
    std::vector<TH1F*>       wpHistsAll[SampleConfig::kNSamples];
    std::vector<std::string> wpLabelsAll[SampleConfig::kNSamples];
    std::vector<int>         wpColsAll[SampleConfig::kNSamples];

    for (size_t iw = 0; iw < workingPoints.size(); ++iw) {
        const std::string& wp = workingPoints[iw];
        std::string fileWP    = PathBuilder::WPtoFileWP(wp);
        std::string labelWP   = PathBuilder::WPtoLabel(wp);

        auto cfSet = CutFlowUtils::LoadAllCutflows(cfg, wp, ctx);
        auto wpCutNames = CutFlowUtils::DefaultCutNames(cfg.channel, labelWP);

        for (int s = 0; s < SampleConfig::kNSamples; ++s) {
            if (!cfSet.IsLoaded(s)) continue;
            // hEff copies the values, so it stays valid when cfSet is deleted
            TH1F* hEff = MakeEffHist(
                cfSet[s], wpCutNames,
                Form("hwpcomp_%s_%s_%s", cfg.key.c_str(), fileWP.c_str(), fileNames[s].c_str()),
                Form("%s;Cut;Efficiency", dispNames[s].c_str()),
                s != SampleConfig::kData);
            if (!hEff) continue;
            wpHistsAll[s].push_back(hEff);
            wpLabelsAll[s].push_back(labelWP);
            wpColsAll[s].push_back(wpColors[iw % wpColors.size()]);
        }
    }

    // One canvas per sample
    for (int s = 0; s < SampleConfig::kNSamples; ++s) {
        std::vector<TH1F*>&       wpHists  = wpHistsAll[s];
        std::vector<std::string>& wpLabels = wpLabelsAll[s];
        std::vector<int>&         wpCols   = wpColsAll[s];

        if (wpHists.empty()) continue;

        DrawEffCanvas(wpHists, wpLabels, wpCols,
                      dispNames[s] + " Selection Efficiency",
                      outDir + "WPComparison_" + fileNames[s] + "_log.pdf",
                      true);
        DrawEffCanvas(wpHists, wpLabels, wpCols,
                      dispNames[s] + " Selection Efficiency",
                      outDir + "WPComparison_" + fileNames[s] + "_linear.pdf",
                      false);

        for (auto* h : wpHists) delete h;
    }

    Log::info(ctx, "Output: " + outDir);
}

}  // namespace EfficiencyImpl


// Entry point called by the Efficiency.cxx of both channels.
// variant: one variant key, or "ALL" for all variants of the channel.
inline void Efficiency_Run(const std::string& variant, const std::string& channel) {
    AtlasUtils::LoadAtlasStyle();
    MassCutConfig::PrintBanner("Efficiency");

    std::vector<std::string> keys = (variant == "ALL")
        ? GetPrimaryVariantKeys(channel)
        : std::vector<std::string>{variant};

    for (const auto& key : keys) {
        VariantConfig cfg = GetVariant(key);
        EfficiencyImpl::RunEfficiency(cfg);
    }
}

#endif // EFFICIENCY_H
