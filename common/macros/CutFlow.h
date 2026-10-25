/**
 * ============================================================================
 *  CutFlow.h - cutflow plots and tables (DB and MJB)
 * ============================================================================
 * ===========================================
 * @author:    Thiziri AMEZZA
 * @email:     thiziri.amezza@cern.ch
 * @date:      18 April 2026
 * ===========================================
 *
 * Reads the cutflow histogram "h_cutflow" written by the slimming producers
 * for data, the Z(bb) signal MC and the main background MC (gamma+jets for DB,
 * multijets for MJB), and makes, for each variant and working point:
 *
 *   1. SignalToBackgroundRatio.pdf      S/B after the last cut, vs WP
 *   2. DataMCRatio.pdf                  data/MC after the last cut, vs WP
 *   3. StackedCutflow_<wp>.pdf          stacked MC + data, with data/MC pad
 *      CutflowTable_<wp>.txt            events per cut, summary and
 *                                       relative efficiencies
 *
 * The efficiency plots (with error bars) are made by Efficiency.h.
 *
 * All outputs go to <PlotDir>/<variant>/CutFlow/.
 *
 * Differences between the channels:
 *
 *   Normalisation of the background MC (NormBinForChannel)
 *     DB:   gamma+jets is scaled to data at the "Veto" cut (bin 10), the last
 *           cut before the b-tagging.
 *     MJB:  no scaling (scale factor = 1). The multijet MC is already weighted
 *           with cross section x luminosity, and data/MC agrees within almost
 *           after the b-tagging.
 *     To scale MJB at a given cut, make NormBinForChannel return that bin.
 *
 *   Cuts: DB has 12, MJB has 8 (CutFlowUtils::DefaultCutNames).
 *
 *   Sample labels (gamma+jets / multijets, Z(bb)+gamma / Z(bb)+jets) come
 *   from SampleConfig::DisplayNames(cfg). The table columns are sized for the
 *   longer MJB labels, so both channels print the same table.
 *
 * The cutflows are made by the producers, so EASYBJETS_MASS_CUT does not
 * appear in them (a message says so when it is on).
 *
 * How to run (from ZbbGamma/macros or ZbbJets/macros):
 *
 *   export EASYBJETS_RUN=run3                      run2 or run3
 *   export EASYBJETS_DATA_DIR=$HOME/lpnatlas_data/ZbbGamma    (ZbbJets for MJB)
 *
 *   root -l -b -q 'CutFlow.cxx()'                  nominal (default)
 *   root -l -b -q 'CutFlow.cxx("DB_bJER10v00")'    one variant
 *   root -l -b -q 'CutFlow.cxx("ALL")'             all variants
 *   (MJB: "MJB_bJER10v00", ...)
 *
 * ==========================================================================
 * @version:  2.1
 * ==========================================================================
 */

#ifndef CUTFLOW_H
#define CUTFLOW_H

#include "VariantConfig.h"
#include "CommonUtils.h"
#include "MassCutConfig.h"
#include <fstream>

namespace CutFlowImpl {

// ============================================================================
// Bin where the background MC is scaled to data.
//   -1      no scaling, the MC weights are kept (MJB)
//   N >= 1  the background MC is scaled so that data/MC = 1 in bin N (DB)
// Why each channel differs: see the header of this file.
// Note: CutFlowUtils::NormBin in CommonUtils.h returns 1 for MJB and is not
// used here; this function is the one that sets the cutflow normalisation.
// ============================================================================
inline int NormBinForChannel(const std::string& channel) {
    if (channel == "MJB") return -1;   // no normalization (MJB design)
    return 10;                         // DB: "Veto" bin (just before b-tag)
}

// ============================================================================
// Scale factor data / background MC in the normalisation bin.
// Returns 1 if the channel has no normalisation or an input is missing.
// ============================================================================
inline double ComputeNormSF(TH1F* dataCF, TH1F* bkgCF, const std::string& channel) {
    if (!dataCF || !bkgCF) return 1.0;
    const int normBin = NormBinForChannel(channel);
    if (normBin < 1) return 1.0;  // no normalization
    if (normBin > dataCF->GetNbinsX() || normBin > bkgCF->GetNbinsX()) return 1.0;
    double dv = dataCF->GetBinContent(normBin);
    double mv = bkgCF ->GetBinContent(normBin);
    return (mv > 0) ? dv / mv : 1.0;
}

// Note added to the log line when the scaling is off.
inline std::string NormSFNote(const std::string& channel) {
    return (NormBinForChannel(channel) < 1) ? "  [disabled, SF=1]" : "";
}

// Writes the cut names as x-axis bin labels.
inline void SetCutLabels(TH1F* h, const std::vector<std::string>& cutNames) {
    if (!h) return;
    for (int bin = 1; bin <= (int)cutNames.size() && bin <= h->GetNbinsX(); ++bin)
        h->GetXaxis()->SetBinLabel(bin, cutNames[bin-1].c_str());
}


// ============================================================================
// 1 and 2. S/B and data/MC after the last cut, as a function of the WP
//   S/B      Z(bb) signal / (total MC - signal)
//   data/MC  data / total MC, error from the data and MC statistics:
//            (data/MC) * sqrt(1/data + (dMC/MC)^2), dMC from the bin errors
// The background MC is first scaled with ComputeNormSF.
// ============================================================================
inline void SBandDataMCvsWP(const VariantConfig& cfg) {
    auto workingPoints  = WPConfig::DefaultWPs();
    auto dispNames      = SampleConfig::DisplayNames(cfg);
    std::string outDir  = PathBuilder::OutputDir(cfg, "CutFlow");
    PathBuilder::MakeDir(outDir);

    const int nWP      = (int)workingPoints.size();
    const int finalBin = CutFlowUtils::NCuts(cfg.channel);

    TH1F* h_sb  = HistUtils::CreateDetached("h_sb_ratio",
        "S/B at Final Selection;Working Point;S/B", nWP, 0, nWP);
    TH1F* h_dmc = HistUtils::CreateDetached("h_dmc_ratio",
        "Data/MC at Final Selection;Working Point;Data/MC", nWP, 0, nWP);

    std::vector<std::string> wpLabels;
    std::vector<double>      sbVals;

    for (int i = 0; i < nWP; ++i) {
        const auto& wp     = workingPoints[i];
        std::string fileWP = PathBuilder::WPtoFileWP(wp);

        h_sb ->GetXaxis()->SetBinLabel(i+1, fileWP.c_str());
        h_dmc->GetXaxis()->SetBinLabel(i+1, fileWP.c_str());

        auto cfSet = CutFlowUtils::LoadAllCutflows(cfg, wp, "CutFlow");
        if (!cfSet.IsLoaded(SampleConfig::kData) || !cfSet.IsLoaded(SampleConfig::kMultijets))
            continue;

        TH1F* h_cf[SampleConfig::kNSamples];
        bool  loaded[SampleConfig::kNSamples];
        cfSet.GetRawArrays(h_cf, loaded);

        double sf = ComputeNormSF(cfSet[SampleConfig::kData],
                                   cfSet[SampleConfig::kMultijets],
                                   cfg.channel);
        cfSet[SampleConfig::kMultijets]->Scale(sf);
        Log::info("CutFlow", Form("  WP %s  %s cutflow norm SF: %.4f%s",
                                   fileWP.c_str(),
                                   dispNames[SampleConfig::kMultijets].c_str(),
                                   sf, NormSFNote(cfg.channel).c_str()));

        double sigVal  = cfSet.IsLoaded(SampleConfig::kZbb)
                         ? cfSet[SampleConfig::kZbb]->GetBinContent(finalBin) : 0;
        double totalMC = SampleConfig::TotalMC(h_cf, loaded, finalBin);
        double bkgVal  = totalMC - sigVal;

        double sbRatio  = (bkgVal  > 0) ? sigVal  / bkgVal  : 0;
        double dataVal  = cfSet[SampleConfig::kData]->GetBinContent(finalBin);
        double dmcRatio = (totalMC > 0) ? dataVal / totalMC : 0;
        // MC statistical error: sum in quadrature of the MC bin errors
        double mcErr2 = 0;
        for (int s = 0; s < SampleConfig::kNSamples; ++s)
            if (s != SampleConfig::kData && loaded[s] && h_cf[s])
                mcErr2 += h_cf[s]->GetBinError(finalBin) * h_cf[s]->GetBinError(finalBin);
        double ratioErr = 0;
        if (dataVal > 0 && totalMC > 0)
            ratioErr = dmcRatio * std::sqrt(1.0 / dataVal + mcErr2 / (totalMC * totalMC));

        h_sb ->SetBinContent(i+1, sbRatio);
        h_dmc->SetBinContent(i+1, dmcRatio);
        h_dmc->SetBinError  (i+1, ratioErr);

        wpLabels.push_back(fileWP);
        sbVals  .push_back(sbRatio);
    }

    {
        CanvasGuard cg("canvasSB", "S/B Ratio", 1000, 800,
                       outDir + "SignalToBackgroundRatio.pdf");
        cg->SetGridy();
        h_sb->SetMarkerStyle(20);
        h_sb->SetMarkerSize(1.5);
        h_sb->SetMarkerColor(kBlue);
        h_sb->SetLineColor(kBlue);
        h_sb->SetLineWidth(2);
        h_sb->Draw("EP");
        AtlasUtils::DrawAtlasLabel(0.2, 0.88);
    }

    std::cout << "\n+---------------+----------------+\n"
              << "| Working Point | S/B Ratio      |\n"
              << "+---------------+----------------+\n";
    for (size_t i = 0; i < wpLabels.size(); ++i)
        std::cout << "| " << std::setw(13) << wpLabels[i] << " | "
                  << std::right << std::setw(14) << std::fixed
                  << std::setprecision(4) << sbVals[i] << " |\n";
    std::cout << "+---------------+----------------+\n";

    {
        CanvasGuard cg("canvasDMC", "Data/MC Ratio", 900, 600,
                       outDir + "DataMCRatio.pdf");
        cg->SetGridy();
        h_dmc->SetMarkerStyle(21);
        h_dmc->SetMarkerColor(kBlue);
        h_dmc->SetLineColor(kBlue);
        h_dmc->SetLineWidth(2);
        h_dmc->GetYaxis()->SetRangeUser(0.5, 1.5);
        h_dmc->Draw("EP");
        TLine* line = new TLine(0, 1.0, nWP, 1.0);
        line->SetLineColor(kRed);
        line->SetLineStyle(2);
        line->SetLineWidth(2);
        line->Draw();
        AtlasUtils::DrawAtlasLabel(0.2, 0.88);
    }

    std::cout << "\n+---------------+------------------------+\n"
              << "| Working Point | Data/MC Ratio          |\n"
              << "+---------------+------------------------+\n";
    for (size_t i = 0; i < wpLabels.size(); ++i)
        std::cout << "| " << std::left  << std::setw(13) << wpLabels[i] << " | "
                  << std::right << std::setw(10) << std::fixed << std::setprecision(4)
                  << h_dmc->GetBinContent(i+1) << " +/- "
                  << std::setprecision(4) << h_dmc->GetBinError(i+1) << " |\n";
    std::cout << "+---------------+------------------------+\n";

    delete h_sb;
    delete h_dmc;
}


// ============================================================================
// 3. For each WP: stacked cutflow with a data/MC pad and the cutflow table
// ============================================================================
inline void StackedCutflows(const VariantConfig& cfg) {
    auto workingPoints = WPConfig::DefaultWPs();
    auto displayNames  = SampleConfig::DisplayNames(cfg);
    auto sampleColors  = SampleConfig::Colors();
    auto fillStyles    = SampleConfig::FillStyles();

    std::string outDir = PathBuilder::OutputDir(cfg, "CutFlow");
    PathBuilder::MakeDir(outDir);

    for (const auto& wp : workingPoints) {
        std::string fileWP  = PathBuilder::WPtoFileWP(wp);
        std::string labelWP = PathBuilder::WPtoLabel(wp);
        auto cutNames       = CutFlowUtils::DefaultCutNames(cfg.channel, labelWP);
        const int nCuts     = (int)cutNames.size();

        Log::info("CutFlow", "Stacked cutflow for WP " + fileWP);

        auto cfSet = CutFlowUtils::LoadAllCutflows(cfg, wp, "CutFlow");
        if (!cfSet.IsLoaded(SampleConfig::kData) || !cfSet.IsLoaded(SampleConfig::kMultijets)) {
            Log::warning("CutFlow", "Missing data or " +
                         displayNames[SampleConfig::kMultijets] + " for WP " + fileWP);
            continue;
        }

        TH1F* h_cf[SampleConfig::kNSamples];
        bool  loaded[SampleConfig::kNSamples];
        cfSet.GetRawArrays(h_cf, loaded);

        for (int s = 0; s < SampleConfig::kNSamples; ++s)
            if (loaded[s] && h_cf[s]) h_cf[s]->SetStats(0);

        double sf = ComputeNormSF(cfSet[SampleConfig::kData],
                                   cfSet[SampleConfig::kMultijets],
                                   cfg.channel);
        cfSet[SampleConfig::kMultijets]->Scale(sf);
        Log::info("CutFlow", Form("  %s cutflow norm SF: %.4f%s",
                                   displayNames[SampleConfig::kMultijets].c_str(),
                                   sf, NormSFNote(cfg.channel).c_str()));

        TH1F* hData = cfSet[SampleConfig::kData];
        TH1F* hZbb  = cfSet[SampleConfig::kZbb];
        TH1F* hMJ   = cfSet[SampleConfig::kMultijets];

        hData->SetMarkerStyle(20);
        hData->SetMarkerSize(1.2);
        hData->SetMarkerColor(kBlack);
        hData->SetLineColor(kBlack);
        hData->SetLineWidth(2);

        for (int s = 0; s < SampleConfig::kNSamples; ++s) {
            if (s == SampleConfig::kData || !loaded[s]) continue;
            h_cf[s]->SetFillColor(sampleColors[s]);
            h_cf[s]->SetFillStyle(fillStyles[s]);
            h_cf[s]->SetLineColor(sampleColors[s]);
            h_cf[s]->SetLineWidth(1);
        }

        for (int s = 0; s < SampleConfig::kNSamples; ++s)
            if (loaded[s]) SetCutLabels(h_cf[s], cutNames);

        // ----------------------------------------------------------------
        // Stacked cutflow canvas with Data/MC ratio panel
        // ----------------------------------------------------------------
        {
            CanvasGuard cg(("canvasStackedCF_" + fileWP),
                ("Stacked Cutflow (" + labelWP + ")"), 1000, 800,
                outDir + "StackedCutflow_" + fileWP + ".pdf");

            TPad* padTop = new TPad("padTop", "padTop", 0, 0.3, 1, 1);
            padTop->SetBottomMargin(0.02);
            padTop->SetLogy();
            padTop->Draw();
            padTop->cd();

            THStack* stack = new THStack(("stackCF_" + fileWP).c_str(),
                ("Cutflow (" + labelWP + ");Cut;Events").c_str());
            if (loaded[SampleConfig::kZbb])
                stack->Add(h_cf[SampleConfig::kZbb]);
            if (loaded[SampleConfig::kMultijets])
                stack->Add(h_cf[SampleConfig::kMultijets]);

            TH1F* hSumMC = SampleConfig::MakeTotalMC(h_cf, loaded, "h_sumMC_" + fileWP);

            double maxY = std::max(hData->GetMaximum(),
                       hSumMC ? hSumMC->GetMaximum() : 0.0);
            stack->Draw("HIST");
            stack->SetMaximum(maxY * 1000);
            stack->SetMinimum(0.1);
            stack->GetXaxis()->SetLabelSize(0);
            stack->GetYaxis()->SetTitle("Events");
            stack->GetYaxis()->SetTitleSize(0.05);

            // MC statistical uncertainty band
            TGraphAsymmErrors* mc_err_band = nullptr;
            if (hSumMC) {
                mc_err_band = HistUtils::MakeUncertaintyBand(hSumMC);
                mc_err_band->Draw("2 SAME");
            }

            HistUtils::StyleDataPoints(hData);
            hData->Draw("EP SAME");

            TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9);
            legend->SetBorderSize(0);
            legend->SetFillStyle(0);
            SampleConfig::AddLegendEntries(legend, h_cf, nullptr, loaded, cfg);
            if (mc_err_band) legend->AddEntry(mc_err_band, "MC Stat. Unc.", "f");
            legend->Draw();

            AtlasUtils::DrawAtlasLabel(0.2, 0.88);

            TLatex wpLbl;
            wpLbl.SetNDC();
            wpLbl.SetTextFont(42);
            wpLbl.SetTextSize(0.04);
            wpLbl.DrawLatex(0.2, 0.79, ("Working Point: " + labelWP).c_str());

            // Lower pad: data / MC
            cg.Canvas()->cd();
            TPad* padBottom = new TPad("padBottom", "padBottom", 0, 0.05, 1, 0.3);
            padBottom->SetTopMargin(0.05);
            padBottom->SetBottomMargin(0.3);
            padBottom->Draw();
            padBottom->cd();

            TH1F* hRatio = nullptr;
            if (hSumMC) {
                hRatio = static_cast<TH1F*>(hData->Clone(("h_ratio_" + fileWP).c_str()));
                hRatio->Divide(hSumMC);
                hRatio->SetTitle("");
                hRatio->GetYaxis()->SetTitle("Data/MC");
                hRatio->GetYaxis()->SetTitleSize(0.10);
                hRatio->GetYaxis()->SetTitleOffset(0.5);
                hRatio->GetYaxis()->SetLabelSize(0.08);
                hRatio->GetYaxis()->SetRangeUser(0.5, 1.5);
                hRatio->GetYaxis()->SetNdivisions(505);
                hRatio->GetXaxis()->SetTitle("Cuts");
                hRatio->GetXaxis()->SetTitleSize(0.10);
                hRatio->GetXaxis()->SetTitleOffset(1.4);
                hRatio->GetXaxis()->SetLabelSize(0.10);
                HistUtils::StyleDataPoints(hRatio);
                hRatio->Draw("EP");

                // Hatched MC uncertainty band around 1
                TGraphAsymmErrors* ratio_err_band = HistUtils::MakeRatioUncertaintyBand(hSumMC);
                ratio_err_band->Draw("2 SAME");

                TLine* line = new TLine(0, 1, nCuts, 1);
                line->SetLineStyle(2);
                line->SetLineColor(kBlack);
                line->SetLineWidth(2);
                line->Draw("SAME");

                hRatio->Draw("EP SAME");  // Draw points on top of band
            }
        }

        // ----------------------------------------------------------------
        // Cutflow table, printed and written to CutflowTable_<wp>.txt
        // ----------------------------------------------------------------
        std::ofstream tableFile(outDir + "CutflowTable_" + fileWP + ".txt");

        const std::string& zbbLabel = displayNames[SampleConfig::kZbb];
        const std::string& mjLabel  = displayNames[SampleConfig::kMultijets];

        auto writeHeader = [&](std::ostream& os) {
            os << "\nCutflow table for working point " << fileWP << ":\n"
               << "+------+----------------------+-----------------"
                  "+--------------------+---------------------+----------------+\n"
               << "| Cut  | Description          | Data            | "
               << std::left  << std::setw(20) << zbbLabel
               << "| " << std::setw(21) << mjLabel
               << "| Data/MC        |\n"
               << "+------+----------------------+-----------------"
                  "+--------------------+---------------------+----------------+\n";
        };
        writeHeader(std::cout);
        writeHeader(tableFile);

        for (int b = 1; b <= nCuts && b <= hData->GetNbinsX(); ++b) {
            double dVal    = hData->GetBinContent(b);
            double sVal    = loaded[SampleConfig::kZbb] ? hZbb->GetBinContent(b) : 0;
            double mVal    = hMJ->GetBinContent(b);
            double totalMC = SampleConfig::TotalMC(h_cf, loaded, b);
            double ratio   = (totalMC > 0) ? dVal / totalMC : 0;

            std::ostringstream line;
            line << "| " << std::setw(4) << b-1 << " | "
                 << std::left  << std::setw(20) << cutNames[b-1] << " | "
                 << std::right << std::setw(9)  << std::fixed << std::setprecision(1) << dVal << " | "
                 << std::setw(12) << std::setprecision(1) << sVal << " | "
                 << std::setw(13) << std::setprecision(1) << mVal << " | "
                 << std::setw(8)  << std::setprecision(3) << ratio << " |";

            std::cout  << line.str() << "\n";
            tableFile  << line.str() << "\n";
        }

        const std::string sep =
            "+------+----------------------+-----------------"
            "+--------------------+---------------------+----------------+\n";
        std::cout << sep;
        tableFile << sep;

        int    fb        = nCuts;
        double finalData = hData->GetBinContent(fb);
        double finalSig  = loaded[SampleConfig::kZbb] ? hZbb->GetBinContent(fb) : 0;
        double finalMJ   = hMJ->GetBinContent(fb);
        double finalTMC  = SampleConfig::TotalMC(h_cf, loaded, fb);
        double finalDMC  = (finalTMC > 0) ? finalData / finalTMC : 0;
        // same S/B definition as SignalToBackgroundRatio.pdf: signal / (total MC - signal)
        double finalBkg  = finalTMC - finalSig;
        double finalSB   = (finalBkg > 0) ? finalSig / finalBkg : 0;
        double overallEff= (hData->GetBinContent(1) > 0)
                           ? finalData / hData->GetBinContent(1) : 0;

        std::ostringstream summary;
        summary << "\nSummary for working point " << fileWP << ":\n"
                << "  - Final data events:        " << finalData << "\n"
                << "  - Final " << zbbLabel << " events:  " << finalSig  << "\n"
                << "  - Final " << mjLabel << " events:    " << finalMJ   << "\n"
                << "  - Final data/MC ratio:      " << std::fixed << std::setprecision(3) << finalDMC << "\n"
                << "  - Final S/B ratio:          " << std::setprecision(4) << finalSB   << "\n"
                << "  - Overall efficiency:       " << std::setprecision(6) << overallEff << "\n";
        std::cout << summary.str();
        tableFile << summary.str();

        // Relative efficiency table (file only): each cut / "Initial"
        double initialData = hData->GetBinContent(1);
        double initialZbb  = loaded[SampleConfig::kZbb] ? hZbb->GetBinContent(1) : 0;
        double initialMJ   = hMJ->GetBinContent(1);

        tableFile << "\nRelative efficiencies (normalized to initial events):\n"
                  << "+------+----------------------+-----------------"
                     "+------------------+-----------------+\n"
                  << "| Cut  | Description          | Data            | "
                  << std::left  << std::setw(17) << zbbLabel
                  << "| " << std::setw(16) << mjLabel
                  << "|\n"
                  << "+------+----------------------+-----------------"
                     "+------------------+-----------------+\n";

        for (int b = 1; b <= nCuts && b <= hData->GetNbinsX(); ++b) {
            double dEff = (initialData > 0) ? hData->GetBinContent(b) / initialData : 0;
            double sEff = (initialZbb  > 0) ? hZbb ->GetBinContent(b) / initialZbb  : 0;
            double mEff = (initialMJ   > 0) ? hMJ  ->GetBinContent(b) / initialMJ   : 0;

            tableFile << "| " << std::setw(4) << b-1 << " | "
                      << std::left  << std::setw(20) << cutNames[b-1] << " | "
                      << std::right << std::setw(9)  << std::fixed
                      << std::setprecision(6) << dEff << " | "
                      << std::setw(9) << sEff << " | "
                      << std::setw(9) << mEff << " |\n";
        }
        tableFile << "+------+----------------------+-----------------"
                     "+------------------+-----------------+\n";
        tableFile.close();

        Log::info("CutFlow", "Table saved: " + outDir + "CutflowTable_" + fileWP + ".txt");
    }
}

}  // namespace CutFlowImpl


// Entry point called by the CutFlow.cxx of both channels.
// variant: one variant key, or "ALL" for all variants of the channel.
inline void CutFlow_Run(const std::string& variant, const std::string& channel) {
    AtlasUtils::LoadAtlasStyle(); 
    MassCutConfig::PrintBanner("CutFlow");
    if (MassCutConfig::IsActive()) {
        std::cout << "[MassCutConfig] CutFlow displays slimmer-side cutflow "
                     "histograms; the post-slimmer EASYBJETS_MASS_CUT is NOT "
                     "shown here as an extra row (cutflows are pre-baked).\n";
    }

    auto runForVariant = [](const std::string& key) {
        VariantConfig cfg = GetVariant(key);
        PathBuilder::MakeDir(PathBuilder::OutputDir(cfg, "CutFlow"));
        Log::info("CutFlow", "=== " + cfg.label + " (" + cfg.key + ") ===");
        CutFlowImpl::SBandDataMCvsWP(cfg);
        CutFlowImpl::StackedCutflows(cfg);
    };

    if (variant == "ALL") {
        for (const auto& key : GetPrimaryVariantKeys(channel)) {
            Log::info("CutFlow", "Running for: " + key);
            runForVariant(key);
        }
    } else {
        runForVariant(variant);
    }

    Log::info("CutFlow", "All analyses complete.");
}

#endif // CUTFLOW_H
