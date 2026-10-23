/**
 * ============================================================================
 *  MassFit_Plots.h - plots of the mass fit
 * ============================================================================
 * ===========================================
 * @author:    Thiziri AMEZZA
 * @email:     thiziri.amezza@cern.ch
 * @date:      18 April 2026
 * ===========================================
 *
 * Per-fit diagnostic plots:
 *   CreateIndividualFitAnalysis   for one background function: the S+B overlay
 *                                 (Mass_SB_Plot.h), the stacked Data/MC with the
 *                                 function, the extracted signal with its fit, and
 *                                 signal_fit_params.txt (read by
 *                                 CrossChecks/MassFit_BackgroundStability.C)
 *   DrawCVScoreVsComplexity       in-sample chi2/ndf and CV score against the
 *                                 order, one pad per family
 *   DrawCVFoldIllustration        the CV folds of one background function
 *   CreatePtBinComparisonPlots    collects the fit_summary.txt of the pT bins
 *
 * The per-fold CV plots are in MassFit_Selection.h, next to the CV itself.
 *
 * Part of the mass fit: it is not run on its own but through Mass_Fit.cxx
 * (how to run and all options: Mass_Fit.h). No option is read here.
 *
 * ==========================================================================
 * @version:  1.0
 * ==========================================================================
 */

#ifndef MASSFIT_PLOTS_H
#define MASSFIT_PLOTS_H

#include "MassFit_Types.h"
#include "Mass_SB_Plot.h"

inline void CreateIndividualFitAnalysis(FitFunction* fit, TH1F* h_data, TH1F* h_mc, TH1F* h_multijets_scaled,
                                        const std::string& outputPath, bool isOptimalFit = false,
                                        const std::string& ptLabel = "") {

    std::string fitDir = outputPath + fit->name + "/";
    gSystem->Exec(("mkdir -p " + fitDir).c_str());

    // Background subtraction to get signal
    TH1F* h_signal_extracted = (TH1F*)h_data->Clone((fit->name + "_signal").c_str());
    h_signal_extracted->SetDirectory(0);
    h_signal_extracted->SetTitle(("Isolated Signal using " + fit->name + ";M^{J} [GeV];Events").c_str());

    for (int i = 1; i <= h_signal_extracted->GetNbinsX(); ++i) {
        double binCenter = h_signal_extracted->GetBinCenter(i);
        double binContent = h_signal_extracted->GetBinContent(i);
        double fitValue = std::max(0.0, fit->func->Eval(binCenter));
        h_signal_extracted->SetBinContent(i, binContent - fitValue);
    }

    // Fit Gaussian to extracted signal - disable diagnostic plot here
    TF1* signal_gauss = PerformOptimizedGaussianFit(h_signal_extracted, (fit->name + "_gauss").c_str(),
                                                    {0, 0}, false, fitDir);

    // ----------------------------------------------------------------------
    // ATLAS-style Signal+Background overlay plot with pull panel.
    // Drawn for EVERY F-test-selected background function (one per family),
    // each into its own fitDir so they do not overlap. The dual-validation
    // winner is additionally flagged via isOptimalFit (adds a "(optimal)"
    // annotation and also writes a copy named signal_background_fit_OPTIMAL.pdf
    // for easy identification).
    {
        std::string topoLabel = (CurrentChannel() == "MJB") ? "Z(#rightarrow b#bar{b}) + jets"
                                                            : "Z(#rightarrow b#bar{b}) + #gamma";
        // Z(bb)+jets (MJB): signal is small on the multijet background, so it is
        // drawn x3 there for visibility (legend shows "Signal ×3"),
        // and zoom the Data/(S+B) ratio panel.
        bool isMJB = (CurrentChannel() == "MJB");
        double sigScale = isMJB ? 3.0 : 1.0;
        double ratioMin = isMJB ? 0.95 : 0.65;
        double ratioMax = isMJB ? 1.05 : 1.35;
        std::string ptShow = isMJB ? "p_{T} > 450 GeV" : "p_{T} > 200 GeV";
        (void)ptLabel;
        std::string variantLabel = CurrentVariantLabel();
        DrawSignalBackgroundPlot(h_data, fit->func, signal_gauss, h_signal_extracted, fitDir, topoLabel,
                                 RunConfig::LumiLabel(), isOptimalFit, ptShow, sigScale, ratioMin, ratioMax,
                                 variantLabel);

        // Cross-check: the same plot but with the ATLAS reference methodology
        if (h_mc)
            DrawSignalBackgroundPlot_MCTemplate(h_data, fit->func, h_mc, fitDir, topoLabel,
                                                RunConfig::LumiLabel(), isOptimalFit, ptShow, sigScale,
                                                ratioMin, ratioMax, variantLabel);
    }

    // Create stacked Data/MC with the background function.
    TCanvas* c_stack = new TCanvas("c_stack", ("Stack - " + fit->name).c_str(), 800, 800);
    auto pads = PadUtils::CreateRatioPads(c_stack, 0.30);
    pads.upper->cd();

    THStack* hs = new THStack("hs", (fit->name + " Background Model;M^{J} [GeV];Events").c_str());
    hs->Add(h_multijets_scaled);
    hs->Add(h_mc);
    hs->SetMaximum(hs->GetMaximum() * 1.25);
    hs->Draw("HIST");
    hs->GetYaxis()->SetTitle("Events");
    hs->GetYaxis()->SetTitleSize(0.05);
    hs->GetXaxis()->SetLabelSize(0); // Hide x-axis labels

    fit->func->SetLineColor(kGreen);
    fit->func->SetLineWidth(2);
    fit->func->Draw("same");
    HistUtils::StyleDataPoints(h_data);
    h_data->Draw("EP same");

    TH1F* h_total = (TH1F*)hs->GetStack()->Last();
    TGraphAsymmErrors* mc_err_band = HistUtils::MakeUncertaintyBand(h_total);
    mc_err_band->Draw("2 same");

    TLegend* leg = new TLegend(0.55, 0.55, 0.9, 0.85);
    leg->SetBorderSize(0);
    leg->AddEntry(h_data, "Data", "lep");
    leg->AddEntry(h_mc, SampleConfig::DisplayNames(CurrentChannel())[SampleConfig::kZbb].c_str(), "f");
    leg->AddEntry(h_multijets_scaled,
                  SampleConfig::DisplayNames(CurrentChannel())[SampleConfig::kMultijets].c_str(), "f");
    leg->AddEntry(mc_err_band, "MC Stat. Unc.", "f");

    double safe_chi2_ndf = (fit->ndf > 0 && std::isfinite(fit->chi2)) ? fit->chi2 / fit->ndf : 0.0;
    std::string fit_label = fit->name + " Fit";
    if (safe_chi2_ndf > 0 && safe_chi2_ndf < 1000) {
        char chi2_buffer[100];
        snprintf(chi2_buffer, sizeof(chi2_buffer), " (#chi^{2}/ndf = %.3f)", safe_chi2_ndf);
        fit_label += chi2_buffer;
    }
    leg->AddEntry(fit->func, fit_label.c_str(), "l");
    leg->Draw();

    AtlasUtils::DrawAtlasLabel();

    pads.lower->cd();
    pads.lower->SetGridy(0);
    TH1F* h_ratio = (TH1F*)h_data->Clone("h_ratio");
    h_ratio->SetDirectory(0);
    h_ratio->SetTitle("");
    h_ratio->Divide(h_total);
    PadUtils::StyleRatioHist(h_ratio, "M^{J} [GeV]", 0.5, 1.5);
    HistUtils::StyleDataPoints(h_ratio);
    h_ratio->Draw("EP");

    TGraphAsymmErrors* ratio_err_band = HistUtils::MakeRatioUncertaintyBand(h_total);
    ratio_err_band->Draw("2 same");
    h_ratio->Draw("EP same");

    TLine* line = PadUtils::DrawRefLine(pads.lower, 1.0);

    SaveCanvasPdf(c_stack, fitDir + "background_model.pdf");

    delete h_ratio;
    delete mc_err_band;
    delete ratio_err_band;
    delete leg;
    delete line;
    delete c_stack;

    // Signal plot
    TCanvas* c_signal = new TCanvas("c_signal", ("Signal - " + fit->name).c_str(), 800, 600);

    TH1F* h_mc_copy = (TH1F*)h_mc->Clone("h_mc_copy_temp");
    h_mc_copy->SetDirectory(0);
    h_mc_copy->SetFillColor(kRed);

    double mc_max = h_mc_copy->GetBinContent(h_mc_copy->GetMaximumBin());
    double signal_max = h_signal_extracted->GetBinContent(h_signal_extracted->GetMaximumBin());
    double overall_max = std::max(mc_max, signal_max);

    h_mc_copy->GetYaxis()->SetRangeUser(0, overall_max * 1.5);
    h_mc_copy->GetXaxis()->SetTitle("M^{J} [GeV]");
    h_mc_copy->GetYaxis()->SetTitle("Events");
    h_mc_copy->Draw("HIST");

    h_signal_extracted->SetMarkerStyle(21);
    h_signal_extracted->SetMarkerColor(kBlack);
    h_signal_extracted->Draw("EP same");

    signal_gauss->SetLineColor(kViolet - 4);
    signal_gauss->SetLineWidth(2);
    signal_gauss->Draw("same");

    // Fit MC
    TH1F* h_mc_for_fit = (TH1F*)h_mc->Clone("h_mc_for_fit_temp");
    h_mc_for_fit->SetDirectory(0);
    TF1* mc_gauss = PerformOptimizedGaussianFit(h_mc_for_fit, "mc_gauss_temp", {0, 0}, false, "");

    mc_gauss->SetLineColor(kCyan);
    mc_gauss->SetLineWidth(2);
    mc_gauss->Draw("same");

    TLegend* leg_sig = new TLegend(0.5, 0.7, 0.9, 0.9);
    leg_sig->SetBorderSize(0);
    leg_sig->SetTextSize(0.02);
    leg_sig->AddEntry(h_signal_extracted, "Isolated Signal (Data)", "lep");
    leg_sig->AddEntry(h_mc_copy, "MC Signal", "f");

    double sig_mean = signal_gauss->GetParameter(1);
    double sig_sigma = signal_gauss->GetParameter(2);
    double sig_chi2_ndf =
        (signal_gauss->GetNDF() > 0) ? signal_gauss->GetChisquare() / signal_gauss->GetNDF() : 0.0;

    if (std::isfinite(sig_mean) && std::isfinite(sig_sigma) && sig_mean > 0 && sig_sigma > 0) {
        leg_sig->AddEntry(
            signal_gauss,
            Form("Data: #mu=%.2f, #sigma=%.2f, #chi^{2}/ndf=%.2f", sig_mean, sig_sigma, sig_chi2_ndf), "l");
    }

    double mc_mean = mc_gauss->GetParameter(1);
    double mc_sigma = mc_gauss->GetParameter(2);
    double mc_chi2_ndf = (mc_gauss->GetNDF() > 0) ? mc_gauss->GetChisquare() / mc_gauss->GetNDF() : 0.0;

    if (std::isfinite(mc_mean) && std::isfinite(mc_sigma) && mc_mean > 0 && mc_sigma > 0) {
        leg_sig->AddEntry(
            mc_gauss, Form("MC: #mu=%.2f, #sigma=%.2f, #chi^{2}/ndf=%.2f", mc_mean, mc_sigma, mc_chi2_ndf),
            "l");
    }
    leg_sig->Draw();

    AtlasUtils::DrawAtlasLabel(0.2, 0.9);
    SaveCanvasPdf(c_signal, fitDir + "extracted_signal.pdf");

    // --- Background-modeling stability: ne file per background function results subfolder;
    //     The MassFit_BackgroundStability.C reads the three functions' files per
    //     (variant, WP, pT bin) and reports the spread of the data peak (mu)
    //     and width (sigma) as the background-modeling systematic for stability tests.
    //     NOTE: mc_gauss is fit to the MC signal template directly (no background
    //     subtraction), so mc_mu/mc_sigma are identical across functions-the function
    //     dependence lives entirely in data_mu/data_sigma (how well we model the background).
    {
        std::ofstream pf((fitDir + "signal_fit_params.txt").c_str());
        if (pf.is_open()) {
            pf << "# Per-function signal-fit results (background-modeling stability)\n";
            pf << "function " << fit->name << "\n";
            pf << "type " << fit->type << "\n";
            pf << "order " << fit->order << "\n";
            pf << "npar " << fit->npar << "\n";
            pf << "is_optimal " << (isOptimalFit ? 1 : 0) << "\n";
            pf << "bkg_chi2_ndf " << fit->chi2_ndf << "\n";
            pf << "spurious_signal " << fit->spurious_signal << "\n";
            pf << "data_mu " << sig_mean << "\n";
            pf << "data_mu_err " << signal_gauss->GetParError(1) << "\n";
            pf << "data_sigma " << sig_sigma << "\n";
            pf << "data_sigma_err " << signal_gauss->GetParError(2) << "\n";
            pf << "data_chi2_ndf " << sig_chi2_ndf << "\n";
            pf << "mc_mu " << mc_mean << "\n";
            pf << "mc_mu_err " << mc_gauss->GetParError(1) << "\n";
            pf << "mc_sigma " << mc_sigma << "\n";
            pf << "mc_sigma_err " << mc_gauss->GetParError(2) << "\n";
            pf << "mc_chi2_ndf " << mc_chi2_ndf << "\n";
            pf << "JMS " << ((mc_mean != 0) ? sig_mean / mc_mean : 0.0) << "\n";
            pf << "JMR " << ((mc_sigma != 0) ? sig_sigma / mc_sigma : 0.0) << "\n";
            pf.close();
        }
    }

    delete leg_sig;
    delete c_signal;
    delete h_mc_copy;
    delete h_mc_for_fit;
    delete mc_gauss;
    delete signal_gauss;
    delete h_signal_extracted;

    // Summary
    std::ofstream summary((fitDir + "summary.txt").c_str());
    if (summary.is_open()) {
        summary << "Analysis using " << fit->name << " as background model" << std::endl;
        summary << "==========================================" << std::endl << std::endl;
        summary << "Background Fit Quality:" << std::endl;
        summary << "  Chi^2/ndf = " << fit->chi2_ndf << std::endl;
        summary << "  Spurious signal = " << fit->spurious_signal << " events" << std::endl;
        summary.close();
    }

    std::cout << "  -> Saved to " << fitDir << std::endl;
}

// ============================================================================
// CV DIAGNOSTIC PLOTS
// ----------------------------------------------------------------------------
// (1) DrawCVScoreVsComplexity:
//     One pad per family: in-sample chi2/ndf (falls with the order) vs the
//     vs out-of-sample CV Chi^2/bin, with the 1-SE band and a star on the chosen
//     function (rises again when the function overfits).
// (2) DrawCVFoldIllustration:
//     A single pad that shows what one CV fold does: the data in the
//     training bins (filled) vs held-out bins (open) and the curve trained on
//     the training bins, with the held-out Chi^2.
// ============================================================================
inline void DrawCVScoreVsComplexity(const std::vector<FitFunction>& allFits, const FitFunction* bestFit,
                                    double cvMin, double cvThreshold, const std::string& outDir,
                                    const std::string& titleTag) {
    const std::vector<std::string> fams = {"PolyExp", "Poly", "Exp"};
    TCanvas* c = new TCanvas("c_cv_complexity", "CV vs complexity", 1500, 520);
    c->Divide(3, 1);
    int pad = 1;
    for (const auto& fam : fams) {
        c->cd(pad++);
        gPad->SetLeftMargin(0.16);
        gPad->SetBottomMargin(0.13);

        std::vector<const FitFunction*> fs;
        for (const auto& f : allFits)
            if (f.type == fam) fs.push_back(&f);
        std::sort(fs.begin(), fs.end(),
                  [](const FitFunction* a, const FitFunction* b) { return a->order < b->order; });

        TGraph* gIn = new TGraph();             // in-sample Chi^2/ndf
        TGraphErrors* gCV = new TGraphErrors(); // out-of-sample CV Chi^2/bin +/- SE
        int pIn = 0, pCV = 0;
        double ymin = 1e30, ymax = -1e30;
        for (const auto* f : fs) {
            if (f->chi2_ndf < 999) {
                gIn->SetPoint(pIn++, f->order, f->chi2_ndf);
                ymin = std::min(ymin, f->chi2_ndf);
                ymax = std::max(ymax, f->chi2_ndf);
            }
            if (f->cv_valid) {
                gCV->SetPoint(pCV, f->order, f->cv_score);
                gCV->SetPointError(pCV, 0.0, f->cv_sem);
                pCV++;
                // Range must include the CV points + their SE, else the
                // CV curve runs off-scale (one family can have large CV).
                ymin = std::min(ymin, f->cv_score - f->cv_sem);
                ymax = std::max(ymax, f->cv_score + f->cv_sem);
            }
        }
        // Always keep the acceptance band edges (best CV .. best CV + 1 SE) in view.
        ymin = std::min(ymin, cvMin);
        ymax = std::max(ymax, cvThreshold);
        if (ymin > ymax) {
            ymin = 0.0;
            ymax = 2.0;
        }
        double span = (ymax - ymin > 1e-9) ? (ymax - ymin) : 1.0;
        double pad_lo = std::max(0.0, ymin - 0.10 * span);
        double pad_hi = ymax + 0.42 * span; // headroom for the family title + legend

        gIn->SetMarkerStyle(24);
        gIn->SetMarkerColor(kGray + 2);
        gIn->SetLineColor(kGray + 2);
        gIn->SetLineStyle(2);
        gIn->SetMarkerSize(1.3);
        gCV->SetMarkerStyle(20);
        gCV->SetMarkerColor(kAzure + 2);
        gCV->SetLineColor(kAzure + 2);
        gCV->SetLineWidth(2);
        gCV->SetMarkerSize(1.3);

        TH1F* fr = new TH1F(Form("cvfr_%s", fam.c_str()), "", 10, 1.5, 5.5);
        fr->SetStats(0);
        fr->SetMinimum(pad_lo);
        fr->SetMaximum(pad_hi);
        fr->GetXaxis()->SetTitle("Polynomial order N");
        fr->GetYaxis()->SetTitle("in-sample #chi^{2}/ndf   &   CV #chi^{2}/bin");
        fr->GetXaxis()->SetTitleSize(0.050);
        fr->GetXaxis()->SetLabelSize(0.045);
        fr->GetXaxis()->SetTitleOffset(1.15);
        fr->GetYaxis()->SetTitleSize(0.046);
        fr->GetYaxis()->SetLabelSize(0.045);
        fr->GetYaxis()->SetTitleOffset(1.60);
        fr->GetXaxis()->SetNdivisions(404);
        fr->Draw();

        // Acceptance band [best CV, best CV + 1 SE] (green), with dashed upper edge.
        TBox* band = new TBox(1.5, cvMin, 5.5, cvThreshold);
        band->SetFillColorAlpha(kGreen + 2, 0.15);
        band->SetLineWidth(0);
        band->Draw("same");
        TLine* lthr = new TLine(1.5, cvThreshold, 5.5, cvThreshold);
        lthr->SetLineColor(kGreen + 2);
        lthr->SetLineStyle(2);
        lthr->SetLineWidth(2);
        lthr->Draw("same");

        if (pIn > 0) gIn->Draw("LP same");
        if (pCV > 0) gCV->Draw("LP same");

        const bool isWinnerPad = (bestFit && bestFit->type == fam && bestFit->cv_valid);
        TGraph* gStar = nullptr;
        if (isWinnerPad) {
            gStar = new TGraph();
            gStar->SetPoint(0, bestFit->order, bestFit->cv_score);
            gStar->SetMarkerStyle(29);
            gStar->SetMarkerSize(2.8);
            gStar->SetMarkerColor(kRed + 1);
            gStar->Draw("P same");
        }

        // Family title (inside the frame, top-left) so each panel is self-identifying.
        TLatex tf;
        tf.SetNDC();
        tf.SetTextFont(62);
        tf.SetTextSize(0.042);
        tf.SetTextColor(kBlack);
        tf.DrawLatex(0.20, 0.89, (fam + " family").c_str());
        // Winning (optimal) family only: note the selected model + its CV compactly, under
        // the title (the red star + legend already flag which point it is).
        if (isWinnerPad) {
            TLatex tsv;
            tsv.SetNDC();
            tsv.SetTextColor(kRed + 1);
            tsv.SetTextSize(0.026);
            tsv.DrawLatex(0.20, 0.83, Form("selected: %s", bestFit->name.c_str()));
            tsv.DrawLatex(0.20, 0.785, Form("CV = %.2f #pm %.2f", bestFit->cv_score, bestFit->cv_sem));
        }

        // Compact legend, top-right. "in-sample"/"CV" already imply training/held-out.
        TLegend* leg = new TLegend(0.52, 0.69, 0.95, 0.90);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.027);
        leg->AddEntry(gIn, "in-sample #chi^{2}/ndf", "lp");
        leg->AddEntry(gCV, "CV #chi^{2}/bin #pm SE", "lp");
        leg->AddEntry(band, "best CV + 1 SE band", "f");
        leg->AddEntry(lthr, "1-SE threshold", "l");
        if (isWinnerPad) leg->AddEntry(gStar, "selected model", "p");
        leg->Draw();
    }
    // No global canvas-top text: each panel carries its own family title, and the
    // winning panel notes the selected model + CV inside it.
    (void)titleTag;
    c->SaveAs((outDir + "cv_score_vs_complexity.pdf").c_str());
    delete c;
}

inline void DrawCVFoldIllustration(TF1* templateFunc, TH1F* h, double xmin, double xmax, int k,
                                   const std::string& outDir, const std::string& fname,
                                   const std::string& titleTag) {
    if (!templateFunc || !h) return;
    std::vector<int> bins;
    for (int b = 1; b <= h->GetNbinsX(); ++b) {
        double xc = h->GetBinCenter(b);
        if (xc < xmin || xc > xmax) continue;
        if (h->GetBinError(b) <= 0) continue;
        bins.push_back(b);
    }
    const int N = (int)bins.size();
    if (N < k) return;

    TGraphErrors* gTrain = new TGraphErrors();
    TGraphErrors* gTest = new TGraphErrors();
    TGraphErrors grFit;
    int pTr = 0, pTe = 0, pf = 0;
    double ymax = 0;
    for (int idx = 0; idx < N; ++idx) {
        int b = bins[idx];
        double x = h->GetBinCenter(b), y = h->GetBinContent(b), e = h->GetBinError(b);
        ymax = std::max(ymax, y + e);
        if (idx % k == 0) {
            gTest->SetPoint(pTe, x, y);
            gTest->SetPointError(pTe, 0, e);
            pTe++;
        } else {
            gTrain->SetPoint(pTr, x, y);
            gTrain->SetPointError(pTr, 0, e);
            pTr++;
            grFit.SetPoint(pf, x, y);
            grFit.SetPointError(pf, 0, e);
            pf++;
        }
    }
    TF1* ff = (TF1*)templateFunc->Clone("cvfold_fit");
    for (int ip = 0; ip < ff->GetNpar(); ++ip) ff->SetParameter(ip, templateFunc->GetParameter(ip));
    grFit.Fit(ff, "Q0");

    double chi2t = 0;
    int nt = 0;
    for (int idx = 0; idx < N; ++idx) {
        if (idx % k != 0) continue;
        int b = bins[idx];
        double e = h->GetBinError(b);
        double r = (h->GetBinContent(b) - ff->Eval(h->GetBinCenter(b))) / e;
        chi2t += r * r;
        ++nt;
    }

    TCanvas* c = new TCanvas("c_cv_fold", "CV fold illustration", 800, 600);
    c->SetLeftMargin(0.16);
    c->SetBottomMargin(0.14);
    c->SetRightMargin(0.05);
    c->SetTopMargin(0.06);
    TH1F* fr = new TH1F("cvfold_fr", "", 10, xmin, xmax);
    fr->SetStats(0);
    fr->SetMinimum(0.0);
    fr->SetMaximum(1.30 * ymax);
    fr->GetXaxis()->SetTitle("m [GeV]");
    fr->GetYaxis()->SetTitle("Background template (events)");
    fr->GetXaxis()->SetTitleSize(0.048);
    fr->GetXaxis()->SetTitleOffset(1.15);
    fr->GetXaxis()->SetLabelSize(0.042);
    fr->GetYaxis()->SetTitleSize(0.048);
    fr->GetYaxis()->SetTitleOffset(1.55);
    fr->GetYaxis()->SetLabelSize(0.042);
    fr->Draw();
    gTrain->SetMarkerStyle(20);
    gTrain->SetMarkerColor(kBlack);
    gTrain->SetMarkerSize(1.0);
    gTrain->SetLineColor(kBlack);
    gTest->SetMarkerStyle(24);
    gTest->SetMarkerColor(kRed + 1);
    gTest->SetMarkerSize(1.3);
    gTest->SetLineColor(kRed + 1);
    ff->SetLineColor(kAzure + 2);
    ff->SetLineWidth(2);
    ff->Draw("L same");
    gTrain->Draw("P same");
    gTest->Draw("P same");

    TLegend* leg = new TLegend(0.42, 0.67, 0.94, 0.87);
    leg->SetBorderSize(1);
    leg->SetFillStyle(1001);
    leg->SetFillColor(kWhite);
    leg->SetTextSize(0.028);
    leg->SetHeader("Illustrative fold (1 of k held out)");
    leg->AddEntry(gTrain, "Training bins (fitted)", "p");
    leg->AddEntry(gTest, "Held-out bins (scored, not fitted)", "p");
    leg->AddEntry(ff, Form("%s fit to training bins", templateFunc->GetName()), "l");
    leg->Draw();
    TLatex tx;
    tx.SetNDC();
    tx.SetTextSize(0.036);
    tx.DrawLatex(0.18, 0.88,
                 Form("held-out #chi^{2}/bin = %.3f  (%d bins)", (nt > 0 ? chi2t / nt : 0.0), nt));
    (void)titleTag;

    c->SaveAs((outDir + fname).c_str());
    delete ff;
    delete c;
}

// ============================================================================
// SUMMARY FUNCTION TO COMPARE PT BINS
// ============================================================================
inline void CreatePtBinComparisonPlots(const std::string& dataset, const std::string& wp,
                                       const PtBinConfig::PtBinDef& ptCfg) {
    VariantConfig vcfg_comp = GetVariant(dataset);
    DatasetConfig config = DatasetConfig::FromVariant(vcfg_comp);

    std::string fileWP = PathBuilder::WPtoFileWP(wp);
    std::string basePath = (RunConfig::PlotDir() + "/") + config.output_subdir + "/Mass_Fit_"
                           + CurrentChannel() + "/" + fileWP + "/";

    // Create comparison canvas
    TCanvas* c_compare = new TCanvas("c_pt_comparison", "pT Bin Comparison", 1800, 600);
    c_compare->Divide(3, 1);

    std::cout << "\n=========================================" << std::endl;
    std::cout << "Creating pT bin comparison plots for " << dataset << " WP " << wp << std::endl;
    std::cout << "=========================================\n" << std::endl;

    // Read summaries and create comparison
    std::ofstream combined_summary((basePath + "pt_bin_comparison_summary.txt").c_str());
    if (combined_summary.is_open()) {
        combined_summary << "PT BIN COMPARISON SUMMARY" << std::endl;
        combined_summary << "Dataset: " << dataset << ", WP: " << wp << std::endl;
        combined_summary << "=================================" << std::endl << std::endl;

        for (int iPt = 0; iPt < ptCfg.nBins; ++iPt) {
            std::string ptBinPath = basePath + ptCfg.fileNames[iPt] + "/";
            std::string summaryFile = ptBinPath + "fit_summary.txt";

            combined_summary << "pT Bin " << (iPt + 1) << ": " << ptCfg.labels[iPt] << std::endl;
            combined_summary << "---------------------------------" << std::endl;

            std::ifstream infile(summaryFile.c_str());
            if (infile.is_open()) {
                std::string line;
                while (std::getline(infile, line)) {
                    combined_summary << "  " << line << std::endl;
                }
                infile.close();
            } else {
                combined_summary << "  No summary available" << std::endl;
            }
            combined_summary << std::endl;
        }

        combined_summary.close();
        std::cout << "Combined summary saved to: " << basePath << "pt_bin_comparison_summary.txt"
                  << std::endl;
    }

    SaveCanvasPdf(c_compare, basePath + "pt_bin_comparison.pdf");
    delete c_compare;
}

#endif // MASSFIT_PLOTS_H
