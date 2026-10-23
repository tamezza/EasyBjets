/**
 * ============================================================================
 *  Mass_SB_Plot.h - signal + background overlay on the data mass spectrum
 * ============================================================================
 * ===========================================
 * @author:    Thiziri AMEZZA
 * @email:     thiziri.amezza@cern.ch
 * @date:      18 April 2026
 * ===========================================
 *
 * Signal + background plot of the large-R jet mass spectrum, in the usual ATLAS
 * style: the full (not subtracted) data with three curves on top and a
 * Data/(S+B) ratio pad below.
 *   background (blue, dashed)   the background function chosen by the F-test and
 *                               k-fold CV, fitted on the sidebands and
 *                               extrapolated under the peak
 *   signal     (orange, dashed) the Z->bb peak fitted on data minus background,
 *                               normalised to the extracted signal yield
 *   S+B        (red, solid)     the sum of the two
 *
 * DrawSignalBackgroundPlot() does not refit anything. It is called from
 * MassFit_Plots.h with the background and signal fits already in memory, so the
 * plot always shows exactly the model the mu and sigma were extracted from.
 *
 * DrawSignalBackgroundPlot_MCTemplate() is a cross-check with the standard ATLAS
 * method: the signal shape is the MC template, and its yield and the background
 * parameters are fitted together to the data over the full range.
 * Output: signal_background_fit_MCtemplate.{pdf,png} next to the first plot.
 *
 * When isOptimal is set, both also write an *_OPTIMAL copy, so the plot of the
 * selected background function is easy to find among all the candidates.
 *
 * Part of the mass fit: it is not run on its own but through Mass_Fit.cxx
 * (how to run and all options: Mass_Fit.h). No option is read here.
 *
 * ==========================================================================
 * @version:  1.0
 * ==========================================================================
 */

#ifndef MASS_SB_PLOT_H
#define MASS_SB_PLOT_H

#include "TH1F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TGraph.h"
#include "TStyle.h"
#include <string>
#include <cmath>

#include "RunEra.h"
#include <vector>

// Gaussian with the amplitude set so its bin sum over +-3 sigma equals
// the yield of h_signal_extracted in the same window.
inline TF1* BuildSignalGaussianForOverlay(TH1F* h_signal_extracted, double mu, double sigma, double xmin,
                                          double xmax, const std::string& name) {
    double lo = mu - 3.0 * sigma;
    double hi = mu + 3.0 * sigma;
    double yield = 0.0;
    double binW = h_signal_extracted->GetBinWidth(1);
    for (int i = 1; i <= h_signal_extracted->GetNbinsX(); ++i) {
        double c = h_signal_extracted->GetBinCenter(i);
        if (c < lo || c > hi) continue;
        yield += h_signal_extracted->GetBinContent(i);
    }
    if (yield < 0) yield = 0;

    // A normalised Gaussian integrates to (A * sigma * sqrt(2 pi)) in "x units".
    // The histogram is in events-per-bin; to match the data drawn as events/bin
    // we want G(m) to give the events in the bin at m, i.e. yield is distributed
    // as a Gaussian in m with the bin width folded in:
    //   G(m) = (yield * binW) / (sigma sqrt(2 pi)) * exp(-(m-mu)^2 / 2 sigma^2)
    // so that sum_i G(m_i) ~= yield over the peak.
    double amp = (sigma > 0) ? (yield * binW) / (sigma * std::sqrt(2.0 * M_PI)) : 0.0;

    TF1* g = new TF1(name.c_str(), "[0]*exp(-0.5*((x-[1])/[2])**2)", xmin, xmax);
    g->SetParameters(amp, mu, sigma);
    return g;
}

// Build the signal overlay from the actual fitted model (Gaussian/double-sided CB..),
// so the orange signal curve AND the red S+B total reflect the model that was really
// fitted - not a Gaussian reconstructed. Falls back to the Gaussian overlay
// if the fitted model is missing or its integral is degenerate.
inline TF1* BuildSignalModelForOverlay(TF1* signal_fit, TH1F* h_signal_extracted, double mu, double sigma,
                                       double xmin, double xmax, const std::string& name) {
    if (!signal_fit || signal_fit->GetNpar() <= 3)
        return BuildSignalGaussianForOverlay(h_signal_extracted, mu, sigma, xmin, xmax, name);

    TF1* f = (TF1*)signal_fit->Clone(name.c_str());
    if (!f) return BuildSignalGaussianForOverlay(h_signal_extracted, mu, sigma, xmin, xmax, name);
    f->SetRange(xmin, xmax);

    const double lo = mu - 3.0 * sigma, hi = mu + 3.0 * sigma;
    double yield = 0.0, sumModel = 0.0;
    for (int i = 1; i <= h_signal_extracted->GetNbinsX(); ++i) {
        double c = h_signal_extracted->GetBinCenter(i);
        if (c < lo || c > hi) continue;
        yield += h_signal_extracted->GetBinContent(i);
        sumModel += f->Eval(c);
    }
    if (yield < 0) yield = 0;
    const double amp0 = f->GetParameter(0);
    if (sumModel > 1e-12 && std::fabs(amp0) > 0) {
        f->SetParameter(0, amp0 * yield / sumModel); // Σ_bins f == yield
        return f;
    }
    delete f; // degenerate -> fall back to a plain Gaussian overlay
    return BuildSignalGaussianForOverlay(h_signal_extracted, mu, sigma, xmin, xmax, name);
}

// Data-driven S+B plot: background from the sideband fit, signal from the
// fit to data minus background.
//
//   h_data              full (not subtracted) data mass histogram
//   f_background        selected background function (sideband fit)
//   signal_gauss        signal fit to data minus background; any of the
//                       supported shapes (Gauss, Bukin, CB, DSCB), not only
//                       a Gaussian. Gives mu, sigma and the shape.
//   h_signal_extracted  data minus background, used to normalise the signal
//   outDir              where the pdf/png are written
//   topologyLabel       process label, e.g. "Z(#rightarrow b#bar{b}) + #gamma"
//   lumiLabel           sqrt(s) and luminosity, default from RunEra.h
//   isOptimal           also write an *_OPTIMAL copy (selected background)
//   ptLabel             pT bin text on the plot, empty for inclusive
//   signalScale         draws the signal xN to make it visible (e.g. x3 for
//                       Z+jets); the S+B curve and the ratio use the true signal
//   ratioMin, ratioMax  y range of the Data/(S+B) pad
//   variantLabel        calibration name on the plot (nominal, bJER10v01, ...)
inline void DrawSignalBackgroundPlot(TH1F* h_data, TF1* f_background, TF1* signal_gauss,
                                     TH1F* h_signal_extracted, const std::string& outDir,
                                     const std::string& topologyLabel,
                                     const std::string& lumiLabel = RunConfig::LumiLabel(),
                                     bool isOptimal = false, const std::string& ptLabel = "",
                                     double signalScale = 1.0, double ratioMin = 0.5, double ratioMax = 1.5,
                                     const std::string& variantLabel = "") {
    if (!h_data || !f_background || !signal_gauss || !h_signal_extracted) {
        std::cerr << "[Mass_SB_Plot] WARN: null input, skipping S+B plot.\n";
        return;
    }

    const double xmin = h_data->GetXaxis()->GetXmin();
    const double xmax = h_data->GetXaxis()->GetXmax();
    const double mu = signal_gauss->GetParameter(1);
    const double sigma = signal_gauss->GetParameter(2);

    // --- The signal overlay from the actual fitted model (orange) ---
    TF1* f_signal = BuildSignalModelForOverlay(signal_gauss, h_signal_extracted, mu, sigma, xmin, xmax,
                                               "f_signal_overlay");

    // --- The total S+B curve (red) as a fine-sampled TGraph ---
    const int nSample = 500;
    std::vector<double> gx(nSample), gy(nSample);
    for (int i = 0; i < nSample; ++i) {
        double x = xmin + (xmax - xmin) * (i + 0.5) / nSample;
        gx[i] = x;
        double bg = std::max(0.0, f_background->Eval(x));
        double sg = f_signal->Eval(x);
        gy[i] = bg + sg;
    }
    TGraph* g_total = new TGraph(nSample, gx.data(), gy.data());

    // Signal-overlay curve, optionally enlarged xsignalScale for visibility
    // (e.g. x3 for Z(bb)+jets, where the signal is tiny on the multijet
    // background).
    TF1* f_signal_draw = f_signal;
    if (signalScale != 1.0) {
        f_signal_draw = (TF1*)f_signal->Clone("f_signal_overlay_scaled");
        f_signal_draw->SetParameter(0, f_signal->GetParameter(0) * signalScale);
    }

    // ======================================================================
    // Canvas: main panel (top) + Data/MC ratio panel (bottom)
    // ======================================================================
    TCanvas* c = new TCanvas("c_sb_plot", "Signal + Background fit", 700, 750);

    TPad* pad1 = new TPad("pad_main", "pad_main", 0.0, 0.30, 1.0, 1.0);
    pad1->SetBottomMargin(0.02);
    pad1->SetLeftMargin(0.14);
    pad1->SetTopMargin(0.06);
    pad1->SetRightMargin(0.05);
    pad1->Draw();

    TPad* pad2 = new TPad("pad_ratio", "pad_ratio", 0.0, 0.0, 1.0, 0.30);
    pad2->SetTopMargin(0.03);
    pad2->SetBottomMargin(0.35);
    pad2->SetLeftMargin(0.14);
    pad2->SetRightMargin(0.05);
    pad2->Draw();

    // ---------------- Main panel ----------------
    pad1->cd();

    h_data->SetMarkerStyle(20);
    h_data->SetMarkerSize(0.8);
    h_data->SetMarkerColor(kBlack);
    h_data->SetLineColor(kBlack);
    h_data->SetStats(0);
    h_data->GetYaxis()->SetTitle("Events");
    h_data->GetYaxis()->SetTitleSize(0.055);
    h_data->GetYaxis()->SetTitleOffset(1.15);
    h_data->GetYaxis()->SetLabelSize(0.045);
    h_data->GetYaxis()->SetNdivisions(510);
    h_data->GetXaxis()->SetLabelSize(0);
    double ymax = h_data->GetMaximum();

    h_data->SetMaximum(1.2 * ymax);
    h_data->SetMinimum(0.0);
    h_data->Draw("E");

    // Background (blue dashed)
    f_background->SetLineColor(kBlue);
    f_background->SetLineStyle(2);
    f_background->SetLineWidth(2);
    f_background->SetNpx(500);
    f_background->Draw("same");

    // Signal (orange dashed) - the (optionally xsignalScale) overlay clone.
    f_signal_draw->SetLineColor(kOrange + 1);
    f_signal_draw->SetLineStyle(2);
    f_signal_draw->SetLineWidth(2);
    f_signal_draw->SetNpx(500);
    f_signal_draw->Draw("same");

    // Signal+Background (red solid)
    g_total->SetLineColor(kRed);
    g_total->SetLineStyle(1);
    g_total->SetLineWidth(2);
    g_total->Draw("L same");

    // Re-draw data on top so the points sit above the curves
    h_data->Draw("E same");

    // ATLAS label + lumi
    TLatex tex;
    tex.SetNDC();
    tex.SetTextFont(42);
    tex.SetTextSize(0.055);
    tex.DrawLatex(0.18, 0.86, "#font[72]{ATLAS} Preliminary");
    tex.SetTextSize(0.035);
    tex.DrawLatex(0.18, 0.79, lumiLabel.c_str());
    // Applied-calibration tag (Nominal / bJR10v00 / bJR10v01 / #mu-in-jet), drawn
    // bold just above the topology label.
    if (!variantLabel.empty()) {
        tex.SetTextSize(0.040);
        tex.DrawLatex(0.18, 0.51, Form("#bf{%s}", variantLabel.c_str()));
        tex.SetTextSize(0.035);
    }
    tex.DrawLatex(0.18, 0.45, topologyLabel.c_str());
    tex.DrawLatex(0.18, 0.4, "#it{R} = 1.0 UFO Soft Drop");
    tex.DrawLatex(0.18, 0.35, "GN2X (0.46% QCD WP)");
    if (!ptLabel.empty()) tex.DrawLatex(0.18, 0.3, ptLabel.c_str());
    // Extracted resolution (width) from the data signal Gaussian fit.
    const double sigmaErr = signal_gauss->GetParError(2);
    tex.DrawLatex(0.18, 0.25, Form("#sigma_{Data} = %.1f #pm %.1f GeV", sigma, sigmaErr));

    // Mark the background-function name (and flag the dual-validation winner).
    std::string bgName = f_background->GetName();
    tex.SetTextSize(0.036);

    // Legend
    TLegend* leg = new TLegend(0.58, 0.66, 0.93, 0.90);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.040);

    // Signal model label.
    int sigNp = signal_gauss->GetNpar();
    std::string sigModel = (sigNp == 7) ? " (DSCB)" : (sigNp == 5) ? " (CB)" : (sigNp == 4) ? " (Bukin)" : "";
    std::string sigLabel = (signalScale != 1.0) ? Form("Signal%s #times %g", sigModel.c_str(), signalScale)
                                                : std::string("Signal") + sigModel;
    leg->AddEntry(h_data, "Data", "lep");
    leg->AddEntry(g_total, "Signal+Backgrounds", "l");
    leg->AddEntry(f_signal_draw, sigLabel.c_str(), "l");
    leg->AddEntry(f_background, "Backgrounds", "l");
    leg->Draw();

    // ---------------- Data / MC ratio panel ----------------
    pad2->cd();

    // Per-bin ratio of data to the total fitted model (background + signal).
    // Error bars carry the data statistical uncertainty, err_data / model.
    TH1F* h_ratio = (TH1F*)h_data->Clone("h_ratio");
    h_ratio->SetDirectory(0);
    h_ratio->Reset();
    for (int i = 1; i <= h_data->GetNbinsX(); ++i) {
        double x = h_data->GetBinCenter(i);
        double obs = h_data->GetBinContent(i);
        double err = h_data->GetBinError(i);
        double exp = std::max(0.0, f_background->Eval(x)) + f_signal->Eval(x);
        double ratio = (exp > 0) ? obs / exp : 0.0;
        double ratioErr = (exp > 0) ? err / exp : 0.0;
        h_ratio->SetBinContent(i, ratio);
        h_ratio->SetBinError(i, ratioErr);
    }

    h_ratio->SetStats(0);
    h_ratio->SetMarkerStyle(20);
    h_ratio->SetMarkerSize(0.8);
    h_ratio->SetMarkerColor(kBlack);
    h_ratio->SetLineColor(kBlack);
    h_ratio->SetMaximum(ratioMax);
    h_ratio->SetMinimum(ratioMin);
    h_ratio->GetYaxis()->SetTitle("Data / (S+B)");
    h_ratio->GetYaxis()->SetNdivisions(505);
    h_ratio->GetYaxis()->SetTitleSize(0.13);
    h_ratio->GetYaxis()->SetTitleOffset(0.42);
    h_ratio->GetYaxis()->SetLabelSize(0.11);
    h_ratio->GetXaxis()->SetTitle("Large-#it{R} (#it{b}#bar{#it{b}})-jet mass [GeV]");
    h_ratio->GetXaxis()->SetTitleSize(0.14);
    h_ratio->GetXaxis()->SetTitleOffset(1.05);
    h_ratio->GetXaxis()->SetLabelSize(0.11);
    h_ratio->Draw("E");

    // Reference line at unity
    TLine* l_unity = new TLine(xmin, 1.0, xmax, 1.0);
    l_unity->SetLineStyle(2);
    l_unity->SetLineColor(kBlack);
    l_unity->SetLineWidth(1);
    l_unity->Draw("same");

    std::string outPath = outDir;
    if (!outPath.empty() && outPath.back() != '/') outPath += "/";
    c->SaveAs((outPath + "signal_background_fit.pdf").c_str());
    c->SaveAs((outPath + "signal_background_fit.png").c_str());
    // For the dual-validation winner, also write an easily-identifiable copy
    // so the optimal plot can be located without inspecting every family folder.
    if (isOptimal) {
        c->SaveAs((outPath + "signal_background_fit_OPTIMAL.pdf").c_str());
        c->SaveAs((outPath + "signal_background_fit_OPTIMAL.png").c_str());
    }

    delete g_total;
    delete f_signal;
    delete h_ratio;
    delete leg;
    delete c;
}

// Cross-check of DrawSignalBackgroundPlot with the standard ATLAS method:
// the signal shape is the MC Z->bb template, and its yield N_sig is fitted
// together with the background parameters in one chi2 fit to the data over
// the full mass range,
//
//   model(m) = N_sig * s(m) * binW + B(m)
//
// with s(m) the MC template normalised to unit area and B the same background
// function as the data-driven plot, starting from its sideband fit values.
// If the sideband subtraction is unbiased, the two plots should agree.
//
// Arguments are the same as DrawSignalBackgroundPlot, except
//   f_background_seed  background fit used as the starting point
//   h_mc_signal        MC Z->bb mass template (signal shape)
//
// Output: signal_background_fit_MCtemplate.{pdf,png}, plus *_OPTIMAL if set.
inline void DrawSignalBackgroundPlot_MCTemplate(TH1F* h_data, TF1* f_background_seed, TH1F* h_mc_signal,
                                                const std::string& outDir, const std::string& topologyLabel,
                                                const std::string& lumiLabel = "", bool isOptimal = false,
                                                const std::string& ptLabel = "", double signalScale = 1.0,
                                                double ratioMin = 0.65, double ratioMax = 1.35,
                                                const std::string& variantLabel = "") {
    if (!h_data || !f_background_seed || !h_mc_signal) {
        std::cerr << "[Mass_SB_Plot] WARN: null input, skipping MC-template S+B plot.\n";
        return;
    }
    const double xmin = h_data->GetXaxis()->GetXmin();
    const double xmax = h_data->GetXaxis()->GetXmax();
    const double binW = h_data->GetBinWidth(1);

    // --- Normalised MC signal template (unit-area density) -------------------
    TH1F* hsig = (TH1F*)h_mc_signal->Clone("h_mc_sig_norm_tmpl");
    hsig->SetDirectory(0);
    double sigArea = hsig->Integral("width");
    if (sigArea <= 0) {
        std::cerr << "[Mass_SB_Plot] WARN: empty MC signal template.\n";
        delete hsig;
        return;
    }
    hsig->Scale(1.0 / sigArea); // now Integral("width") == 1

    // --- Background clone (form + seed params from the sideband fit) ---------
    TF1* fbkg = (TF1*)f_background_seed->Clone("f_bkg_mctmpl");
    const int nbkg = fbkg->GetNpar();

    // --- Seeds: split total data yield into bkg (from seed) + signal ---------
    double Ndata = h_data->Integral();
    double Nbkg_seed = 0.0;
    for (int i = 1; i <= h_data->GetNbinsX(); ++i)
        Nbkg_seed += std::max(0.0, fbkg->Eval(h_data->GetBinCenter(i)));
    double Nsig_seed = std::max(1.0, Ndata - Nbkg_seed);

    // --- Combined S+B model: p[0]=N_sig (events), p[1..]=background params ---
    auto modelFn = [hsig, fbkg, binW, nbkg](double* x, double* p) -> double {
        for (int i = 0; i < nbkg; ++i) fbkg->SetParameter(i, p[1 + i]);
        double dens = hsig->Interpolate(x[0]);
        if (dens < 0) dens = 0.0;
        double bkg = fbkg->Eval(x[0]);
        if (bkg < 0) bkg = 0.0;
        return p[0] * dens * binW + bkg;
    };
    TF1* fSB = new TF1("f_sb_mctmpl", modelFn, xmin, xmax, 1 + nbkg);
    fSB->SetParName(0, "N_sig");
    fSB->SetParameter(0, Nsig_seed);
    fSB->SetParLimits(0, 0.0, 1.5 * Ndata);
    for (int i = 0; i < nbkg; ++i) fSB->SetParameter(1 + i, f_background_seed->GetParameter(i));

    h_data->Fit(fSB, "RQN0", "", xmin, xmax);


    const double Nsig_fit = fSB->GetParameter(0);
    const double chi2 = fSB->GetChisquare();
    const int ndf = fSB->GetNDF();

    for (int i = 0; i < nbkg; ++i) fbkg->SetParameter(i, fSB->GetParameter(1 + i));

    // --- Signal overlay curve (optionally xsignalScale for visibility) ------
    const double sDraw = Nsig_fit * signalScale;
    auto sigFn = [hsig, binW, sDraw](double* x, double*) -> double {
        double d = hsig->Interpolate(x[0]);
        if (d < 0) d = 0.0;
        return sDraw * d * binW;
    };
    TF1* f_sig_draw = new TF1("f_sig_mctmpl_draw", sigFn, xmin, xmax, 0);

    // --- Canvas (same layout as the data-driven plot) -----------------------
    TCanvas* c = new TCanvas("c_sb_mctmpl", "MC-template S+B fit", 700, 750);
    TPad* pad1 = new TPad("pad_main_t", "", 0.0, 0.30, 1.0, 1.0);
    pad1->SetBottomMargin(0.02);
    pad1->SetLeftMargin(0.14);
    pad1->SetTopMargin(0.06);
    pad1->SetRightMargin(0.05);
    pad1->Draw();
    TPad* pad2 = new TPad("pad_ratio_t", "", 0.0, 0.0, 1.0, 0.30);
    pad2->SetTopMargin(0.03);
    pad2->SetBottomMargin(0.35);
    pad2->SetLeftMargin(0.14);
    pad2->SetRightMargin(0.05);
    pad2->Draw();

    pad1->cd();
    h_data->SetMarkerStyle(20);
    h_data->SetMarkerSize(0.8);
    h_data->SetMarkerColor(kBlack);
    h_data->SetLineColor(kBlack);
    h_data->SetStats(0);
    h_data->GetYaxis()->SetTitle("Events");
    h_data->GetYaxis()->SetTitleSize(0.055);
    h_data->GetYaxis()->SetTitleOffset(1.15);
    h_data->GetYaxis()->SetLabelSize(0.045);
    h_data->GetYaxis()->SetNdivisions(510);
    h_data->GetXaxis()->SetLabelSize(0);
    double ymax = h_data->GetMaximum();
    h_data->SetMaximum(1.05 * ymax);
    h_data->SetMinimum(0.0);
    h_data->Draw("E");

    fbkg->SetLineColor(kBlue);
    fbkg->SetLineStyle(2);
    fbkg->SetLineWidth(2);
    fbkg->SetNpx(500);
    fbkg->Draw("same");
    f_sig_draw->SetLineColor(kOrange + 1);
    f_sig_draw->SetLineStyle(2);
    f_sig_draw->SetLineWidth(2);
    f_sig_draw->SetNpx(500);
    f_sig_draw->Draw("same");
    fSB->SetLineColor(kRed);
    fSB->SetLineStyle(1);
    fSB->SetLineWidth(2);
    fSB->SetNpx(500);
    fSB->Draw("same");
    h_data->Draw("E same");

    TLatex tex;
    tex.SetNDC();
    tex.SetTextFont(42);
    tex.SetTextSize(0.055);
    tex.DrawLatex(0.18, 0.86, "#font[72]{ATLAS} Preliminary");
    tex.SetTextSize(0.042);
    tex.DrawLatex(0.18, 0.79, lumiLabel.c_str());

    if (!variantLabel.empty()) {
        tex.SetTextSize(0.046);
        tex.DrawLatex(0.18, 0.50, Form("#bf{%s}", variantLabel.c_str()));
        tex.SetTextSize(0.042);
    }
    tex.DrawLatex(0.18, 0.44, topologyLabel.c_str());
    tex.DrawLatex(0.18, 0.38, "#it{R} = 1.0 UFO Soft Drop");
    tex.DrawLatex(0.18, 0.32, "GN2X (0.46% QCD WP)");
    if (!ptLabel.empty()) tex.DrawLatex(0.18, 0.26, ptLabel.c_str());
    tex.SetTextSize(0.034);
    tex.DrawLatex(0.18, 0.20, "MC-template S+B fit to data");
    tex.DrawLatex(0.18, 0.15, Form("N_{sig} = %.0f,  #chi^{2}/ndf = %.1f/%d", Nsig_fit, chi2, ndf));

    TLegend* leg = new TLegend(0.58, 0.66, 0.93, 0.90);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.040);
    std::string sigLabel =
        (signalScale != 1.0) ? Form("Signal (MC) #times %g", signalScale) : std::string("Signal (MC)");
    leg->AddEntry(h_data, "Data", "lep");
    leg->AddEntry(fSB, "Signal+Backgrounds", "l");
    leg->AddEntry(f_sig_draw, sigLabel.c_str(), "l");
    leg->AddEntry(fbkg, "Backgrounds", "l");
    leg->Draw();

    pad2->cd();
    TH1F* h_ratio = (TH1F*)h_data->Clone("h_ratio_mctmpl");
    h_ratio->SetDirectory(0);
    h_ratio->Reset();
    for (int i = 1; i <= h_data->GetNbinsX(); ++i) {
        double x = h_data->GetBinCenter(i), obs = h_data->GetBinContent(i), err = h_data->GetBinError(i);
        double exp = std::max(0.0, fSB->Eval(x));
        h_ratio->SetBinContent(i, exp > 0 ? obs / exp : 0.0);
        h_ratio->SetBinError(i, exp > 0 ? err / exp : 0.0);
    }
    h_ratio->SetStats(0);
    h_ratio->SetMarkerStyle(20);
    h_ratio->SetMarkerSize(0.8);
    h_ratio->SetMarkerColor(kBlack);
    h_ratio->SetLineColor(kBlack);
    h_ratio->SetMaximum(ratioMax);
    h_ratio->SetMinimum(ratioMin);
    h_ratio->GetYaxis()->SetTitle("Data / (S+B)");
    h_ratio->GetYaxis()->SetNdivisions(505);
    h_ratio->GetYaxis()->SetTitleSize(0.13);
    h_ratio->GetYaxis()->SetTitleOffset(0.42);
    h_ratio->GetYaxis()->SetLabelSize(0.11);
    h_ratio->GetXaxis()->SetTitle("Large-#it{R} (#it{b}#bar{#it{b}})-jet mass [GeV]");
    h_ratio->GetXaxis()->SetTitleSize(0.14);
    h_ratio->GetXaxis()->SetTitleOffset(1.05);
    h_ratio->GetXaxis()->SetLabelSize(0.11);
    h_ratio->Draw("E");
    TLine* l_unity = new TLine(xmin, 1.0, xmax, 1.0);
    l_unity->SetLineStyle(2);
    l_unity->SetLineColor(kBlack);
    l_unity->SetLineWidth(1);
    l_unity->Draw("same");

    std::string outPath = outDir;
    if (!outPath.empty() && outPath.back() != '/') outPath += "/";
    c->SaveAs((outPath + "signal_background_fit_MCtemplate.pdf").c_str());
    c->SaveAs((outPath + "signal_background_fit_MCtemplate.png").c_str());
    if (isOptimal) {
        c->SaveAs((outPath + "signal_background_fit_MCtemplate_OPTIMAL.pdf").c_str());
        c->SaveAs((outPath + "signal_background_fit_MCtemplate_OPTIMAL.png").c_str());
    }

    delete fSB;
    delete f_sig_draw;
    delete fbkg;
    delete hsig;
    delete h_ratio;
    delete l_unity;
    delete leg;
    delete c;
}

#endif // MASS_SB_PLOT_H
