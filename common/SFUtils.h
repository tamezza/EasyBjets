/**
 * =============================================================================================
 *  SFUtils.h — mass-dependent scale factor for the main background
 * =============================================================================================
 * ===========================================
 * @author:    Thiziri AMEZZA
 * @email:     thiziri.amezza@cern.ch
 * @date:      18 April 2026 
 * ===========================================
 *
 * The main background (gamma+jets for DB, multijets for MJB) is normalised
 * to data in the mass sidebands on either side of the Z peak. A linear
 * dependence,
 *     SF(m) = p0 + p1 * m,
 * is fitted to the Data/MC ratio of the sideband bins and then applied
 *   - bin by bin to mass histograms (ApplySFBinByBin), which is the same
 *     thing as weighting events when the histogram axis is the mass itself,
 *   - event by event to anything else, e.g. the response R (GetEventWeight).
 * ==========================================================================
 * @version:  1.0
 * ==========================================================================
 */
#ifndef SFUTILS_H
#define SFUTILS_H

#include <vector>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>
#include "TString.h"

#include "TH1F.h"
#include "TF1.h"
#include "TFitResult.h"   
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TLine.h"
#include "TBox.h"
#include "TLegend.h"
#include "TTree.h"

#include "CommonUtils.h" // Log::, Validate::, AtlasUtils::

namespace SFUtils {

// ============================================================================
// LinearSF — Structure to hold mass-dependent scale factor parameters
// ============================================================================
// SF(m) = p0 + p1 × m
//
// The scale factor corrects the γ+jets MC mismodelling to match data in the
// mass sidebands (outside the Z peak). A linear dependence on mass accounts
// for shape differences between γ+jets MC and data.
// ============================================================================
struct LinearSF {
    double p0 = 1.0;        // Intercept
    double p1 = 0.0;        // Slope (per GeV)
    double p0_err = 0.0;    
    double p1_err = 0.0;    
    double cov_p0_p1 = 0.0; // Covariance between intercept and slope - needed for Error(m)
    double chi2 = 0.0;      
    int ndf = 0;            
    int nPoints = 0;        
    bool valid = false;     
    
    // Evaluate SF at a given mass
    double operator()(double m) const { 
        return p0 + p1 * m; 
    }
    
    // sigma^2(SF) = sigma^2(p0) + m^2 sigma^2(p1) + 2 m cov(p0,p1)
    double Error(double m) const {
        double var = p0_err * p0_err
                   + m * m * p1_err * p1_err
                   + 2.0 * m * cov_p0_p1;
        return (var > 0.0) ? std::sqrt(var) : 0.0;
    }
    
    // Print summary
    void Print(const std::string& ctx = "SFUtils") const {
        if (valid) {
            Log::info(ctx, Form("SF(m) = %.4f (± %.4f) + %.8f (± %.8f) × m",
                      p0, p0_err, p1, p1_err));
            double denom = p0_err * p1_err;
            double corr = (denom > 0.0) ? cov_p0_p1 / denom : 0.0;
            Log::info(ctx, Form("  cov(p0,p1) = %.4g,  corr = %.3f", cov_p0_p1, corr));
            Log::info(ctx, Form("  χ²/ndf = %.2f/%d = %.2f,  nPoints = %d",
                      chi2, ndf, (ndf > 0 ? chi2/ndf : 0.0), nPoints));
        } else {
            Log::warning(ctx, "SF fit failed — using SF = 1.0");
        }
    }
};


// Fit the Data/MC ratio in [sb1_lo, sb1_hi] and [sb2_lo, sb2_hi].
// Bins with no MC are skipped. Needs at least three points.
inline LinearSF FitMassDependentSF(
    TH1F* hData, 
    TH1F* hBkg,
    double sb1_lo, 
    double sb1_hi,
    double sb2_lo, 
    double sb2_hi,
    const std::string& diagCanvasPath = "",
    const std::string& extraLabel = "",
    const std::string& bkgLegendLabel = "#gamma+jets")  
{
    const std::string ctx = "SFUtils";
    LinearSF sf;
    
    if (!hData || !hBkg) {
        Log::warning(ctx, "Null histogram(s) — SF = 1.0");
        return sf;
    }
    
    if (hData->GetEntries() == 0 || hBkg->GetEntries() == 0) {
        Log::warning(ctx, "Empty histogram(s) — SF = 1.0");
        return sf;
    }

    std::vector<double> x_vals, y_vals, y_errs;
    
    auto addSidebandPoints = [&](double lo, double hi) {
        int binLo = hBkg->FindBin(lo);
        int binHi = hBkg->FindBin(hi);
        // Avoid double-counting bin edges
        if (hBkg->GetBinCenter(binHi) > hi) binHi--;
        
        for (int b = binLo; b <= binHi; ++b) {
            double xc = hBkg->GetBinCenter(b);
            double d  = hData->GetBinContent(b);
            double mc = hBkg->GetBinContent(b);
            
            if (mc <= 0) continue;  
            
            double ratio = d / mc;
            
            // Error propagation: 
            double d_err  = hData->GetBinError(b);
            double mc_err = hBkg->GetBinError(b);
            double ratio_err = std::sqrt(
                std::pow(d_err / mc, 2) +
                std::pow(d * mc_err / (mc * mc), 2)
            );
            
            x_vals.push_back(xc);
            y_vals.push_back(ratio);
            y_errs.push_back(ratio_err);
        }
    };
    
    addSidebandPoints(sb1_lo, sb1_hi);
    addSidebandPoints(sb2_lo, sb2_hi);
    
    sf.nPoints = static_cast<int>(x_vals.size());
    
    if (sf.nPoints < 3) {
        Log::warning(ctx, Form("Only %d sideband points — need ≥3 for fit. SF = 1.0", sf.nPoints));
        return sf;
    }
    
    TGraphErrors* gr = new TGraphErrors(sf.nPoints);
    for (int i = 0; i < sf.nPoints; ++i) {
        gr->SetPoint(i, x_vals[i], y_vals[i]);
        gr->SetPointError(i, 0, y_errs[i]);
    }
    
    double fitRangeLo = std::min(sb1_lo, sb2_lo);
    double fitRangeHi = std::max(sb1_hi, sb2_hi);
    
    TF1* fit = new TF1("sf_linear_fit", "[0] + [1]*x", fitRangeLo, fitRangeHi);
    fit->SetParameters(1.0, 0.0);
    fit->SetParNames("p0", "p1");
    
    // Fit with quiet option
    auto fitResult = gr->Fit(fit, "RQS");
    int fitStatus = fitResult;
    
    if (fitStatus != 0) {
        Log::warning(ctx, Form("Fit failed with status %d — SF = 1.0", fitStatus));
        delete gr;
        delete fit;
        return sf;
    }
    
    sf.p0 = fit->GetParameter(0);
    sf.p1 = fit->GetParameter(1);
    sf.p0_err = fit->GetParError(0);
    sf.p1_err = fit->GetParError(1);
    if (fitResult.Get() && fitResult->NPar() >= 2) {
        sf.cov_p0_p1 = fitResult->CovMatrix(0, 1);
    }
    sf.chi2 = fit->GetChisquare();
    sf.ndf = fit->GetNDF();
    sf.valid = true;
    sf.Print(ctx);
    
    // Optional: save diagnostic canvas
    if (!diagCanvasPath.empty()) {
        TCanvas* c = new TCanvas("sf_fit_canvas", "Sideband SF Fit", 900, 700);
        c->SetLeftMargin(0.12);
        c->SetRightMargin(0.05);
        c->SetBottomMargin(0.12);
        c->SetTopMargin(0.08);
        
        gr->SetTitle(Form(";M^{J} [GeV];Data / %s MC", bkgLegendLabel.c_str()));
        gr->SetMarkerStyle(20);
        gr->SetMarkerSize(1.0);
        gr->SetMarkerColor(kBlue+1);
        gr->SetLineColor(kBlue+1);
        gr->SetLineWidth(2);
        
        // Set axis ranges
        gr->GetXaxis()->SetLimits(fitRangeLo - 10, fitRangeHi + 10);
        double yMin = 0.0, yMax = 2.5;
        for (int i = 0; i < sf.nPoints; ++i) {
            double y = y_vals[i];
            double ye = y_errs[i];
            if (y - ye < yMin) yMin = y - ye - 0.1;
            if (y + ye > yMax) yMax = y + ye + 0.1;
        }
        gr->GetYaxis()->SetRangeUser(std::max(0.0, yMin), yMax);
        
        gr->GetXaxis()->SetTitleSize(0.045);
        gr->GetYaxis()->SetTitleSize(0.045);
        gr->GetXaxis()->SetLabelSize(0.04);
        gr->GetYaxis()->SetLabelSize(0.04);
        gr->GetYaxis()->SetTitleOffset(1.1);
        
        gr->Draw("AP");
        
        // Draw fit line
        fit->SetLineColor(kRed);
        fit->SetLineWidth(2);
        fit->Draw("same");
        
        // Draw sideband regions as shaded boxes
        double boxY1 = c->GetUymin();
        double boxY2 = c->GetUymax();
        
        TBox* box1 = new TBox(sb1_lo, 0, sb1_hi, yMax);
        box1->SetFillColorAlpha(kGreen, 0.15);
        box1->SetLineColor(kGreen+2);
        box1->SetLineStyle(2);
        box1->Draw("same");
        
        TBox* box2 = new TBox(sb2_lo, 0, sb2_hi, yMax);
        box2->SetFillColorAlpha(kGreen, 0.15);
        box2->SetLineColor(kGreen+2);
        box2->SetLineStyle(2);
        box2->Draw("same");
        
        // Draw Z peak exclusion region
        TBox* boxZ = new TBox(sb1_hi, 0, sb2_lo, yMax);
        boxZ->SetFillColorAlpha(kRed, 0.08);
        boxZ->SetLineColor(kRed);
        boxZ->SetLineStyle(2);
        boxZ->Draw("same");
        
        // Redraw graph on top
        gr->Draw("P same");
        
        // Draw unity line
        TLine* unity = new TLine(fitRangeLo - 10, 1.0, fitRangeHi + 10, 1.0);
        unity->SetLineColor(kGray+1);
        unity->SetLineStyle(2);
        unity->Draw("same");
        
        // Labels
        TLatex lx;
        lx.SetNDC();
        lx.SetTextFont(42);
        lx.SetTextSize(0.025);
        
        lx.DrawLatex(0.15, 0.80, Form("SF(m) = %.4f + %.8f #times m", sf.p0, sf.p1));
        lx.DrawLatex(0.15, 0.77, Form("#chi^{2}/ndf = %.2f/%d = %.2f", 
                     sf.chi2, sf.ndf, (sf.ndf > 0 ? sf.chi2/sf.ndf : 0.0)));
        lx.DrawLatex(0.15, 0.74, Form("Sidebands: [%.0f, %.0f] #cup [%.0f, %.0f] GeV",
                     sb1_lo, sb1_hi, sb2_lo, sb2_hi));
        
        if (!extraLabel.empty()) {
            lx.DrawLatex(0.15, 0.71, extraLabel.c_str());
        }
        
        // Legend
        TLegend* leg = new TLegend(0.55, 0.75, 0.92, 0.90);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.025);
        leg->AddEntry(gr, Form("Data / %s (sidebands)", bkgLegendLabel.c_str()), "lep");
        leg->AddEntry(fit, "Linear fit: p_{0} + p_{1} #times m", "l");
        leg->AddEntry(box1, "Sideband regions", "f");
        leg->AddEntry(boxZ, "Z peak exclusion", "f");
        leg->Draw();
        
        AtlasUtils::DrawAtlasLabel(0.15, 0.87);
        
        c->SaveAs(diagCanvasPath.c_str());
        Log::info(ctx, "SF diagnostic plot saved: " + diagCanvasPath);
        
        delete c;
    }
    
    delete gr;
    delete fit;
    
    return sf;
}


// ============================================================================
// FitMassDependentSF_PtBinned — Fit SF per pT bin
// ============================================================================
// For analyses that use pT-dependent SF returns a vector of LinearSF, one per pT bin.
//
// INPUTS:
//   tree      - TTree to read events from
//   nPtBins   - Number of pT bins
//   ptBins    - Array of pT bin edges [nPtBins+1]
//   branchPt  - Branch name for jet pT
//   branchM   - Branch name for jet mass
//   branchW   - Branch name for weight
//   branchIdx - Branch name for selected jet index
//   isData    - If true, weight = 1
//   sb1_lo, sb1_hi, sb2_lo, sb2_hi — Sideband boundaries
//   diagDir   - Directory for diagnostic plots (empty = skip)
// ============================================================================
inline std::vector<LinearSF> FitMassDependentSF_PtBinned(
    TTree* treeData,
    TTree* treeBkg,
    int nPtBins,
    const double* ptBins,
    const std::string& branchPt,
    const std::string& branchM,
    const std::string& branchW,
    const std::string& branchIdx,
    double sb1_lo, double sb1_hi,
    double sb2_lo, double sb2_hi,
    const std::string& diagDir = "",
    const std::string& labelPrefix = "",
    const std::string& bkgLegendLabel = "#gamma+jets")  
{
    const std::string ctx = "SFUtils::PtBinned";
    std::vector<LinearSF> sfVec(nPtBins);
    
    if (!treeData || !treeBkg) {
        Log::warning(ctx, "Null tree(s) — returning SF = 1.0 for all bins");
        return sfVec;
    }
    
    // Create mass histograms per pT bin
    std::vector<TH1F*> h_data(nPtBins);
    std::vector<TH1F*> h_bkg(nPtBins);
    
    for (int b = 0; b < nPtBins; ++b) {
        h_data[b] = new TH1F(Form("h_data_sf_ptbin%d", b), "", 20, 50, 150);
        h_bkg[b]  = new TH1F(Form("h_bkg_sf_ptbin%d", b),  "", 20, 50, 150);
        h_data[b]->Sumw2();
        h_bkg[b]->Sumw2();
    }
    
    // Branch variables
    std::vector<float>* ljet_pt = nullptr;
    std::vector<float>* ljet_m = nullptr;
    double weight = 1.0;
    int selIdx = -1;
    
    // Lambda to find pT bin
    auto findPtBin = [&](double pt) -> int {
        for (int b = 0; b < nPtBins; ++b) {
            if (pt >= ptBins[b] && pt < ptBins[b+1]) return b;
        }
        return -1;
    };
    
    // ---- Validate and set branches for DATA tree ----
    bool dataHasPt  = Validate::HasBranch(treeData, branchPt);
    bool dataHasM   = Validate::HasBranch(treeData, branchM);
    bool dataHasIdx = Validate::HasBranch(treeData, branchIdx);
    
    if (!dataHasIdx) {
        Log::warning(ctx, "Data tree missing '" + branchIdx + "' — cannot compute SF");
        for (int b = 0; b < nPtBins; ++b) { delete h_data[b]; delete h_bkg[b]; }
        return sfVec;
    }
    
    std::string dataPtBranch  = dataHasPt  ? branchPt : "ljet_pt";
    std::string dataMBranch   = dataHasM   ? branchM  : "ljet_m";
    
    if (!dataHasPt) Log::info(ctx, "Data: using fallback branch 'ljet_pt'");
    if (!dataHasM)  Log::info(ctx, "Data: using fallback branch 'ljet_m'");
    
    treeData->SetBranchAddress(dataPtBranch.c_str(), &ljet_pt);
    treeData->SetBranchAddress(dataMBranch.c_str(), &ljet_m);
    treeData->SetBranchAddress(branchIdx.c_str(), &selIdx);
    
    Long64_t nData = treeData->GetEntries();
    for (Long64_t i = 0; i < nData; ++i) {
        treeData->GetEntry(i);
        if (selIdx < 0 || !ljet_pt || !ljet_m ||
            selIdx >= (int)ljet_pt->size() || selIdx >= (int)ljet_m->size()) continue;
        
        double pt = (*ljet_pt)[selIdx];
        double m  = (*ljet_m)[selIdx];
        int bin = findPtBin(pt);
        if (bin >= 0) h_data[bin]->Fill(m);
    }
    
    // ---- Validate and set branches for BACKGROUND tree ----
    bool bkgHasPt  = Validate::HasBranch(treeBkg, branchPt);
    bool bkgHasM   = Validate::HasBranch(treeBkg, branchM);
    bool bkgHasW   = Validate::HasBranch(treeBkg, branchW);
    bool bkgHasIdx = Validate::HasBranch(treeBkg, branchIdx);
    
    if (!bkgHasPt || !bkgHasM || !bkgHasIdx) {
        Log::warning(ctx, "Bkg tree missing required branches — cannot compute SF");
        for (int b = 0; b < nPtBins; ++b) { delete h_data[b]; delete h_bkg[b]; }
        return sfVec;
    }
    
    treeBkg->SetBranchAddress(branchPt.c_str(), &ljet_pt);
    treeBkg->SetBranchAddress(branchM.c_str(), &ljet_m);
    treeBkg->SetBranchAddress(branchIdx.c_str(), &selIdx);
    if (bkgHasW) {
        treeBkg->SetBranchAddress(branchW.c_str(), &weight);
    } else {
        Log::warning(ctx, "Bkg tree missing '" + branchW + "' — using weight = 1.0");
        weight = 1.0;
    }
    
    Long64_t nBkg = treeBkg->GetEntries();
    for (Long64_t i = 0; i < nBkg; ++i) {
        treeBkg->GetEntry(i);
        if (selIdx < 0 || !ljet_pt || !ljet_m ||
            selIdx >= (int)ljet_pt->size() || selIdx >= (int)ljet_m->size()) continue;
        
        double pt = (*ljet_pt)[selIdx];
        double m  = (*ljet_m)[selIdx];
        int bin = findPtBin(pt);
        if (bin >= 0) h_bkg[bin]->Fill(m, weight);
    }
    
    // Fit SF for each pT bin
    for (int b = 0; b < nPtBins; ++b) {
        std::string diagPath = "";
        if (!diagDir.empty()) {
            diagPath = diagDir + Form("/sf_fit_ptbin%d.pdf", b);
        }
        std::string label = labelPrefix;
        if (!label.empty()) label += ", ";
        label += Form("p_{T} #in [%.0f, %.0f] GeV", ptBins[b], ptBins[b+1]);
        
        sfVec[b] = FitMassDependentSF(h_data[b], h_bkg[b],
                                       sb1_lo, sb1_hi, sb2_lo, sb2_hi,
                                       diagPath, label, bkgLegendLabel);
        
        Log::info(ctx, Form("pT bin %d [%.0f, %.0f]: SF = %.4f + %.8f × m",
                  b, ptBins[b], ptBins[b+1], sfVec[b].p0, sfVec[b].p1));
    }
    
    // Cleanup
    for (int b = 0; b < nPtBins; ++b) {
        delete h_data[b];
        delete h_bkg[b];
    }
    
    return sfVec;
}


// Scale each bin by SF(bin centre).
inline void ApplySFBinByBin(TH1F* h, const LinearSF& sf) {
    if (!h) return;
    if (!sf.valid) {
        Log::warning("SFUtils", "SF not valid — histogram unchanged");
        return;
    }
    
    for (int bin = 1; bin <= h->GetNbinsX(); ++bin) {
        double m = h->GetBinCenter(bin);
        double s = sf(m);
        if (s < 0) s = 0;
        
        h->SetBinContent(bin, h->GetBinContent(bin) * s);
        h->SetBinError(bin, h->GetBinError(bin) * s);
    }
    
    Log::info("SFUtils", Form("SF applied bin-by-bin to '%s'", h->GetName()));
}


// SF uncertainty added in quadrature to the bin error.
inline void ApplySFBinByBinWithError(TH1F* h, const LinearSF& sf) {
    if (!h) return;
    if (!sf.valid) {
        Log::warning("SFUtils", "SF not valid — histogram unchanged");
        return;
    }
    
    for (int bin = 1; bin <= h->GetNbinsX(); ++bin) {
        double m = h->GetBinCenter(bin);
        double s = sf(m);
        double s_err = sf.Error(m);
        
        double content = h->GetBinContent(bin);
        double stat_err = h->GetBinError(bin);

        double new_content = content * s; 
        double new_err = std::sqrt(
            std::pow(stat_err * s, 2) +
            std::pow(content * s_err, 2)
        );
        
        h->SetBinContent(bin, new_content);
        h->SetBinError(bin, new_err);
    }
    
    Log::info("SFUtils", Form("SF applied bin-by-bin (with SF error) to '%s'", h->GetName()));
}


// ============================================================================
// GetEventWeight — Compute per-event corrected weight
// ============================================================================
// For R_DB histograms where x-axis != mass, must apply SF per-event.
//
// USAGE:
//   double w_corr = SFUtils::GetEventWeight(w_raw, m_event, sf);
//   h_R_DB->Fill(R_DB, w_corr);
// ============================================================================
inline double GetEventWeight(double rawWeight, double mass, const LinearSF& sf) {
    if (!sf.valid) return rawWeight;
    double s = sf(mass);
    if (s < 0) s = 0;
    return rawWeight * s;
}


// ============================================================================
// GetEventWeight_PtBinned — Per-event weight with pT-dependent SF
// ============================================================================
inline double GetEventWeight_PtBinned(
    double rawWeight, 
    double mass, 
    double pt,
    const std::vector<LinearSF>& sfVec,
    int nPtBins,
    const double* ptBins)
{
    // Find pT bin
    int bin = -1;
    for (int b = 0; b < nPtBins; ++b) {
        if (pt >= ptBins[b] && pt < ptBins[b+1]) {
            bin = b;
            break;
        }
    }
    
    if (bin < 0 || bin >= (int)sfVec.size() || !sfVec[bin].valid) {
        return rawWeight;
    }
    
    return rawWeight * sfVec[bin](mass);
}


// ============================================================================
// PrintSFSummary — Print summary table of SF values
// ============================================================================
inline void PrintSFSummary(
    const std::vector<LinearSF>& sfVec,
    int nPtBins,
    const double* ptBins,
    const std::string& title = "Scale Factor Summary")
{
    std::cout << "\n";
    std::cout << "|======================================================================|\n";
    std::cout << "|  " << std::left << std::setw(68) << title << "|\n";
    std::cout << "|======================================================================|\n";
    std::cout << "| " << std::left
              << std::setw(18) << "pT Range [GeV]"
              << std::setw(14) << "p0"
              << std::setw(14) << "p1"
              << std::setw(12) << "chi2/ndf"
              << std::setw(10) << "nPoints"
              << "|\n";
    std::cout << "|======================================================================|\n";
    
    for (int b = 0; b < nPtBins && b < (int)sfVec.size(); ++b) {
        const auto& sf = sfVec[b];
        std::cout << "| " << std::left 
                  << std::setw(18) << Form("[%.0f, %.0f)", ptBins[b], ptBins[b+1])
                  << std::fixed << std::setprecision(4)
                  << std::setw(14) << sf.p0
                  << std::setw(14) << sf.p1
                  << std::setw(12) << (sf.ndf > 0 ? sf.chi2/sf.ndf : 0.0)
                  << std::setw(10) << sf.nPoints
                  << "|\n";
    }
    
    std::cout << "|======================================================================|\n";
}

// Save LinearSF to a text-file handoff between macros (e.g. Mass_Distributions writes it,
// Response_Fit_Inclusive reads it instead of refitting).
inline void SaveLinearSF(const LinearSF& sf,
                          const std::string& path,
                          const std::string& label = "",
                          double sb1_low  = 0.0, double sb1_high = 0.0,
                          double sb2_low  = 0.0, double sb2_high = 0.0)
{
    std::ofstream out(path);
    if (!out.is_open()) {
        Log::error("SFUtils", "Cannot write SF to: " + path);
        return;
    }
    out << "# LinearSF — " << label << "\n"
        << "format_version 2\n"
        << "p0       " << std::setprecision(10) << sf.p0     << "\n"
        << "p1       " << std::setprecision(10) << sf.p1     << "\n"
        << "p0_err   " << std::setprecision(10) << sf.p0_err << "\n"
        << "p1_err   " << std::setprecision(10) << sf.p1_err << "\n"
        << "chi2     " << sf.chi2   << "\n"
        << "ndf      " << sf.ndf    << "\n"
        << "nPoints  " << sf.nPoints<< "\n"
        << "valid    " << sf.valid  << "\n"
        << "sb1_low  " << sb1_low   << "\n"
        << "sb1_high " << sb1_high  << "\n"
        << "sb2_low  " << sb2_low   << "\n"
        << "sb2_high " << sb2_high  << "\n";
    out.close();
    Log::info("SFUtils", "SF saved to: " + path);
}

// Load LinearSF from a text-file written by SaveLinearSF.
//
// When expected sidebands are given, a file fitted with different windows
// (or an old file without the stamp) comes back with valid = false so the
// caller refits instead of using a stale SF.
inline LinearSF LoadLinearSF(const std::string& path,
                              double expected_sb1_low  = 0.0,
                              double expected_sb1_high = 0.0,
                              double expected_sb2_low  = 0.0,
                              double expected_sb2_high = 0.0)
{
    LinearSF sf;
    std::ifstream in(path);
    if (!in.is_open()) {
        Log::warning("SFUtils", "Cannot read SF from: " + path + " — using SF=1");
        return sf;
    }
    int    format_version = 1;
    double sb1_low = 0, sb1_high = 0, sb2_low = 0, sb2_high = 0;
    bool   has_sb_stamp = false;

    std::string key;
    while (in >> key) {
        if      (key == "p0")             in >> sf.p0;
        else if (key == "p1")             in >> sf.p1;
        else if (key == "p0_err")         in >> sf.p0_err;
        else if (key == "p1_err")         in >> sf.p1_err;
        else if (key == "chi2")           in >> sf.chi2;
        else if (key == "ndf")            in >> sf.ndf;
        else if (key == "nPoints")        in >> sf.nPoints;
        else if (key == "valid")          { int v; in >> v; sf.valid = (v != 0); }
        else if (key == "format_version") in >> format_version;
        else if (key == "sb1_low")        { in >> sb1_low;  has_sb_stamp = true; }
        else if (key == "sb1_high")       { in >> sb1_high; has_sb_stamp = true; }
        else if (key == "sb2_low")        { in >> sb2_low;  has_sb_stamp = true; }
        else if (key == "sb2_high")       { in >> sb2_high; has_sb_stamp = true; }
        else { std::string dummy; std::getline(in, dummy); } 
    }
    in.close();

    const bool checkRequested =
        (expected_sb1_low != 0.0 || expected_sb1_high != 0.0 ||
         expected_sb2_low != 0.0 || expected_sb2_high != 0.0);

    if (checkRequested) {
        if (!has_sb_stamp) {
            Log::warning("SFUtils", "Cache " + path +
                " has no sideband stamp (format v1) — invalidating to force refit.");
            sf.valid = false;
        } else {
            const double tol = 1e-6;
            const bool mismatch =
                std::abs(sb1_low  - expected_sb1_low ) > tol ||
                std::abs(sb1_high - expected_sb1_high) > tol ||
                std::abs(sb2_low  - expected_sb2_low ) > tol ||
                std::abs(sb2_high - expected_sb2_high) > tol;
            if (mismatch) {
                Log::warning("SFUtils", Form(
                    "Cache %s sideband mismatch — invalidating to force refit.\n"
                    "  cached:   sb1=[%.2f,%.2f]  sb2=[%.2f,%.2f]\n"
                    "  expected: sb1=[%.2f,%.2f]  sb2=[%.2f,%.2f]",
                    path.c_str(),
                    sb1_low, sb1_high, sb2_low, sb2_high,
                    expected_sb1_low, expected_sb1_high,
                    expected_sb2_low, expected_sb2_high));
                sf.valid = false;
            }
        }
    }

    sf.Print("LoadLinearSF");
    return sf;
}

}  // namespace SFUtils

#endif  // SFUTILS_H
