/**
 * ============================================================================
 *  Mass_Fit.h - driver of the mass peak fit
 * ============================================================================
 * ===========================================
 * @author:    Thiziri AMEZZA
 * @email:     thiziri.amezza@cern.ch
 * @date:      18 April 2026
 * ===========================================
 *
 * For each variant, working point and pT bin:
 *   1. fill the data, signal MC and background MC mass histograms with the
 *      binning from BinningConfig.h (EASYBJETS_BINNING);
 *   2. normalise the background MC to data in the sidebands with a linear SF(m)
 *      (same method as SFUtils, done inline here);
 *   3. pick the background function on the data sidebands with the F-test and
 *      k-fold CV (MassFit_Selection.h), which also produces the per-function plots;
 *   4. subtract it from the data and fit the Z peak of what remains, and of the
 *      signal MC (MassFit_Signal.h): Gaussian by default, Bukin / CB / DSCB with
 *      EASYBJETS_FIT_MODEL;
 *   5. write fit_summary.txt with the peak positions, widths and their ratios.
 * The other mass-fit headers are listed in MassFit_Types.h.
 *
 * How to run (from ZbbGamma/macros for DB, from ZbbJets/macros for MJB):
 *   export EASYBJETS_RUN=run3
 *   export EASYBJETS_DATA_DIR=$HOME/lpnatlas_data/ZbbGamma    (ZbbJets for MJB)
 *   root 'Mass_Fit.cxx("ALL","inclusive")'   every variant, one pT bin
 *   root 'Mass_Fit.cxx("ALL","ptbins")'      every variant, pT bins of the channel
 *   root 'Mass_Fit.cxx("DB_bJER10v00")'      one variant (MJB_... in ZbbJets), inclusive
 * Then run Mass_Ratios_PeakWidth.cxx with the same options.
 *
 * Output: <PlotDir>/<variant>/Mass_Fit_<channel>/<wp>/<ptbin>_<binning>/
 *   fit_summary.txt, ftest_results.txt, cross_validation/, one folder per
 *   background function with its plots
 *
 * Options (environment variables, read once per ROOT session):
 *   EASYBJETS_RUN             run2 | run3 (default)
 *   EASYBJETS_DATA_DIR        folder holding NtupleSlim_Run2/3 of this channel
 *                             (default ../../DataSets/<ZbbGamma|ZbbJets>)
 *   EASYBJETS_BINNING         mass histogram binning: fixed (default) | fd | scott |
 *                             sturges | peakaware; tags the folder (<ptbin>_fixed, ...)
 *   EASYBJETS_GAUSS_FIT       peak fit method: iterative (default) | windowshop
 *   EASYBJETS_FIT_MODEL       peak shape: gauss (default) | bukin | cb | dscb | all;
 *                             a non-default model writes to Run{2,3}Plots_<model>
 *   EASYBJETS_BKG_SELECT      background choice: cv (default) | dof
 *   EASYBJETS_CV_REGION       CV region: sidebands (default) | full
 *   EASYBJETS_SHOW_SPURIOUS   1 / on to print the spurious-signal test
 *   EASYBJETS_SIDEBAND1_LOW, EASYBJETS_SIDEBAND1_HIGH,
 *   EASYBJETS_SIDEBAND2_LOW, EASYBJETS_SIDEBAND2_MAX,
 *   EASYBJETS_MASS_XMIN, EASYBJETS_MASS_XMAX
 *                             override the sidebands and mass range (GeV)
 * EASYBJETS_MASS_CUT is not used by the fit, which needs the full mass range
 * (sidebands included).
 *
 * ==========================================================================
 * @version:  1.0
 * ==========================================================================
 */

#ifndef MASS_ANALYSIS_H
#define MASS_ANALYSIS_H

#include "VariantConfig.h"
#include "CommonUtils.h"
#include "BinningConfig.h"
#include "Math/QuantFuncMathCore.h"
#include "TGraphErrors.h"
#include "TLegendEntry.h"
#include "TStopwatch.h"
#include <algorithm>
#include <vector>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <map>
#include <string>
#include "Mass_SB_Plot.h"
#include "MassFit_Types.h"
#include "MassFit_Models.h"
#include "MassFit_Signal.h"
#include "MassFit_Plots.h"
#include "MassFit_Selection.h"

// ============================================================================
// MAIN ANALYSIS FUNCTION - WITH PT BINS
// ============================================================================
inline void MassIsolationExtractSig_Generic(const std::string& dataset, const std::string& wp,
                                            const PtBinConfig::PtBinDef& ptCfg) {
    VariantConfig vcfg = GetVariant(dataset);
    DatasetConfig config = DatasetConfig::FromVariant(vcfg);
    CurrentVariantLabel() = vcfg.label;

    std::string fileWP = PathBuilder::WPtoFileWP(wp);
    std::string labelWP_display = PathBuilder::WPtoLabel(wp);

    std::cout << "\n*******************************************" << std::endl;
    std::cout << "Processing " << config.name << " - WP: " << labelWP_display << std::endl;
    std::cout << "  Mass branch:  " << config.mass_branch << std::endl;
    std::cout << "  pT branch:    " << config.pt_branch << std::endl;
    std::cout << "  Mass range:   [" << config.xmin << ", " << config.xmax << "] GeV" << std::endl;
    std::cout << "  Sidebands:    [" << config.sideband1_low << ", " << config.sideband1_high << "] + ["
              << config.sideband2_low << ", " << config.sideband2_max << "] GeV" << std::endl;
    std::cout << "  Binning:      [" << config.min_bins << ", " << config.max_bins << "]" << std::endl;
    std::cout << "*******************************************\n" << std::endl;

    TStopwatch timer;
    timer.Start();

    std::string mcPath = PathBuilder::SamplePath(vcfg, SampleConfig::FileNames(vcfg)[SampleConfig::kZbb], wp);
    std::string dataPath = PathBuilder::SamplePath(vcfg, "data", wp);
    std::string multijetsPath =
        PathBuilder::SamplePath(vcfg, SampleConfig::FileNames(vcfg)[SampleConfig::kMultijets], wp);

    TFile* fmc = TFile::Open(mcPath.c_str(), "READ");
    TFile* fdata = TFile::Open(dataPath.c_str(), "READ");
    TFile* fmultijets = TFile::Open(multijetsPath.c_str(), "READ");

    if (!fmc || fmc->IsZombie() || !fdata || fdata->IsZombie() || !fmultijets || fmultijets->IsZombie()) {
        std::cerr << "Error: Could not open input files" << std::endl;
        if (fmc) delete fmc;
        if (fdata) delete fdata;
        if (fmultijets) delete fmultijets;
        return;
    }

    TTree* tree_mc = (TTree*)fmc->Get("nominal");
    TTree* tree_data = (TTree*)fdata->Get("nominal");
    TTree* tree_multijets = (TTree*)fmultijets->Get("nominal");

    if (!tree_mc || !tree_data || !tree_multijets) {
        std::cerr << "Error: Could not retrieve trees" << std::endl;
        fmc->Close();
        fdata->Close();
        fmultijets->Close();
        delete fmc;
        delete fdata;
        delete fmultijets;
        return;
    }

    // Set up branch addresses
    std::vector<float>* ljet_mass = nullptr;
    std::vector<float>* ljet_pt = nullptr;
    double weight, weight_multijets;
    int selected_ljet_index = -1;

    // MC: variant-specific branches (e.g. ljet_bJR10v00_mass, ljet_bJR10v00_pt)
    tree_mc->SetBranchAddress(config.mass_branch.c_str(), &ljet_mass);
    tree_mc->SetBranchAddress(config.pt_branch.c_str(), &ljet_pt);
    tree_mc->SetBranchAddress("weight", &weight);
    tree_mc->SetBranchAddress("selected_ljet_index", &selected_ljet_index);

    // Data falls back to the nominal branches if the variant ones are missing.
    bool dataHasVariantMass = Validate::HasBranch(tree_data, config.mass_branch);
    bool dataHasVariantPt = Validate::HasBranch(tree_data, config.pt_branch);
    tree_data->SetBranchAddress(dataHasVariantMass ? config.mass_branch.c_str() : "ljet_m", &ljet_mass);
    tree_data->SetBranchAddress(dataHasVariantPt ? config.pt_branch.c_str() : "ljet_pt", &ljet_pt);
    tree_data->SetBranchAddress("selected_ljet_index", &selected_ljet_index);
    std::cout << "  Data mass branch: " << (dataHasVariantMass ? config.mass_branch : "ljet_m [fallback]")
              << std::endl;
    std::cout << "  Data pT branch:   " << (dataHasVariantPt ? config.pt_branch : "ljet_pt [fallback]")
              << std::endl;

    // gamma+jets MC: variant-specific branches
    tree_multijets->SetBranchAddress(config.mass_branch.c_str(), &ljet_mass);
    tree_multijets->SetBranchAddress(config.pt_branch.c_str(), &ljet_pt);
    tree_multijets->SetBranchAddress("weight", &weight_multijets);
    tree_multijets->SetBranchAddress("selected_ljet_index", &selected_ljet_index);

    // ============================================================================
    // LOOP OVER PT BINS
    // ============================================================================
    for (int iPt = 0; iPt < ptCfg.nBins; ++iPt) {
        double ptMin = ptCfg.bins[iPt];
        double ptMax = ptCfg.bins[iPt + 1];

        std::cout << "\n=========================================" << std::endl;
        std::cout << "Processing pT bin " << (iPt + 1) << "/" << ptCfg.nBins << ": " << ptCfg.labels[iPt]
                  << std::endl;
        std::cout << "=========================================\n" << std::endl;

        // Create output directory for this pT bin. The mode tag is appended
        // so runs under different EASYBJETS_BINNING values don't overwrite
        // each other's PDFs - you can run all five modes in sequence and
        // the outputs sit side by side for comparison.
        std::string outputDirPath = (RunConfig::PlotDir() + "/") + config.output_subdir + "/Mass_Fit_"
                                    + CurrentChannel() + "/" + fileWP + "/" + ptCfg.fileNames[iPt] + "_"
                                    + BinningConfig::ModeTag() + "/";
        gSystem->Exec(("mkdir -p " + outputDirPath).c_str());

        // Calculate optimal binning for this pT bin using the selected mode
        // (fixed/fd/scott/sturges/peakaware/variable - see common/BinningConfig.h).
        std::cout << "Computing binning for pT bin " << ptCfg.labels[iPt] << "..." << std::endl;
        TH1F* h_temp_mc = new TH1F("h_temp_mc", "temp", 100, config.xmin, config.xmax);
        h_temp_mc->Sumw2();

        Long64_t nentries = tree_mc->GetEntries();
        for (Long64_t i = 0; i < nentries; i++) {
            tree_mc->GetEntry(i);
            if (selected_ljet_index >= 0 && selected_ljet_index < (int)ljet_mass->size()) {
                float jet_pt = (*ljet_pt)[selected_ljet_index];
                if (jet_pt >= ptMin && jet_pt < ptMax) {
                    h_temp_mc->Fill((*ljet_mass)[selected_ljet_index], weight);
                }
            }
        }

        BinningConfig::Result bin_result =
            BinningConfig::Compute(h_temp_mc, CurrentChannel(), config.xmin, config.xmax,
                                   //config.min_bins, config.max_bins
                                   config.GetMinBins(iPt), config.GetMaxBins(iPt));
        delete h_temp_mc;

        int optimal_bins = bin_result.nbins; // optimal number of bins
        std::cout << "Using " << optimal_bins << " bins for pT range " << ptCfg.labels[iPt] << " ["
                  << bin_result.caption << "]" << std::endl;

        // Create histograms for this pT bin. BinningConfig::MakeHist picks
        // between uniform and variable-width constructors based on whether
        // bin_result.edges is populated.
        std::string hist_suffix = "_" + std::string(ptCfg.fileNames[iPt]);
        TH1F* h_mc = BinningConfig::MakeHist("h_mc_mass" + hist_suffix, "Large-R Jet Mass", bin_result,
                                             config.xmin, config.xmax);
        TH1F* h_data = BinningConfig::MakeHist("h_data_mass" + hist_suffix, "Large-R Jet Mass", bin_result,
                                               config.xmin, config.xmax);
        TH1F* h_multijets = BinningConfig::MakeHist("h_multijets_mass" + hist_suffix, "Large-R Jet Mass",
                                                    bin_result, config.xmin, config.xmax);

        h_mc->Sumw2();
        h_multijets->Sumw2();

        // Set styles
        h_data->SetMarkerStyle(21);
        h_data->SetMarkerSize(0.8);
        h_data->SetMarkerColor(kBlack);
        h_data->SetLineColor(kBlack);
        h_multijets->SetFillColor(kBlue);
        h_multijets->SetLineColor(kBlue);
        h_mc->SetFillColor(kRed);
        h_mc->SetLineColor(kRed);

        // Fill histograms with pT selection
        std::cout << "Filling histograms for pT bin " << ptCfg.labels[iPt] << "..." << std::endl;

        // Fill MC
        nentries = tree_mc->GetEntries();
        for (Long64_t i = 0; i < nentries; i++) {
            tree_mc->GetEntry(i);
            if (selected_ljet_index >= 0 && selected_ljet_index < (int)ljet_mass->size()) {
                float jet_pt = (*ljet_pt)[selected_ljet_index];
                if (jet_pt >= ptMin && jet_pt < ptMax) {
                    h_mc->Fill((*ljet_mass)[selected_ljet_index], weight);
                }
            }
        }

        // Fill Data
        nentries = tree_data->GetEntries();
        for (Long64_t i = 0; i < nentries; i++) {
            tree_data->GetEntry(i);
            if (selected_ljet_index >= 0 && selected_ljet_index < (int)ljet_mass->size()) {
                float jet_pt = (*ljet_pt)[selected_ljet_index];
                if (jet_pt >= ptMin && jet_pt < ptMax) {
                    h_data->Fill((*ljet_mass)[selected_ljet_index]);
                }
            }
        }

        // Fill Multijets
        nentries = tree_multijets->GetEntries();
        for (Long64_t i = 0; i < nentries; i++) {
            tree_multijets->GetEntry(i);
            if (selected_ljet_index >= 0 && selected_ljet_index < (int)ljet_mass->size()) {
                float jet_pt = (*ljet_pt)[selected_ljet_index];
                if (jet_pt >= ptMin && jet_pt < ptMax) {
                    h_multijets->Fill((*ljet_mass)[selected_ljet_index], weight_multijets);
                }
            }
        }

        std::cout << "Events in this pT bin:" << std::endl;
        std::cout << "  MC: " << h_mc->Integral() << std::endl;
        std::cout << "  Data: " << h_data->Integral() << std::endl;
        std::cout << "  Multijets: " << h_multijets->Integral() << std::endl;

        // Check if we have enough statistics
        if (h_data->Integral() < 100) {
            std::cout << "WARNING: Low statistics in pT bin " << ptCfg.labels[iPt] << ". Skipping analysis."
                      << std::endl;
            delete h_mc;
            delete h_data;
            delete h_multijets;
            continue;
        }


        // Use configuration-specific sideband boundaries
        int binLow1 = h_multijets->FindBin(config.sideband1_low);
        int binHigh1 = h_multijets->FindBin(config.sideband1_high);
        int binLow2 = h_multijets->FindBin(config.sideband2_low);
        int binHigh2 = h_multijets->FindBin(config.sideband2_max);

        // Calculate scale factor for this pT bin
        std::vector<double> x_values, y_values, y_errors;

        for (int bin = binLow1; bin <= binHigh1; ++bin) {
            double bin_center = h_multijets->GetBinCenter(bin);
            double data_content = h_data->GetBinContent(bin);
            double mc_content = h_multijets->GetBinContent(bin);

            if (mc_content > 0) {
                x_values.push_back(bin_center);
                y_values.push_back(data_content / mc_content);
                double data_error = h_data->GetBinError(bin);
                double mc_error = h_multijets->GetBinError(bin);
                double ratio_error = sqrt(pow(data_error / mc_content, 2)
                                          + pow(data_content * mc_error / (mc_content * mc_content), 2));
                y_errors.push_back(ratio_error);
            }
        }

        for (int bin = binLow2; bin <= binHigh2; ++bin) {
            double bin_center = h_multijets->GetBinCenter(bin);
            double data_content = h_data->GetBinContent(bin);
            double mc_content = h_multijets->GetBinContent(bin);

            if (mc_content > 0) {
                x_values.push_back(bin_center);
                y_values.push_back(data_content / mc_content);
                double data_error = h_data->GetBinError(bin);
                double mc_error = h_multijets->GetBinError(bin);
                double ratio_error = sqrt(pow(data_error / mc_content, 2)
                                          + pow(data_content * mc_error / (mc_content * mc_content), 2));
                y_errors.push_back(ratio_error);
            }
        }

        TGraphErrors* ratio_graph = new TGraphErrors(x_values.size());
        for (size_t i = 0; i < x_values.size(); ++i) {
            ratio_graph->SetPoint(i, x_values[i], y_values[i]);
            ratio_graph->SetPointError(i, 0, y_errors[i]);
        }

        TF1* fit_func = new TF1("fit_func", "[0] + [1] * x", config.sideband1_low, config.sideband2_max);
        fit_func->SetParameters(1.0, 0.0);
        ratio_graph->Fit("fit_func", "RQ");

        double p0 = fit_func->GetParameter(0);
        double p1 = fit_func->GetParameter(1);
        double p0_err = fit_func->GetParError(0);
        double p1_err = fit_func->GetParError(1);

        // Save scale factor fit plot
        TCanvas* c_sf = new TCanvas("c_sf", "Scale Factor Fit", 800, 600);
        ratio_graph->SetTitle(Form("Data/Multijets Ratio - %s;M^{J} [GeV];Data/MC", ptCfg.labels[iPt]));
        ratio_graph->SetMarkerStyle(21);
        ratio_graph->SetMarkerColor(kBlue);
        ratio_graph->Draw("AP");
        fit_func->SetLineColor(kRed);
        fit_func->Draw("same");

        TLatex sf_label;
        sf_label.SetNDC();
        sf_label.SetTextSize(0.035);
        sf_label.DrawLatex(0.2, 0.8, Form("SF(m) = %.4f + %.6f #times m", p0, p1));
        sf_label.DrawLatex(0.2, 0.75,
                           Form("#chi^{2}/ndf = %.2f/%d", fit_func->GetChisquare(), fit_func->GetNDF()));
        AtlasUtils::DrawAtlasLabel(0.2, 0.9);

        SaveCanvasPdf(c_sf, outputDirPath + "scale_factor_fit.pdf");
        delete c_sf;

        // Scale multijets histogram
        TH1F* h_multijets_scaled = (TH1F*)h_multijets->Clone("h_multijets_scaled");
        for (int bin = 1; bin <= h_multijets_scaled->GetNbinsX(); ++bin) {
            double bin_center = h_multijets_scaled->GetBinCenter(bin);
            double content = h_multijets_scaled->GetBinContent(bin);
            double error = h_multijets_scaled->GetBinError(bin);
            double scale_factor = p0 + p1 * bin_center;
            h_multijets_scaled->SetBinContent(bin, content * scale_factor);
            h_multijets_scaled->SetBinError(bin, error * scale_factor);
        }

        // Create sideband histogram
        TH1F* h_data_sidebands = (TH1F*)h_data->Clone("h_data_sidebands");
        for (int i = binHigh1 + 1; i < binLow2; ++i) {
            h_data_sidebands->SetBinContent(i, 0);
            h_data_sidebands->SetBinError(i, 0);
        }

        // ============================================================================
        // PERFORM F-TEST BACKGROUND SELECTION AND FULL ANALYSIS
        // ============================================================================

        std::cout << "\nPerforming F-test background selection..." << std::endl;
        FitFunction* optimalFit = SelectOptimalBackgroundFit(
            h_data_sidebands, h_data, h_mc, h_multijets_scaled, config.xmin, config.xmax, binHigh1 + 1,
            binLow2 - 1, outputDirPath, ptCfg.labels[iPt]);


        TF1* polyExpFit = optimalFit->func;
        double chiSquare = polyExpFit->GetChisquare();
        int ndf = polyExpFit->GetNDF();
        double chi2NDF = (ndf > 0) ? chiSquare / ndf : 0.0;

        // Extract signal
        TH1F* h_signal = (TH1F*)h_data->Clone("h_signal");
        h_signal->SetTitle("Isolated Signal (Data);M^{J} [GeV];Events");


        for (int i = 1; i <= h_signal->GetNbinsX(); ++i) {
            double binCenter = h_signal->GetBinCenter(i);
            double binContent = h_signal->GetBinContent(i);
            double fitValue = std::max(0.0, polyExpFit->Eval(binCenter));
            h_signal->SetBinContent(i, binContent - fitValue);
        }

        // Summary plots - shared PadUtils ratio layout (same as Mass_Distributions)
        TCanvas* stack_canvas = new TCanvas("stack_canvas", "Mass Distribution", 800, 800);
        auto pads = PadUtils::CreateRatioPads(stack_canvas, 0.30);
        pads.upper->cd();

        THStack* hs = new THStack("hs", "Large-R Jet Mass;M^{J} [GeV];Events");
        hs->Add(h_multijets_scaled);
        hs->Add(h_mc);
        hs->SetMaximum(hs->GetMaximum() * 1.25);
        hs->Draw("HIST");
        hs->GetYaxis()->SetTitle("Events");
        hs->GetYaxis()->SetTitleSize(0.05);
        hs->GetXaxis()->SetLabelSize(0);

        double fitOrigXmin = polyExpFit->GetXmin();
        double fitOrigXmax = polyExpFit->GetXmax();
        double fitDrawMax = fitOrigXmax;
        for (double x = fitOrigXmin; x <= fitOrigXmax; x += 0.5) {
            if (polyExpFit->Eval(x) < 0) {
                fitDrawMax = x;
                break;
            }
        }
        polyExpFit->SetRange(fitOrigXmin, fitDrawMax);
        polyExpFit->SetLineColor(kGreen);
        polyExpFit->SetLineWidth(2);
        polyExpFit->Draw("same");
        polyExpFit->SetRange(fitOrigXmin, fitOrigXmax);
        HistUtils::StyleDataPoints(h_data);
        h_data->Draw("EP same");

        TH1F* h_mc_total = (TH1F*)hs->GetStack()->Last();
        TGraphAsymmErrors* err_band = HistUtils::MakeUncertaintyBand(h_mc_total);
        err_band->Draw("2 same");

        TLegend* legend = new TLegend(0.55, 0.55, 0.9, 0.9);
        legend->SetBorderSize(0);
        legend->SetTextSize(0.028);
        legend->AddEntry(h_data, "Data", "lep");
        legend->AddEntry(h_mc, SampleConfig::DisplayNames(CurrentChannel())[SampleConfig::kZbb].c_str(), "f");
        legend->AddEntry(h_multijets_scaled,
                         SampleConfig::DisplayNames(CurrentChannel())[SampleConfig::kMultijets].c_str(), "f");
        legend->AddEntry(err_band, "MC Stat. Unc.", "f");

        double safe_chi2ndf = (ndf > 0 && std::isfinite(chiSquare)) ? chiSquare / ndf : 0.0;
        std::string fit_entry = optimalFit->name + " Fit";
        if (safe_chi2ndf > 0 && safe_chi2ndf < 1000) {
            char chi2_str[100];
            snprintf(chi2_str, sizeof(chi2_str), " (#chi^{2}/NDF = %.2f)", safe_chi2ndf);
            fit_entry += chi2_str;
        }
        legend->AddEntry(polyExpFit, fit_entry.c_str(), "l");
        legend->Draw();

        // Add pT bin label
        TLatex ptLabel;
        ptLabel.SetNDC();
        ptLabel.SetTextSize(0.035);
        ptLabel.DrawLatex(0.2, 0.80, ptCfg.labels[iPt]);

        // Binning mode caption - makes it possible to tell which binning
        // produced this plot at a glance.
        TLatex binLabel;
        binLabel.SetNDC();
        binLabel.SetTextSize(0.030);
        binLabel.SetTextColor(kGray + 2);
        binLabel.DrawLatex(0.2, 0.76, bin_result.caption.c_str());

        AtlasUtils::DrawAtlasLabel(0.2, 0.9);

        // Ratio panel - shared style (large, readable Data/MC labels).
        pads.lower->cd();
        pads.lower->SetGridy(0);
        TH1F* h_ratio = (TH1F*)h_data->Clone("h_ratio");
        h_ratio->SetDirectory(0);
        h_ratio->SetTitle("");
        h_ratio->Divide(h_mc_total);
        PadUtils::StyleRatioHist(h_ratio, "M^{J} [GeV]", 0.5, 1.5);
        HistUtils::StyleDataPoints(h_ratio);
        h_ratio->Draw("EP");

        TGraphAsymmErrors* ratio_err_band = HistUtils::MakeRatioUncertaintyBand(h_mc_total);
        ratio_err_band->Draw("2 same");
        h_ratio->Draw("EP same");

        TLine* line = PadUtils::DrawRefLine(pads.lower, 1.0);

        SaveCanvasPdf(stack_canvas, outputDirPath + "Mass_stack_scaled.pdf");
        delete h_ratio;
        delete line;
        delete err_band;
        delete ratio_err_band;
        delete legend;
        delete stack_canvas;

        // Fit Gaussians to signal
        std::cout << "\nFitting Gaussians to extracted signal..." << std::endl;
        TH1F* h_mc_copy = (TH1F*)h_mc->Clone("h_mc_for_fit");
        h_mc_copy->SetDirectory(0);

        TF1* mc_gauss_fit = PerformOptimizedGaussianFit(h_mc_copy, "mc_gauss_fit");
        TF1* data_gauss_fit =
            PerformOptimizedGaussianFit(h_signal, "data_gauss_fit", {0, 0}, true, outputDirPath);

        // Create signal comparison plot
        TCanvas* c_signal = new TCanvas("c_signal", "Signal Comparison", 800, 600);

        double mc_max = h_mc->GetBinContent(h_mc->GetMaximumBin());
        double signal_max = h_signal->GetBinContent(h_signal->GetMaximumBin());
        double overall_max = std::max(mc_max, signal_max);

        h_mc->SetFillColor(kRed);
        h_mc->GetYaxis()->SetRangeUser(0, overall_max * 1.3);
        h_mc->GetXaxis()->SetTitle("M^{J} [GeV]");
        h_mc->GetYaxis()->SetTitle("Events");
        h_mc->Draw("HIST");

        h_signal->SetMarkerStyle(21);
        h_signal->SetMarkerColor(kBlack);
        h_signal->SetLineColor(kBlack);
        h_signal->Draw("EP same");

        mc_gauss_fit->SetLineColor(kCyan);
        mc_gauss_fit->SetLineWidth(2);
        mc_gauss_fit->Draw("same");

        data_gauss_fit->SetLineColor(kViolet - 4);
        data_gauss_fit->SetLineWidth(2);
        data_gauss_fit->Draw("same");

        TLegend* leg_signal = new TLegend(0.50, 0.65, 0.95, 0.9);
        leg_signal->SetBorderSize(0);
        leg_signal->SetTextSize(0.025);
        leg_signal->AddEntry(h_signal, "Isolated Signal (Data)", "lep");
        leg_signal->AddEntry(h_mc, "MC Signal", "f");

        double mc_mean = mc_gauss_fit->GetParameter(1);
        double mc_sigma = mc_gauss_fit->GetParameter(2);
        double mc_chi2_ndf =
            (mc_gauss_fit->GetNDF() > 0) ? mc_gauss_fit->GetChisquare() / mc_gauss_fit->GetNDF() : 0;

        double data_mean = data_gauss_fit->GetParameter(1);
        double data_sigma = data_gauss_fit->GetParameter(2);
        double data_chi2_ndf =
            (data_gauss_fit->GetNDF() > 0) ? data_gauss_fit->GetChisquare() / data_gauss_fit->GetNDF() : 0;

        if (std::isfinite(mc_mean) && std::isfinite(mc_sigma) && mc_mean > 0 && mc_sigma > 0) {
            leg_signal->AddEntry(
                mc_gauss_fit,
                Form("MC: #mu = %.2f, #sigma = %.2f, #chi^{2}/ndf = %.2f", mc_mean, mc_sigma, mc_chi2_ndf),
                "l");
        }

        if (std::isfinite(data_mean) && std::isfinite(data_sigma) && data_mean > 0 && data_sigma > 0) {
            leg_signal->AddEntry(data_gauss_fit,
                                 Form("Data: #mu = %.2f, #sigma = %.2f, #chi^{2}/ndf = %.2f", data_mean,
                                      data_sigma, data_chi2_ndf),
                                 "l");
        }
        leg_signal->Draw();

        AtlasUtils::DrawAtlasLabel();
        TLatex pt_label2;
        pt_label2.SetNDC();
        pt_label2.SetTextSize(0.035);
        pt_label2.DrawLatex(0.18, 0.76, ptCfg.labels[iPt]);

        SaveCanvasPdf(c_signal, outputDirPath + "Isolated_Signal_GaussFits.pdf");
        delete c_signal;
        delete leg_signal;

        // Save summary
        double estimateMJBackground = h_multijets_scaled->Integral(binHigh1 + 1, binLow2 - 1);
        double dataSignal = h_data->Integral(binHigh1 + 1, binLow2 - 1);
        double signalEvents = dataSignal - estimateMJBackground;

        std::ofstream summary((outputDirPath + "fit_summary.txt").c_str());
        if (summary.is_open()) {
            summary << "F-Test Mass Analysis Summary" << std::endl;
            summary << "pT Bin: " << ptCfg.labels[iPt] << std::endl;
            summary << "WP: " << labelWP_display << std::endl;
            summary << "=================================" << std::endl << std::endl;
            summary << "Events in pT bin:" << std::endl;
            summary << "  MC: " << h_mc->Integral() << std::endl;
            summary << "  Data: " << h_data->Integral() << std::endl;
            summary << "  Multijets: " << h_multijets->Integral() << std::endl << std::endl;

            summary << "DATA/MC FIT PARAMETERS:" << std::endl;
            summary << "  Function: SF(m) = p0 + p1 * m" << std::endl;
            summary << "  Fit Range: [" << std::fixed << std::setprecision(2) << config.sideband1_low << ", "
                    << config.sideband2_max << "] GeV" << std::endl;
            summary << "  Number of Points: " << ratio_graph->GetN() << std::endl;
            summary << "  Parameters:" << std::endl;
            summary << "    p0 = " << std::setprecision(6) << p0 << " +/- " << fit_func->GetParError(0)
                    << std::endl;
            summary << "    p1 = " << std::scientific << std::setprecision(6) << p1 << " +/- "
                    << fit_func->GetParError(1) << " GeV^-1" << std::endl;
            summary << "  Fit Quality:" << std::endl;
            summary << "    Chi^2 = " << std::fixed << std::setprecision(4) << fit_func->GetChisquare()
                    << std::endl;
            summary << "    NDF = " << fit_func->GetNDF() << std::endl;
            summary << "    Chi^2/NDF = " << std::setprecision(4)
                    << fit_func->GetChisquare() / fit_func->GetNDF() << std::endl;
            summary << "    Prob(Chi^2,NDF) = " << std::setprecision(4) << fit_func->GetProb() << std::endl
                    << std::endl;
            summary << "OPTIMAL BACKGROUND FUNCTION:" << std::endl;
            summary << "  Function: " << optimalFit->name << std::endl;
            summary << "  Chi^2/ndf = " << std::fixed << std::setprecision(4) << chi2NDF << std::endl;
            summary << "  Spurious signal = " << optimalFit->spurious_signal << " events" << std::endl
                    << std::endl;
            summary << "Signal Extraction:" << std::endl;
            summary << "  Data Events in Signal Region: " << dataSignal << std::endl;
            summary << "  Estimated Background: " << estimateMJBackground << std::endl;
            summary << "  Extracted Signal Events: " << signalEvents << std::endl << std::endl;
            summary << "Gaussian Fits:" << std::endl;
            summary << "  MC: mu = " << std::setprecision(3) << mc_gauss_fit->GetParameter(1) << " +/- "
                    << mc_gauss_fit->GetParError(1) << " GeV" << std::endl;
            summary << "      sigma = " << mc_gauss_fit->GetParameter(2) << " +/- "
                    << mc_gauss_fit->GetParError(2) << " GeV" << std::endl;
            summary << "  Data: mu = " << data_gauss_fit->GetParameter(1) << " +/- "
                    << data_gauss_fit->GetParError(1) << " GeV" << std::endl;
            summary << "        sigma = " << data_gauss_fit->GetParameter(2) << " +/- "
                    << data_gauss_fit->GetParError(2) << " GeV" << std::endl;
            summary << "  Mass Shift (Data-MC): Deltamu = "
                    << (data_gauss_fit->GetParameter(1) - mc_gauss_fit->GetParameter(1)) << " GeV"
                    << std::endl;

            double sigma_data = data_gauss_fit->GetParameter(2);
            double sigma_mc = mc_gauss_fit->GetParameter(2);
            double sigma_diff = sigma_data - sigma_mc;
            double sigma_quadrature = TMath::Sqrt(TMath::Abs(sigma_data * sigma_data - sigma_mc * sigma_mc));
            // JMR = (sigma_data − sigma_MC) × √|sigma^2_data − sigma^2_MC|
            double jmr = sigma_diff * sigma_quadrature;

            summary << "  Data/MC Ratios:" << std::endl;
            summary << "    mu_Data/mu_MC = " << std::setprecision(6)
                    << (data_gauss_fit->GetParameter(1) / mc_gauss_fit->GetParameter(1)) << std::endl;
            summary << "    sigma_Data/sigma_MC = " << std::setprecision(6) << (sigma_data / sigma_mc)
                    << std::endl;
            summary << "    Deltasigma = " << std::setprecision(4) << sigma_diff << " GeV" << std::endl;
            summary << "    JMR = (sigma_data−sigma_MC)×√|sigma^2_data−sigma^2_MC| = " << std::setprecision(4)
                    << jmr << " GeV^2" << std::endl;
            summary.close();
        }

        std::cout << "\n=========================================" << std::endl;
        std::cout << "Summary for pT bin " << ptCfg.labels[iPt] << std::endl;
        std::cout << "=========================================" << std::endl;
        std::cout << "Optimal background: " << optimalFit->name << std::endl;
        std::cout << "Background Chi^2/ndf = " << chi2NDF << std::endl;
        std::cout << "Signal events: " << signalEvents << std::endl;
        std::cout << "=========================================\n" << std::endl;

        if (h_signal) {
            h_signal->SetDirectory(0);
            h_signal->GetListOfFunctions()->Clear();
            delete h_signal;
            h_signal = nullptr;
        }
        if (h_mc_copy) {
            h_mc_copy->SetDirectory(0);
            delete h_mc_copy;
            h_mc_copy = nullptr;
        }
        if (mc_gauss_fit) {
            delete mc_gauss_fit;
            mc_gauss_fit = nullptr;
        }
        if (data_gauss_fit) {
            delete data_gauss_fit;
            data_gauss_fit = nullptr;
        }
        if (optimalFit) {
            if (optimalFit->func) delete optimalFit->func;
            delete optimalFit;
            optimalFit = nullptr;
        }

        h_mc->SetDirectory(0);
        h_data->SetDirectory(0);
        h_multijets->SetDirectory(0);
        h_multijets_scaled->SetDirectory(0);
        h_data_sidebands->SetDirectory(0);

        delete h_mc;
        delete h_data;
        delete h_multijets;
        delete h_multijets_scaled;
        delete h_data_sidebands;
        delete ratio_graph;
        delete fit_func;

        std::cout << "Completed pT bin " << ptCfg.labels[iPt] << std::endl;
    }


    fmc->Close();
    fdata->Close();
    fmultijets->Close();
    delete fmc;
    delete fdata;
    delete fmultijets;

    timer.Stop();
    std::cout << "\nCompleted all pT bins for " << config.name << " WP " << labelWP_display << " in ";
    timer.Print();
}

inline void Mass_Analysis_Run(const std::string& variant, const std::string& ptMode,
                              const std::string& channel) {
    CurrentChannel() = channel;
    gROOT->SetBatch(kTRUE);
    AtlasUtils::LoadAtlasStyle();

    // Select pT binning from VariantConfig.h
    auto ptCfg = PtBinConfig::GetBins(ptMode, CurrentChannel());
    Log::info("MassAnalysis",
              Form("pT mode: %s (%d bin%s)", ptMode.c_str(), ptCfg.nBins, ptCfg.nBins > 1 ? "s" : ""));

    auto workingPoints = AnalysisWPConfig::StandardWPs();

    std::vector<std::string> datasets;
    if (variant == "ALL") {
        datasets = GetPrimaryVariantKeys(CurrentChannel());
    } else {
        datasets.push_back(variant);
    }

    for (const auto& dataset : datasets) {
        for (const auto& wp : workingPoints) {
            std::cout << "\n***********************************************" << std::endl;
            std::cout << "Starting analysis: " << dataset << " - WP " << wp << "%" << std::endl;
            std::cout << "pT mode: " << ptMode << " (" << ptCfg.nBins << " bins)" << std::endl;
            std::cout << "***********************************************\n" << std::endl;

            MassIsolationExtractSig_Generic(dataset, wp, ptCfg);
            if (ptCfg.nBins > 1) CreatePtBinComparisonPlots(dataset, wp, ptCfg);
        }
    }

    std::cout << "\nAll datasets, working points, and pT bins processed!" << std::endl;
}

#endif // MASS_ANALYSIS_H