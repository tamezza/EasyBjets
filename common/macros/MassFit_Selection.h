/**
 * ============================================================================
 *  MassFit_Selection.h - choice of the background function
 * ============================================================================
 * ===========================================
 * @author:    Thiziri AMEZZA
 * @email:     thiziri.amezza@cern.ch
 * @date:      18 April 2026
 * ===========================================
 *
 * SelectOptimalBackgroundFit() does, for one variant, WP and pT bin:
 *   1. fit the 12 candidates (poly, exp, polyexp; orders 2 to 5) to the data
 *      sidebands;
 *   2. within each family, an F-test between consecutive orders keeps the lowest
 *      order that the next one does not improve significantly (p < 0.05); this
 *      gives one finalist per family;
 *   3. fit the finalists to the MC background template as well, and run a k-fold
 *      cross-validation (k = 5) of every candidate on that template: the usable
 *      bins are split into k interleaved folds, each fold is predicted by a fit
 *      to the other k-1, and the score is the mean held-out chi2 per bin with its
 *      standard error over the folds;
 *   4. among the finalists within one standard error of the best CV score, pick
 *      the one with the fewest parameters (ties: better average of the data and
 *      MC chi2/ndf).
 * If the CV is not possible (too few bins), the lowest combined data and MC
 * chi2/ndf is used (the CV is still computed and plotted).
 *
 * Written to the output folder: ftest_results.txt, ftest_all_results.txt,
 * all_fit_parameters.txt, optimal_fit_parameters.txt, the two comparison
 * canvases, and cross_validation/ with the CV tables and plots. The per-function
 * plots and signal_fit_params.txt come from CreateIndividualFitAnalysis
 * (MassFit_Plots.h), called for each finalist.
 *
 * Part of the mass fit: it is not run on its own but through Mass_Fit.cxx
 * (how to run and all options: Mass_Fit.h).
 *
 * Options read here:
 *   EASYBJETS_BKG_SELECT   cv (default): the CV rule above
 *                          dof: the finalist with the fewest parameters, no CV
 *   EASYBJETS_CV_REGION    sidebands (default): CV on the sidebands of the MC
 *                          template, the same region as the F-test
 *                          full: CV on the whole mass range (cross-check)
 *   EASYBJETS_SHOW_SPURIOUS   1 / on to print the spurious-signal test
 *
 * ==========================================================================
 * @version:  1.0
 * ==========================================================================
 */

#ifndef MASSFIT_SELECTION_H
#define MASSFIT_SELECTION_H

#include "MassFit_Types.h"
#include "MassFit_Models.h"
#include "MassFit_Plots.h"
#include "TLegend.h"

// F statistic for the nested pair (np parameters, chi2_p) vs (nq, chi2_q) on n bins.
inline double CalculateFTest(double chi2_p, double chi2_q, int np, int nq, int n) {
    int ndf_p = n - np;
    int ndf_q = n - nq;
    if (ndf_q <= 0 || ndf_p <= 0 || ndf_p <= ndf_q) return -1;
    double F = ((chi2_p - chi2_q) / (ndf_p - ndf_q)) / (chi2_q / ndf_q);
    return F;
}

inline double GetFCritical(int df1, int df2, double alpha = 0.05) {
    return ROOT::Math::fdistribution_quantile(1.0 - alpha, df1, df2);
}

// Events left in the signal window after subtracting the background function from the data.
// This is a diagnostic, not used in the selection (old checks not used in our analysis now,
// left here in case we want to revisit this).
inline double CalculateSpuriousSignal(TH1F* h_data, TH1F* h_mc_template, TF1* fitFunc, int binSignalLow,
                                      int binSignalHigh) {
    TH1F* h_subtracted = (TH1F*)h_data->Clone("h_temp_subtracted");

    for (int i = 1; i <= h_subtracted->GetNbinsX(); ++i) {
        double binCenter = h_subtracted->GetBinCenter(i);
        double binContent = h_subtracted->GetBinContent(i);
        double fitValue = std::max(0.0, fitFunc->Eval(binCenter));
        h_subtracted->SetBinContent(i, binContent - fitValue);
        if (h_subtracted->GetBinContent(i) < 0) {
            h_subtracted->SetBinContent(i, 0);
        }
    }

    double spurious = h_subtracted->Integral(binSignalLow, binSignalHigh);
    delete h_subtracted;
    return spurious;
}

// ============================================================================
// k-FOLD CROSS-VALIDATION
// ----------------------------------------------------------------------------
// Evaluates how well a candidate background function predicts unseen bins,
// rather than how well it fits the data used for training. This helps identify
// overfitting that may not be apparent from the in-sample chi2 or F-test.
//
// The usable bins of `h` in [xmin, xmax] are split into k interleaved folds.
// For each fold, the function is fitted to the other k-1 folds and evaluated
// on the held-out bins. Interleaved folds cover the full mass range, making
// each test an interpolation rather than an extrapolation test.
//
// Each training fit uses a TGraphErrors containing only the training bins,
// starting from the parameter values of the converged full-range fit.
//
// CVResult is defined in MassFit_Types.h.
// ============================================================================
// Fits the function on k-1 folds and evaluates its chi2 on the held-out fold.

inline CVResult KFoldCV(TF1* templateFunc, TH1F* h, double xmin, double xmax, int k = 5) {
    CVResult R;
    if (!templateFunc || !h) return R;

    // Collect bins within the fit range that have a non-zero uncertainty.
    std::vector<int> bins;
    for (int b = 1; b <= h->GetNbinsX(); ++b) {
        double xc = h->GetBinCenter(b);
        if (xc < xmin || xc > xmax) continue;
        if (h->GetBinError(b) <= 0) continue;
        bins.push_back(b);
    }
    const int N = (int)bins.size();
    const int npar = templateFunc->GetNpar();
    // Require enough training bins to constrain the fit parameters.
    if (k < 2 || N < k || (N - N / k) < npar + 1) return R;

    std::vector<double> foldChi2(k, 0.0);
    std::vector<int> foldN(k, 0);

    for (int f = 0; f < k; ++f) {
        // Build the training graph, excluding the current fold.
        TGraphErrors gr;
        int p = 0;
        for (int idx = 0; idx < N; ++idx) {
            if (idx % k == f) continue;
            int b = bins[idx];
            gr.SetPoint(p, h->GetBinCenter(b), h->GetBinContent(b));
            gr.SetPointError(p, 0.0, h->GetBinError(b));
            ++p;
        }
        // Start from the parameters of the full-range fit.
        TF1* ff = (TF1*)templateFunc->Clone(Form("%s_cv%d", templateFunc->GetName(), f));
        for (int ip = 0; ip < ff->GetNpar(); ++ip) ff->SetParameter(ip, templateFunc->GetParameter(ip));
        gr.Fit(ff, "Q0");

        // Evaluate chi2 on the held-out bins.
        double chi2t = 0.0;
        int nt = 0;
        for (int idx = 0; idx < N; ++idx) {
            if (idx % k != f) continue;
            int b = bins[idx];
            double e = h->GetBinError(b);
            double r = (h->GetBinContent(b) - ff->Eval(h->GetBinCenter(b))) / e;
            chi2t += r * r;
            ++nt;
        }
        foldChi2[f] = chi2t;
        foldN[f] = nt;
        delete ff;
    }

    double totChi2 = 0;
    int totN = 0;
    std::vector<double> perFold(k, 0.0);
    for (int f = 0; f < k; ++f) {
        totChi2 += foldChi2[f];
        totN += foldN[f];
        perFold[f] = (foldN[f] > 0) ? foldChi2[f] / foldN[f] : 0.0;
    }
    if (totN == 0) return R;
    // Use the mean chi2 per bin across folds as the CV score, with the
    // standard error calculated from the same per-fold values. This keeps
    // the score and its uncertainty consistent when fold sizes differ.
    double mean = 0;
    for (double v : perFold) mean += v;
    mean /= k;
    double var = 0;
    for (double v : perFold) var += (v - mean) * (v - mean);
    var /= (k > 1 ? (k - 1) : 1);
    R.score = mean;
    R.sem = std::sqrt(var / k);
    R.nTest = totN;
    R.valid = true;
    return R;
}

// ============================================================================
// CROSS-VALIDATION DIAGNOSTICS & MODEL SELECTION
// ----------------------------------------------------------------------------
// Evaluates candidate models on held-out folds to assess predictive performance
// and help identify overfitting. Computes per-fold metrics, generates
// diagnostic plots, and applies the 1-standard-error (1-SE) rule for selection.
//
// Selection procedure:
//  1. Fit each model on k-1 folds and evaluate it on the held-out fold.
//  2. Compare models using their out-of-sample (OOS) chi2 per bin.
//  3. Among models within one standard error of the minimum CV score, prefer
//     the model with the fewest parameters.
//
// This provides a complementary measure of model quality to the in-sample
// fit and F-tests, helping balance predictive performance and model complexity.
// ============================================================================
struct CVFoldRow {
    int nTrain = 0, nTest = 0;
    double chi2 = 0, chi2perbin = 0;
};

struct CVFoldsDetail {
    std::vector<CVFoldRow> folds;
    double score = 1e9, sem = 0;
    bool valid = false;
};

inline CVFoldsDetail KFoldCVDetailed(TF1* templateFunc, TH1F* h, double xmin, double xmax, int k = 5) {
    CVFoldsDetail D;
    if (!templateFunc || !h) return D;
    std::vector<int> bins;
    for (int b = 1; b <= h->GetNbinsX(); ++b) {
        double xc = h->GetBinCenter(b);
        if (xc < xmin || xc > xmax) continue;
        if (h->GetBinError(b) <= 0) continue;
        bins.push_back(b);
    }
    const int N = (int)bins.size(), npar = templateFunc->GetNpar();
    if (k < 2 || N < k || (N - N / k) < npar + 1) return D;
    std::vector<double> perFold(k, 0.0);
    double totChi2 = 0;
    int totN = 0;
    for (int f = 0; f < k; ++f) {
        TGraphErrors gr;
        int p = 0, ntr = 0;
        for (int idx = 0; idx < N; ++idx) {
            if (idx % k == f) continue;
            int b = bins[idx];
            gr.SetPoint(p, h->GetBinCenter(b), h->GetBinContent(b));
            gr.SetPointError(p, 0.0, h->GetBinError(b));
            ++p;
            ++ntr;
        }
        TF1* ff = (TF1*)templateFunc->Clone(Form("%s_cvd%d", templateFunc->GetName(), f));
        for (int ip = 0; ip < ff->GetNpar(); ++ip) ff->SetParameter(ip, templateFunc->GetParameter(ip));
        gr.Fit(ff, "Q0");
        double chi2t = 0;
        int nt = 0;
        for (int idx = 0; idx < N; ++idx) {
            if (idx % k != f) continue;
            int b = bins[idx];
            double e = h->GetBinError(b);
            double r = (h->GetBinContent(b) - ff->Eval(h->GetBinCenter(b))) / e;
            chi2t += r * r;
            ++nt;
        }
        CVFoldRow row;
        row.nTrain = ntr;
        row.nTest = nt;
        row.chi2 = chi2t;
        row.chi2perbin = (nt > 0 ? chi2t / nt : 0.0);
        D.folds.push_back(row);
        perFold[f] = row.chi2perbin;
        totChi2 += chi2t;
        totN += nt;
        delete ff;
    }
    if (totN == 0) return D;
    double mean = 0;
    for (double v : perFold) mean += v;
    mean /= k;
    double var = 0;
    for (double v : perFold) var += (v - mean) * (v - mean);
    var /= (k > 1 ? (k - 1) : 1);
    D.score = mean;
    D.sem = std::sqrt(var / k);
    D.valid = true;
    return D;
}

// Returns the valid candidate with the lowest in-sample chi2/ndf.
// Used to compare the in-sample ranking with the cross-validation result.
// what a selection without CV would pick, shown on the decision plot when
// it differs from the choice.
inline const FitFunction* InSampleChampion(const std::vector<FitFunction>& allFits) {
    const FitFunction* champ = nullptr;
    for (const auto& fn : allFits)
        if (fn.chi2_ndf < 999 && (!champ || fn.chi2_ndf < champ->chi2_ndf)) champ = &fn;
    return champ;
}

// Writes the per-fold results, the 1-SE selection and a comparison between
// the in-sample best fit and the final selection.
inline void WriteCVFoldsTxt(const std::string& path, const std::vector<FitFunction*>& finalists,
                            const FitFunction* chosen, const std::vector<FitFunction>& allFits, TH1F* h,
                            double xmin, double xmax, int k, double cvMin, double cvThreshold) {
    std::ofstream o(path);
    if (!o.is_open()) return;
    o << "k-fold cross-validation, per-fold breakdown (k=" << k << ")\n";
    o << "Each fold is predicted by a fit to the other k-1 folds and scored by its held-out chi2\n";
    o << "on the held-out fold (out-of-sample). mean = (1/k) sum(fold chi2/bin),\n";
    o << "SE = sqrt(var/k). The 1-SE rule below turns these into one choice.\n";
    o << "================================================================\n\n";
    for (auto* ff : finalists) {
        CVFoldsDetail d = KFoldCVDetailed(ff->func, h, xmin, xmax, k);
        o << ff->name << "   (npar=" << ff->npar << ", in-sample chi2/ndf=" << std::setprecision(4)
          << ff->chi2_ndf << ")" << (ff == chosen ? "   <<< CHOSEN" : "") << "\n";
        if (!d.valid) {
            o << "   [CV unavailable - too few bins]\n\n";
            continue;
        }
        o << "   fold | #train | #test | held-out chi2 | chi2/bin\n";
        o << "   -----+--------+-------+---------------+---------\n";
        for (size_t f = 0; f < d.folds.size(); ++f)
            o << "   " << std::setw(4) << (f + 1) << " | " << std::setw(6) << d.folds[f].nTrain << " | "
              << std::setw(5) << d.folds[f].nTest << " | " << std::setw(13) << std::setprecision(4)
              << d.folds[f].chi2 << " | " << std::setprecision(4) << d.folds[f].chi2perbin << "\n";
        o << "   --> CV mean = " << std::setprecision(4) << d.score << "    SE = " << std::setprecision(3)
          << d.sem << "\n\n";
    }
    o << "----------------------------------------------------------------\n";
    o << "1-SE RULE -> SINGLE CHOICE\n";
    o << "  best (lowest) CV among finalists = " << std::setprecision(4) << cvMin << "\n";
    o << "  1-SE acceptance band = [" << std::setprecision(4) << cvMin << ", " << cvThreshold << "]\n";
    o << "  finalists inside the band (statistically equivalent out-of-sample):\n";
    for (auto* ff : finalists)
        if (ff->cv_valid && ff->cv_score <= cvThreshold)
            o << "      " << ff->name << "  (CV=" << std::setprecision(4) << ff->cv_score
              << ", npar=" << ff->npar << ")\n";
    o << "  -> keep the FEWEST-parameter one (ties: best combined data+MC chi2/ndf): "
      << (chosen ? chosen->name : "?") << "\n\n";
    o << "----------------------------------------------------------------\n";
    o << "OVERFITTING CHECK (does in-sample disagree with out-of-sample?)\n";
    const FitFunction* champ = InSampleChampion(allFits);
    if (champ)
        o << "  in-sample champion (lowest chi2/ndf - what a greedy in-sample / over-eager\n"
          << "    F-test rule gravitates to): " << champ->name << "  (npar=" << champ->npar
          << ", in-sample chi2/ndf=" << std::setprecision(4) << champ->chi2_ndf
          << (champ->cv_valid ? (", CV=" + std::to_string(champ->cv_score)) : "") << ")\n";
    if (chosen)
        o << "  selected (CV + 1-SE):         " << chosen->name << "  (npar=" << chosen->npar
          << ", in-sample chi2/ndf=" << std::setprecision(4) << chosen->chi2_ndf
          << (chosen->cv_valid ? (", CV=" + std::to_string(chosen->cv_score)) : "") << ")\n";
    if (champ && chosen) {
        if (champ == chosen)
            o << "  => No demotion: the in-sample-best function is also the simplest within\n"
              << "     the 1-SE CV band.\n";
        else if (chosen->npar <= champ->npar && champ->cv_valid && chosen->cv_valid
                 && chosen->cv_score <= champ->cv_score + 1e-9)
            o << "  => OVERFITTING AVOIDED: " << champ->name << " has the better IN-sample fit ("
              << std::setprecision(4) << champ->chi2_ndf << ") but its OUT-of-sample CV (" << champ->cv_score
              << ") is NOT better than " << chosen->name << " (CV=" << chosen->cv_score << "), while costing "
              << (champ->npar - chosen->npar) << " more parameter(s). The 1-SE rule keeps the simpler "
              << chosen->name << ".\n";
        else
            o << "  => Selected differs from the in-sample champion; see the CV band above.\n";
    }
    o.close();
}

// Draws one fold on the current pad: training bins (black), held-out bins
// (red, open), the fit to the training bins, and the held-out chi2/bin.
// Saves both a combined figure and individual plots for each fold.
inline void DrawCVFoldsAll(TF1* templateFunc, TH1F* h, double xmin, double xmax, int k,
                           const std::string& outPath, const std::string& titleTag) {
    if (!templateFunc || !h) return;
    (void)titleTag;
    std::vector<int> bins;
    for (int b = 1; b <= h->GetNbinsX(); ++b) {
        double xc = h->GetBinCenter(b);
        if (xc < xmin || xc > xmax) continue;
        if (h->GetBinError(b) <= 0) continue;
        bins.push_back(b);
    }
    const int N = (int)bins.size();
    if (N < k) return;
    double ymax = 0;
    for (int idx = 0; idx < N; ++idx) {
        int b = bins[idx];
        ymax = std::max(ymax, h->GetBinContent(b) + h->GetBinError(b));
    }
    const std::string fname = templateFunc->GetName();

    // Fit each fold once and store the training points, held-out points,
    // fitted parameters and held-out chi2/bin for plotting.
    struct FoldData {
        std::vector<double> trx, trY, trE, tex, teY, teE;
        std::vector<double> par;
        double chi2bin = 0;
    };
    std::vector<FoldData> folds(k);
    for (int f = 0; f < k; ++f) {
        FoldData& fd = folds[f];
        TGraphErrors grfit;
        int pf = 0;
        for (int idx = 0; idx < N; ++idx) {
            int b = bins[idx];
            double x = h->GetBinCenter(b), y = h->GetBinContent(b), e = h->GetBinError(b);
            if (idx % k == f) {
                fd.tex.push_back(x);
                fd.teY.push_back(y);
                fd.teE.push_back(e);
            } else {
                fd.trx.push_back(x);
                fd.trY.push_back(y);
                fd.trE.push_back(e);
                grfit.SetPoint(pf, x, y);
                grfit.SetPointError(pf, 0, e);
                ++pf;
            }
        }
        TF1* ff = (TF1*)templateFunc->Clone(Form("cvfitf%d", f));
        for (int ip = 0; ip < ff->GetNpar(); ++ip) ff->SetParameter(ip, templateFunc->GetParameter(ip));
        grfit.Fit(ff, "Q0");
        double chi2t = 0;
        int nt = 0;
        for (size_t j = 0; j < fd.tex.size(); ++j) {
            double r = (fd.teY[j] - ff->Eval(fd.tex[j])) / fd.teE[j];
            chi2t += r * r;
            ++nt;
        }
        fd.chi2bin = nt > 0 ? chi2t / nt : 0.0;
        fd.par.assign(ff->GetNpar(), 0.0);
        for (int ip = 0; ip < ff->GetNpar(); ++ip) fd.par[ip] = ff->GetParameter(ip);
        delete ff;
    }

    // Draw one fold, including the fitted function and the training and
    // held-out data points.
    auto renderFold = [&](int f, double txtSize) {
        const FoldData& fd = folds[f];
        TH1F* fr = new TH1F(Form("cvfr_%d_%p", f, (void*)gPad), "", 10, xmin, xmax);
        fr->SetStats(0);
        fr->SetMinimum(0);
        fr->SetMaximum(1.30 * ymax);
        fr->SetTitle(Form("fold %d / %d;m [GeV];Background template (events)", f + 1, k));
        fr->GetYaxis()->SetTitleOffset(1.6);
        fr->Draw();
        TF1* ff = (TF1*)templateFunc->Clone(Form("cvdraw_%d_%p", f, (void*)gPad));
        for (size_t ip = 0; ip < fd.par.size(); ++ip) ff->SetParameter((int)ip, fd.par[ip]);
        ff->SetLineColor(kAzure + 2);
        ff->SetLineWidth(2);
        ff->Draw("L same");
        TGraphErrors* gtr = new TGraphErrors();
        for (size_t j = 0; j < fd.trx.size(); ++j) {
            gtr->SetPoint((int)j, fd.trx[j], fd.trY[j]);
            gtr->SetPointError((int)j, 0, fd.trE[j]);
        }
        gtr->SetMarkerStyle(20);
        gtr->SetMarkerSize(0.7);
        gtr->SetMarkerColor(kBlack);
        gtr->SetLineColor(kBlack);
        gtr->Draw("P same");
        TGraphErrors* gte = new TGraphErrors();
        for (size_t j = 0; j < fd.tex.size(); ++j) {
            gte->SetPoint((int)j, fd.tex[j], fd.teY[j]);
            gte->SetPointError((int)j, 0, fd.teE[j]);
        }
        gte->SetMarkerStyle(24);
        gte->SetMarkerSize(1.1);
        gte->SetMarkerColor(kRed + 1);
        gte->SetLineColor(kRed + 1);
        gte->Draw("P same");
        TLatex tx;
        tx.SetNDC();
        tx.SetTextSize(txtSize);
        tx.DrawLatex(0.20, 0.84, Form("held-out #chi^{2}/bin = %.2f", fd.chi2bin));
    };
    // Proxy for the colour legend (train=black filled, held-out=red open, curve=blue).
    auto makeLegend = [&](double x1, double y1, double x2, double y2, double ts, bool header) {
        TGraph* lTr = new TGraph();
        lTr->SetMarkerStyle(20);
        lTr->SetMarkerColor(kBlack);
        lTr->SetMarkerSize(1.2);
        TGraph* lTe = new TGraph();
        lTe->SetMarkerStyle(24);
        lTe->SetMarkerColor(kRed + 1);
        lTe->SetMarkerSize(1.4);
        TLine* lFt = new TLine();
        lFt->SetLineColor(kAzure + 2);
        lFt->SetLineWidth(3);
        TLegend* leg = new TLegend(x1, y1, x2, y2);
        leg->SetBorderSize(1);
        leg->SetFillStyle(1001);
        leg->SetFillColor(kWhite);
        leg->SetTextSize(ts);
        if (header) leg->SetHeader(Form("k-fold CV: %s", fname.c_str()));
        leg->AddEntry(lTr, "Training bins", "p");
        leg->AddEntry(lTe, "Held-out bins", "p");
        leg->AddEntry(lFt, "Fit to training bins", "l");
        leg->Draw();
    };

    // Legend distinguishing training bins, held-out bins and the fitted function.
    int cols = (int)std::ceil(std::sqrt((double)k)), rows = (int)std::ceil((double)k / cols);
    int slots = cols * rows;
    TCanvas* c = new TCanvas("c_cv_allfolds", "CV folds", 360 * cols, 300 * rows);
    c->Divide(cols, rows);
    for (int f = 0; f < k; ++f) {
        c->cd(f + 1);
        gPad->SetLeftMargin(0.17);
        gPad->SetBottomMargin(0.15);
        renderFold(f, 0.062);
    }
    if (k < slots) {
        c->cd(k + 1);
        makeLegend(0.10, 0.45, 0.92, 0.80, 0.075, true);
    } else {
        c->cd(1);
        makeLegend(0.56, 0.68, 0.93, 0.88, 0.040, false);
    }
    TLatex tt;
    tt.SetNDC();
    tt.SetTextSize(0.020);
    tt.SetTextColor(kGray + 3);
    c->cd(0);
    tt.DrawLatex(0.01, 0.985, Form("k-fold cross-validation: %s", fname.c_str()));
    c->SaveAs(outPath.c_str());
    delete c;

    // Save all folds on one canvas, using a spare panel for the legend
    std::string base = outPath, ext = ".pdf";
    std::string::size_type dot = base.rfind('.');
    if (dot != std::string::npos) {
        ext = base.substr(dot);
        base = base.substr(0, dot);
    }
    for (int f = 0; f < k; ++f) {
        TCanvas* cf = new TCanvas(Form("c_cv_fold%d", f), "CV fold", 640, 520);
        cf->SetLeftMargin(0.17);
        cf->SetBottomMargin(0.14);
        cf->SetRightMargin(0.05);
        cf->SetTopMargin(0.08);
        renderFold(f, 0.045);
        makeLegend(0.60, 0.66, 0.94, 0.90, 0.034, true);
        cf->SaveAs(Form("%s_fold%d%s", base.c_str(), f + 1, ext.c_str()));
        delete cf;
    }
}

// Plots the finalists' CV scores and uncertainties against model complexity.
// The 1-SE band identifies candidates with comparable predictive performance;
// the selected model and the in-sample best fit are highlighted separately.
// Finalists' CV +- SE against their number of parameters, the 1-SE band,
// the chosen one (red star) and the in-sample best if different (grey diamond).
inline void DrawCVDecision(const std::vector<FitFunction*>& finalists, const FitFunction* chosen,
                           const std::vector<FitFunction>& allFits, double cvMin, double cvThreshold,
                           const std::string& outPath, const std::string& titleTag) {
    if (finalists.empty()) return;
    const FitFunction* champ = InSampleChampion(allFits);
    int nmin = 1000000, nmax = -1000000;
    double ylo = 1e9, yhi = -1e9;
    for (auto* ff : finalists) {
        if (!ff->cv_valid) continue;
        nmin = std::min(nmin, ff->npar);
        nmax = std::max(nmax, ff->npar);
        ylo = std::min(ylo, ff->cv_score - ff->cv_sem);
        yhi = std::max(yhi, ff->cv_score + ff->cv_sem);
    }
    if (champ && champ->cv_valid) {
        nmin = std::min(nmin, champ->npar);
        nmax = std::max(nmax, champ->npar);
        ylo = std::min(ylo, champ->cv_score);
        yhi = std::max(yhi, champ->cv_score);
    }
    if (nmin > nmax) {
        nmin = 2;
        nmax = 6;
    }
    if (ylo > yhi) {
        ylo = 0;
        yhi = 2;
    }
    const double fxlo = nmin - 1.0, fxhi = nmax + 1.0;
    double ypad = 0.18 * (yhi - ylo) + 1e-6;
    double ymin = ylo - ypad, ymax = yhi + 3.0 * ypad;
    TCanvas* c = new TCanvas("c_cv_decision", "CV decision", 760, 580);
    c->SetLeftMargin(0.15);
    c->SetBottomMargin(0.14);
    c->SetRightMargin(0.05);
    c->SetTopMargin(0.06);
    TH1F* fr = new TH1F("cvdecfr", "", 10, fxlo, fxhi);
    fr->SetStats(0);
    fr->SetMinimum(ymin);
    fr->SetMaximum(ymax);
    fr->SetTitle("");
    fr->GetXaxis()->SetTitle("Number of parameters (model DoF)");
    fr->GetYaxis()->SetTitle("CV #chi^{2}/bin");
    fr->GetXaxis()->SetTitleSize(0.050);
    fr->GetXaxis()->SetTitleOffset(1.18);
    fr->GetXaxis()->SetLabelSize(0.046);
    fr->GetYaxis()->SetTitleSize(0.054);
    fr->GetYaxis()->SetTitleOffset(1.10);
    fr->GetYaxis()->SetLabelSize(0.046);
    fr->GetXaxis()->SetNdivisions(nmax - nmin + 2, kFALSE);
    fr->Draw();

    for (int n = nmin; n <= nmax; ++n) {
        TLine* gl = new TLine(n, ymin, n, ymax);
        gl->SetLineColor(kGray);
        gl->SetLineStyle(3);
        gl->SetLineWidth(1);
        gl->Draw();
    }
    TBox* band = new TBox(fxlo, cvMin, fxhi, cvThreshold);
    band->SetFillColorAlpha(kGreen + 2, 0.18);
    band->SetLineColor(kGreen + 2);
    band->Draw("same");
    TLine* lthr = new TLine(fxlo, cvThreshold, fxhi, cvThreshold);
    lthr->SetLineColor(kGreen + 2);
    lthr->SetLineStyle(2);
    lthr->SetLineWidth(2);
    lthr->Draw();

    // finalists with the same npar are spread sideways so the error bars do not overlap
    const double dxoff = 0.24;
    std::vector<double> foff(finalists.size(), 0.0);
    for (size_t i = 0; i < finalists.size(); ++i) {
        if (!finalists[i]->cv_valid) continue;
        int n = 0, slot = 0;
        for (size_t j = 0; j < finalists.size(); ++j) {
            if (!finalists[j]->cv_valid || finalists[j]->npar != finalists[i]->npar) continue;
            if (j < i) ++slot;
            ++n;
        }
        foff[i] = (slot - (n - 1) / 2.0) * dxoff;
    }
    auto offByName = [&](const std::string& nm) -> double {
        for (size_t i = 0; i < finalists.size(); ++i)
            if (finalists[i]->name == nm) return foff[i];
        return 0.0;
    };
    const double labOff = 0.04 * (yhi - ylo);
    for (size_t i = 0; i < finalists.size(); ++i) {
        auto* ff = finalists[i];
        if (!ff->cv_valid) continue;
        double xx = ff->npar + foff[i];
        TGraphErrors* g = new TGraphErrors();
        g->SetPoint(0, xx, ff->cv_score);
        g->SetPointError(0, 0, ff->cv_sem);
        g->SetMarkerStyle(20);
        g->SetMarkerSize(1.5);
        g->SetMarkerColor(kBlack);
        g->SetLineColor(kBlack);
        g->SetLineWidth(2);
        g->Draw("P same");
        TLatex t;
        t.SetTextSize(0.032);
        if (foff[i] < -0.01) {
            t.SetTextAlign(32);
            t.DrawLatex(xx - 0.10, ff->cv_score, ff->name.c_str());
        } else if (foff[i] > 0.01) {
            t.SetTextAlign(12);
            t.DrawLatex(xx + 0.10, ff->cv_score, ff->name.c_str());
        } else {
            t.SetTextAlign(21);
            t.DrawLatex(xx, ff->cv_score + ff->cv_sem + labOff, ff->name.c_str());
        }
    }
    if (champ && chosen && champ != chosen && champ->cv_valid) {
        TGraph* gc = new TGraph();
        gc->SetPoint(0, champ->npar + offByName(champ->name), champ->cv_score);
        gc->SetMarkerStyle(34);
        gc->SetMarkerSize(2.0);
        gc->SetMarkerColor(kGray + 2);
        gc->Draw("P same");
    }
    if (chosen && chosen->cv_valid) {
        TGraph* gs = new TGraph();
        gs->SetPoint(0, chosen->npar + offByName(chosen->name), chosen->cv_score);
        gs->SetMarkerStyle(29);
        gs->SetMarkerSize(2.9);
        gs->SetMarkerColor(kRed + 1);
        gs->Draw("P same");
    }
    // Legend (top-right): explains the markers, the 1-SE band and its threshold.
    TGraph* legStar = new TGraph();
    legStar->SetMarkerStyle(29);
    legStar->SetMarkerColor(kRed + 1);
    legStar->SetMarkerSize(2.6);
    TGraph* legDot = new TGraph();
    legDot->SetMarkerStyle(20);
    legDot->SetMarkerColor(kBlack);
    legDot->SetMarkerSize(1.5);
    TGraph* legDia = new TGraph();
    legDia->SetMarkerStyle(34);
    legDia->SetMarkerColor(kGray + 2);
    legDia->SetMarkerSize(2.0);
    TLegend* leg = new TLegend(0.60, 0.72, 0.96, 0.955);
    leg->SetBorderSize(1);
    leg->SetFillStyle(1001);
    leg->SetFillColor(kWhite);
    leg->SetTextSize(0.028);
    leg->AddEntry(legStar, "Optimal model (simplest in band)", "p");
    leg->AddEntry(legDot, "Candidate fits (CV #chi^{2} #pm SE)", "p");
    if (champ && chosen && champ != chosen && champ->cv_valid)
        leg->AddEntry(legDia, "In-sample best (demoted)", "p");
    leg->AddEntry(band, "best CV + 1 SE band", "f");
    leg->AddEntry(lthr, "1-SE threshold", "l");
    leg->Draw();
    (void)titleTag;
    c->SaveAs(outPath.c_str());
    delete c;
}

// ============================================================================
// Selection
// ============================================================================
// Returns a copy of the chosen FitFunction with its own TF1; the caller owns
// both. The spurious signal of the finalists is filled in on the way.
FitFunction* SelectOptimalBackgroundFit(TH1F* h_sideband, TH1F* h_full_data, TH1F* h_mc_template,
                                        TH1F* h_multijets_scaled, double xmin, double xmax, int binSignalLow,
                                        int binSignalHigh, const std::string& outputPath,
                                        const std::string& ptLabel = "") {

    std::cout << "\n========================================" << std::endl;
    std::cout << "F-TEST BACKGROUND FUNCTION SELECTION" << std::endl;
    std::cout << "========================================\n" << std::endl;

    // --- 1. all candidates on the data sidebands --------------------------

    std::vector<FTestComparison> ftest_results;
    std::vector<FitFunction> allFits;

    // Create polynomial, exponential and polynomial-exponential candidates
    // for orders 2 through 5.
    std::cout << "Creating candidate background functions (orders 2-5)...\n" << std::endl;
    for (int order = 2; order <= 5; ++order) {
        std::string name = "polyexp_" + std::to_string(order);
        TF1* func = CreatePolyExpFit(order, xmin, xmax, name);
        allFits.push_back(FitFunction(func, name, "PolyExp", order));

        name = "poly_" + std::to_string(order);
        func = CreatePolyFit(order, xmin, xmax, name);
        allFits.push_back(FitFunction(func, name, "Poly", order));

        name = "exp_" + std::to_string(order);
        func = CreateExpFit(order, xmin, xmax, name);
        allFits.push_back(FitFunction(func, name, "Exp", order));
    }

    std::cout << "Step 1: Fitting all candidate functions to sideband regions...\n" << std::endl;
    std::cout << std::setw(20) << "Function" << std::setw(12) << "Type" << std::setw(8) << "Order"
              << std::setw(10) << "NPar" << std::setw(12) << "Chi^2" << std::setw(8) << "NDF" << std::setw(12)
              << "Chi^2/NDF" << std::endl;
    std::cout << std::string(82, '-') << std::endl;


    // Fit all functions
    for (auto& fit : allFits) {
        h_sideband->Fit(fit.func, "RQ0");
        fit.chi2 = fit.func->GetChisquare();
        fit.ndf = fit.func->GetNDF();
        fit.npar = fit.func->GetNpar();

        // ========== VALIDATION SECTION ==========
        bool fit_valid = true;
        if (!std::isfinite(fit.chi2) || fit.chi2 > 1e6 || fit.ndf <= 0) {
            fit_valid = false;
        }

        // Validate all parameters are finite
        for (int i = 0; i < fit.npar; ++i) {
            if (!std::isfinite(fit.func->GetParameter(i))) {
                fit_valid = false;
                break;
            }
        }

        //fit.chi2_ndf = fit_valid ? ((fit.ndf > 0) ? fit.chi2 / fit.ndf : 999) : 999;
        if (!fit_valid) {
            std::cout << std::setw(20) << fit.name << std::setw(12) << fit.type << std::setw(8) << fit.order
                      << std::setw(10) << fit.npar << "      **FIT FAILED**" << std::endl;
            fit.chi2_ndf = 999;
            continue;
        }
        // ========== END VALIDATION SECTION ==========

        fit.chi2_ndf = (fit.ndf > 0) ? fit.chi2 / fit.ndf : 999;

        std::cout << std::setw(20) << fit.name << std::setw(12) << fit.type << std::setw(8) << fit.order
                  << std::setw(10) << fit.npar << std::setw(12) << std::fixed << std::setprecision(2)
                  << fit.chi2 << std::setw(8) << fit.ndf << std::setw(12) << std::setprecision(4)
                  << fit.chi2_ndf << std::endl;
    }

    // --- 2. F-test within each family ------------------------------------
    std::cout << std::endl;

    std::cout << "\nStep 2: F-tests between consecutive orders of the same type...\n" << std::endl;

    std::vector<FitFunction*> selectedFits;
    std::vector<std::string> types = {"PolyExp", "Poly", "Exp"};

    for (const auto& type : types) {
        std::cout << "\n" << std::string(60, '=') << std::endl;
        std::cout << "Testing " << type << " Functions (Orders 2-5)" << std::endl;
        std::cout << std::string(60, '=') << std::endl;

        std::vector<FitFunction*> typeFits;
        for (auto& fit : allFits) {
            if (fit.type == type) typeFits.push_back(&fit);
        }
        std::sort(typeFits.begin(), typeFits.end(),
                  [](FitFunction* a, FitFunction* b) { return a->order < b->order; });

        FitFunction* selectedFit = typeFits[0];
        std::cout << "\nStarting with simplest model: " << selectedFit->name
                  << " (Chi^2/ndf = " << selectedFit->chi2_ndf << ")\n"
                  << std::endl;

        for (size_t i = 1; i < typeFits.size(); ++i) {
            FitFunction* simpler = typeFits[i - 1];
            FitFunction* complex = typeFits[i];

            // Skip if either fit failed
            if (simpler->chi2_ndf >= 999 || complex->chi2_ndf >= 999) {
                std::cout << "Skipping F-test: " << simpler->name << " vs " << complex->name
                          << " (one or both fits failed)" << std::endl;
                break; // Stop testing higher orders for this function type
            }

            FTestComparison comparison;
            comparison.simpler_model = simpler->name;
            comparison.complex_model = complex->name;
            comparison.simpler_npar = simpler->npar;
            comparison.complex_npar = complex->npar;
            comparison.simpler_chi2 = simpler->chi2;
            comparison.complex_chi2 = complex->chi2;
            comparison.simpler_ndf = simpler->ndf;
            comparison.complex_ndf = complex->ndf;

            int n_bins = simpler->ndf + simpler->npar;
            double F = CalculateFTest(simpler->chi2, complex->chi2, simpler->npar, complex->npar, n_bins);

            comparison.df1 = complex->npar - simpler->npar;
            comparison.df2 = complex->ndf;
            comparison.F_statistic = F;
            comparison.F_critical = GetFCritical(comparison.df1, comparison.df2);
            comparison.F_critical = GetFCritical(comparison.df1, comparison.df2, 0.05);

            comparison.p_value = 1.0 - TMath::FDistI(F, comparison.df1, comparison.df2);
            comparison.significant = (comparison.p_value < 0.05);

            std::cout << "F-test: " << simpler->name << " vs " << complex->name << std::endl;
            std::cout << "  H0: Simpler model (" << simpler->name << ") is adequate" << std::endl;
            std::cout << "  H1: Complex model (" << complex->name << ") provides significant improvement"
                      << std::endl;
            std::cout << "  ---" << std::endl;
            std::cout << "  DeltaChi^2 = " << (simpler->chi2 - complex->chi2) << std::endl;
            std::cout << "  Deltadf = " << (simpler->ndf - complex->ndf)
                      << " (additional parameters: " << comparison.df1 << ")" << std::endl;
            std::cout << "  F-statistic = " << std::setprecision(4) << F << std::endl;
            std::cout << "  F-critical (Alpha=0.05, df1=" << comparison.df1 << ", df2=" << comparison.df2
                      << ") = " << std::setprecision(4) << comparison.F_critical << std::endl;
            std::cout << "  p-value = " << std::setprecision(6) << comparison.p_value << std::endl;

            if (comparison.significant) {
                std::cout << "  *** SIGNIFICANT (p < 0.05): Reject H0 ***" << std::endl;
                std::cout << "  -> Complex model (" << complex->name
                          << ") provides statistically significant improvement" << std::endl;
                comparison.decision = "Reject H0 - Use complex model";
                selectedFit = complex;
            } else {
                std::cout << "  NOT SIGNIFICANT (p ≥ 0.05): Accept H0" << std::endl;
                std::cout << "  -> Simpler model (" << simpler->name << ") is adequate" << std::endl;
                comparison.decision = "Accept H0 - Use simpler model";
            }
            std::cout << std::endl;

            ftest_results.push_back(comparison);

            if (!comparison.significant) {
                break;
            }
        }

        std::cout << ">>> SELECTED for " << type << ": " << selectedFit->name
                  << " (Chi^2/ndf = " << selectedFit->chi2_ndf << ") <<<\n"
                  << std::endl;
        selectedFits.push_back(selectedFit);
    }

    // --- 3. finalists on the MC template, CV of every candidate ---------
    std::cout << "\n" << std::string(60, '=') << std::endl;

    std::cout << "\n" << std::string(60, '=') << std::endl;
    std::cout << "Step 3: Dual Background Validation" << std::endl;
    std::cout << std::string(60, '=') << std::endl;
    std::cout << "Testing selected functions against both data sidebands AND MC background...\n" << std::endl;

    // Fit each selected function to the MC background template (h_multijets_scaled)
    std::cout << "3a: Fitting to MC Background Template ("
              << SampleConfig::DisplayNames(CurrentChannel())[SampleConfig::kMultijets] << "):\n"
              << std::endl;
    std::cout << std::setw(20) << "Function" << std::setw(15) << "Data Chi^2/NDF" << std::setw(15)
              << "MC Chi^2/NDF" << std::setw(15) << "Combined Score" << std::endl;
    std::cout << std::string(65, '-') << std::endl;

    // Store MC Chi^2/ndf values separately
    std::vector<double> mc_chi2_ndf_values;

    for (auto fit : selectedFits) {
        // Already have fit to data sidebands (fit->chi2, fit->chi2_ndf)
        double data_chi2_ndf = fit->chi2_ndf;

        // Now fit the same function to MC background template
        TF1* mc_fit = (TF1*)fit->func->Clone((fit->name + "_mc").c_str());

        // Reset parameters to initial values for MC fit
        for (int i = 0; i < mc_fit->GetNpar(); ++i) {
            mc_fit->SetParameter(i, fit->func->GetParameter(i));
        }

        // Fit to MC background (using same sideband regions)
        h_multijets_scaled->Fit(mc_fit, "RQ0");

        double mc_chi2 = mc_fit->GetChisquare();
        int mc_ndf = mc_fit->GetNDF();
        double mc_chi2_ndf = (mc_ndf > 0) ? mc_chi2 / mc_ndf : 999;

        // Store in separate vector instead of spurious_signal
        mc_chi2_ndf_values.push_back(mc_chi2_ndf);

        // Combined score: average of the two Chi^2/ndf values
        double combined_score = (data_chi2_ndf + mc_chi2_ndf) / 2.0;

        std::cout << std::setw(20) << fit->name << std::setw(15) << std::setprecision(4) << data_chi2_ndf
                  << std::setw(15) << std::setprecision(4) << mc_chi2_ndf << std::setw(15)
                  << std::setprecision(4) << combined_score << std::endl;

        delete mc_fit;
    }

    std::cout << "\n3b) Selection Logic:\n" << std::endl;
    std::cout << "Goal: Pick the function that GENERALISES best to the background shape" << std::endl;
    std::cout << "Method: k-fold cross-validation on the MC template (sideband region) + 1-standard-error"
              << std::endl;
    std::cout << "        parsimony rule (simplest model within 1 SE of the best CV score)" << std::endl;
    std::cout << std::string(60, '-') << std::endl;

    // k-fold cross-validation on the MC background template
    // CV scores each finalist by its OUT-OF-SAMPLE Chi^2 on held-out bins; the
    // 1-standard-error parsimony rule then keeps the simplest function whose CV score
    // is decisively better out-of-sample within 1 SE of the best
    // CV region on the MC background template - set via env var with EASYBJETS_CV_REGION:
    //   sidebands (default) : restrict CV to the sideband region (signal region window
    //                         zeroed)-apples-to-apples with the F-test, which
    //                         fits the data sidebands.
    //   full                : use the whole [xmin,xmax] MC template, incl. under
    //                         the peak (MC is signal-free there).
    // We always clone the (signal-free, high-stat) MC template; in "sidebands"
    // mode we zero the signal-region bins so KFoldCV's non-zero-error filter
    // keeps only the sideband bins.
    const int kCVFolds = 5;
    bool cvSidebandsOnly = true;
    if (const char* e = std::getenv("EASYBJETS_CV_REGION")) {
        std::string s(e);
        for (auto& c : s) c = (char)std::tolower((unsigned char)c);
        if (s == "full" || s == "fullrange" || s == "all") cvSidebandsOnly = false;
    }
    TH1F* h_cv = (TH1F*)h_multijets_scaled->Clone("h_mc_cv");
    if (cvSidebandsOnly)
        for (int b = binSignalLow; b <= binSignalHigh; ++b) {
            h_cv->SetBinContent(b, 0.0);
            h_cv->SetBinError(b, 0.0);
        }
    std::cout << "  [CV region] " << (cvSidebandsOnly ? "sidebands (signal window excluded)" : "full range")
              << "   (set EASYBJETS_CV_REGION=full|sidebands)" << std::endl;
    // Compute CV for ALL candidates (stored in the struct)-used both for the
    // decision (finalists) and for the U-curve diagnostic (every order).
    for (auto& f : allFits) {
        if (f.chi2_ndf >= 999) continue;
        CVResult r = KFoldCV(f.func, h_cv, xmin, xmax, kCVFolds);
        f.cv_score = r.score;
        f.cv_sem = r.sem;
        f.cv_valid = r.valid;
    }

    std::vector<CVResult> cv(selectedFits.size());
    std::cout << "\n3c) k-fold cross-validation on MC template, sideband region (k=" << kCVFolds << "):\n"
              << std::endl;
    std::cout << std::setw(20) << "Function" << std::setw(8) << "NPar" << std::setw(18)
              << "in-sample Chi^2/ndf" << std::setw(16) << "CV Chi^2/bin" << std::setw(12) << "CV +/-SE"
              << std::endl;
    std::cout << std::string(74, '-') << std::endl;
    bool anyCV = false;
    for (size_t i = 0; i < selectedFits.size(); ++i) {
        cv[i].score = selectedFits[i]->cv_score;
        cv[i].sem = selectedFits[i]->cv_sem;
        cv[i].valid = selectedFits[i]->cv_valid;
        anyCV = anyCV || cv[i].valid;
        std::cout << std::setw(20) << selectedFits[i]->name << std::setw(8) << selectedFits[i]->npar
                  << std::setw(18) << std::setprecision(4) << selectedFits[i]->chi2_ndf;
        if (cv[i].valid)
            std::cout << std::setw(16) << std::setprecision(4) << cv[i].score << std::setw(12)
                      << std::setprecision(3) << cv[i].sem << std::endl;
        else
            std::cout << std::setw(16) << "n/a" << std::setw(12) << "(few bins)" << std::endl;
    }

    FitFunction* bestFit = selectedFits[0];
    int best_idx = 0;
    double cvMin = 0.0, cvThreshold = 0.0; // for the U-curve diagnostic

    if (anyCV) {
        // Best (lowest) CV score among valid candidates.
        int cvBest = -1;
        for (size_t i = 0; i < cv.size(); ++i)
            if (cv[i].valid && (cvBest < 0 || cv[i].score < cv[cvBest].score)) cvBest = (int)i;
        cvMin = cv[cvBest].score;
        cvThreshold = cv[cvBest].score + cv[cvBest].sem; // 1-SE band
        // So among candidates within 1 SE of the best CV, take the
        // fewest parameters; for EQUAL parameter counts the CV scores are within
        // noise of each other, so break the tie by the better combined data+MC
        // Chi^2/ndf-NOT by the (noisy) sub-SE CV difference.
        auto combinedFit = [&](int i) {
            return (selectedFits[i]->chi2_ndf + mc_chi2_ndf_values[i]) / 2.0;
        };
        int chosen = cvBest;
        for (size_t i = 0; i < cv.size(); ++i) {
            if (!cv[i].valid || cv[i].score > cvThreshold) continue;
            if (selectedFits[i]->npar < selectedFits[chosen]->npar
                || (selectedFits[i]->npar == selectedFits[chosen]->npar
                    && combinedFit((int)i) < combinedFit(chosen)))
                chosen = (int)i;
        }
        best_idx = chosen;
        bestFit = selectedFits[chosen];
        std::cout << "\n  Best CV: " << selectedFits[cvBest]->name << " (" << std::setprecision(4)
                  << cv[cvBest].score << ")"
                  << "  ->  1-SE threshold = " << std::setprecision(4) << cvThreshold << "\n"
                  << "  Pick (within 1 SE; fewest params, then best combined Chi^2/ndf): " << bestFit->name
                  << std::endl;
    } else {
        // Fallback: too few bins for CV -> keep the old lowest-combined-score rule.
        std::cout << "\n  [CV unavailable-falling back to lowest combined Chi^2/ndf]" << std::endl;
        double best_combined_score = (bestFit->chi2_ndf + mc_chi2_ndf_values[0]) / 2.0;
        for (size_t i = 1; i < selectedFits.size(); ++i) {
            double combined_score = (selectedFits[i]->chi2_ndf + mc_chi2_ndf_values[i]) / 2.0;
            if (combined_score < best_combined_score) {
                bestFit = selectedFits[i];
                best_combined_score = combined_score;
                best_idx = (int)i;
            }
        }
    }
    // Kept for the printout / fit_summary below (display only-not the decision).
    double best_combined_score = (bestFit->chi2_ndf + mc_chi2_ndf_values[best_idx]) / 2.0;

    // --- 4. decision -----------------------------------------------------
    // EASYBJETS_BKG_SELECT = cv   (default) : k-fold CV + 1-SE parsimony (above)
    //                        = dof            : smallest DoF (fewest parameters)
    //                                           among the F-test finalists, ties
    //                                           broken by best combined data+MC Chi^2/ndf
    // CV is always computed & shown (table + plots); only the DECISION changes.
    std::string bkgSelMode = "cv";
    if (const char* e = std::getenv("EASYBJETS_BKG_SELECT")) {
        std::string s(e);
        for (auto& c : s) c = (char)std::tolower((unsigned char)c);
        if (s == "dof" || s == "mindof" || s == "parsimony" || s == "ftest") bkgSelMode = "dof";
    }
    if (bkgSelMode == "dof") {
        auto comb = [&](int i) {
            return (selectedFits[i]->chi2_ndf + mc_chi2_ndf_values[i]) / 2.0;
        };
        int chosen = 0;
        for (size_t i = 1; i < selectedFits.size(); ++i)
            if (selectedFits[i]->npar < selectedFits[chosen]->npar
                || (selectedFits[i]->npar == selectedFits[chosen]->npar && comb((int)i) < comb(chosen)))
                chosen = (int)i;
        best_idx = chosen;
        bestFit = selectedFits[chosen];
        best_combined_score = comb(chosen);
        std::cout << "\n  [SELECT=dof] Smallest-DoF F-test finalist "
                     "(ties -> best combined chi2/ndf): "
                  << bestFit->name << "  (npar=" << bestFit->npar << ")" << std::endl;
    } else {
        std::cout << "\n  [SELECT=cv] using k-fold CV + 1-SE pick: " << bestFit->name
                  << "   (set EASYBJETS_BKG_SELECT=dof to switch to smallest-DoF)" << std::endl;
    }

    // Group all CV outputs (table + plots) in a cross_validation/subfolder.
    // The table lists EVERY candidate (all orders), so the in-sample
    // vs out-of-sample (overfitting) trend is visible alongside the U-curve plot.
    {
        std::string cvDir = outputPath + "cross_validation/";
        gSystem->Exec(("mkdir -p " + cvDir).c_str());

        std::ofstream cv_file((cvDir + "cv_results.txt").c_str());
        if (cv_file.is_open()) {
            cv_file << "k-FOLD CROSS-VALIDATION (k=" << kCVFolds
                    << ") ON MC BACKGROUND TEMPLATE (SIDEBAND REGION)\n";
            cv_file << "Decision: 1-SE rule-simplest model within 1 SE of the best CV;\n";
            cv_file << "          equal-parameter ties broken by the better combined data+MC chi2/ndf.\n";
            cv_file << "  * = family finalist (survived the within-family F-test)\n";
            cv_file << "  <<< = selected (optimal) function\n";
            cv_file << "==========================================================================\n\n";
            cv_file << std::setw(20) << "Function" << std::setw(8) << "NPar" << std::setw(18)
                    << "in-sample chi2/ndf" << std::setw(16) << "CV chi2/bin" << std::setw(12) << "CV SE"
                    << std::setw(8) << "final" << std::setw(8) << "opt" << "\n";
            cv_file << std::string(90, '-') << "\n";
            for (const auto& f : allFits) {
                bool isFinalist = false;
                for (auto* sf : selectedFits)
                    if (sf == &f) isFinalist = true;
                bool isChosen = (&f == bestFit);
                cv_file << std::setw(20) << f.name << std::setw(8) << f.npar;
                if (f.chi2_ndf < 999)
                    cv_file << std::setw(18) << std::setprecision(4) << f.chi2_ndf;
                else
                    cv_file << std::setw(18) << "FAILED";
                if (f.cv_valid)
                    cv_file << std::setw(16) << std::setprecision(4) << f.cv_score << std::setw(12)
                            << std::setprecision(3) << f.cv_sem;
                else
                    cv_file << std::setw(16) << "n/a" << std::setw(12) << "-";
                cv_file << std::setw(8) << (isFinalist ? "*" : "") << std::setw(8) << (isChosen ? "<<<" : "")
                        << "\n";
            }
            cv_file << "\n1-SE band: [" << std::setprecision(4) << cvMin << ", " << cvThreshold << "]\n";
            cv_file << "Selected: " << bestFit->name << "  (out-of-sample CV, parsimony-first)\n";
            cv_file.close();
        }

        // CV diagnostic plots (always on).
        std::string titleTag = outputPath;
        DrawCVScoreVsComplexity(allFits, bestFit, cvMin, cvThreshold, cvDir, titleTag);
        if (bestFit->cv_valid)
            DrawCVFoldIllustration(bestFit->func, h_cv, xmin, xmax, kCVFolds, cvDir,
                                   "cv_folds_" + bestFit->name + ".pdf", titleTag);

        // Detailed CV diagnostics for the finalists + chosen (always on):
        //   - cv_folds.txt          per-fold tables + 1-SE decision + overfitting check
        //   - cv_decision.pdf       finalists' CV+/-SE vs DoF, 1-SE band, pick, demotion
        //   - cv_folds_all_<f>.pdf  all k folds (train vs held-out) for the chosen
        // All use h_cv = MC template restricted to the sideband region.
        WriteCVFoldsTxt(cvDir + "cv_folds.txt", selectedFits, bestFit, allFits, h_cv, xmin, xmax, kCVFolds,
                        cvMin, cvThreshold);
        DrawCVDecision(selectedFits, bestFit, allFits, cvMin, cvThreshold, cvDir + "cv_decision.pdf",
                       titleTag);
        // cv_folds_all_<name>.pdf for EVERY finalist (chosen + the other candidates),
        // not just the optimal one, so all per-type winners can be compared.
        for (auto* ff : selectedFits)
            if (ff && ff->cv_valid)
                DrawCVFoldsAll(ff->func, h_cv, xmin, xmax, kCVFolds,
                               cvDir + "cv_folds_all_" + ff->name + ".pdf", titleTag);
    }

    // Calculate spurious signal for all selected fits - this kept in case we want to use it in the future
    for (size_t i = 0; i < selectedFits.size(); ++i) {
        selectedFits[i]->spurious_signal = CalculateSpuriousSignal(
            h_full_data, h_mc_template, selectedFits[i]->func, binSignalLow, binSignalHigh);
    }

    // Get the MC Chi^2/ndf for the best fit
    double mc_chi2_ndf_best = mc_chi2_ndf_values[best_idx];

    std::cout << "\n" << std::string(60, '=') << std::endl;
    std::cout << "OPTIMAL FUNCTION SELECTED" << std::endl;
    std::cout << std::string(60, '=') << std::endl;
    std::cout << "Function:              " << bestFit->name << std::endl;
    std::cout << "Type:                  " << bestFit->type << std::endl;
    std::cout << "Order:                 " << bestFit->order << std::endl;
    std::cout << "N Parameters:          " << bestFit->npar << std::endl;
    std::cout << "---" << std::endl;
    std::cout << "Data Sidebands Chi^2:     " << std::setprecision(2) << bestFit->chi2 << std::endl;
    std::cout << "Data Sidebands NDF:    " << bestFit->ndf << std::endl;
    std::cout << "Data Sidebands Chi^2/NDF: " << std::setprecision(4) << bestFit->chi2_ndf << std::endl;
    std::cout << "---" << std::endl;
    std::cout << "MC Background Chi^2/NDF:  " << std::setprecision(4) << mc_chi2_ndf_best << std::endl;
    std::cout << "Combined Score:        " << std::setprecision(4) << best_combined_score << std::endl;
    if (ShowSpuriousSignal()) {
        std::cout << "---" << std::endl;
        std::cout << "Spurious Signal:       " << std::setprecision(2) << bestFit->spurious_signal
                  << " events (reference)" << std::endl;
    }
    std::cout << std::string(60, '=') << "\n" << std::endl;

    std::cout << "\nInterpretation:" << std::endl;
    std::cout << "  -> " << bestFit->name << " fits data sidebands with Chi^2/ndf = " << std::setprecision(3)
              << bestFit->chi2_ndf << std::endl;
    std::cout << "  -> " << bestFit->name << " fits MC background with Chi^2/ndf = " << std::setprecision(3)
              << mc_chi2_ndf_best << std::endl;
    std::cout << "  -> This function best describes BOTH the data and MC background shape" << std::endl;

    // Create individual analysis for each selected function
    std::cout << "\n" << std::string(60, '=') << std::endl;
    std::cout << "Step 3b: Creating Individual Fit Analyses" << std::endl;
    std::cout << std::string(60, '=') << std::endl;


    // Create individual analyses
    for (auto fit : selectedFits) {
        CreateIndividualFitAnalysis(fit, h_full_data, h_mc_template, h_multijets_scaled, outputPath,
                                    /*isOptimalFit=*/(fit == bestFit), ptLabel);
    }

    // Save F-test results
    std::ofstream ftest_file((outputPath + "ftest_results.txt").c_str());
    if (ftest_file.is_open()) {
        ftest_file << "F-TEST BACKGROUND FUNCTION SELECTION RESULTS" << std::endl;
        ftest_file << "============================================\n" << std::endl;

        ftest_file << "ALL CANDIDATE FUNCTIONS:\n" << std::endl;
        ftest_file << std::setw(20) << "Function" << std::setw(12) << "Type" << std::setw(8) << "Order"
                   << std::setw(10) << "NPar" << std::setw(12) << "Chi^2" << std::setw(8) << "NDF"
                   << std::setw(12) << "Chi^2/NDF" << std::endl;
        ftest_file << std::string(82, '-') << std::endl;

        for (const auto& fit : allFits) {
            ftest_file << std::setw(20) << fit.name << std::setw(12) << fit.type << std::setw(8) << fit.order
                       << std::setw(10) << fit.npar;

            // Check if fit failed
            if (fit.chi2_ndf >= 999 || !std::isfinite(fit.chi2)) {
                ftest_file << "      **FIT FAILED**" << std::endl;
            } else {
                ftest_file << std::setw(12) << std::fixed << std::setprecision(2) << fit.chi2 << std::setw(8)
                           << fit.ndf << std::setw(12) << std::setprecision(4) << fit.chi2_ndf << std::endl;
            }
        }


        ftest_file << "\n\nF-TEST COMPARISONS:\n" << std::endl;
        ftest_file << std::string(80, '=') << std::endl;

        for (size_t i = 0; i < ftest_results.size(); ++i) {
            const auto& comp = ftest_results[i];
            ftest_file << "\nComparison " << (i + 1) << ": " << comp.simpler_model << " vs "
                       << comp.complex_model << std::endl;
            ftest_file << std::string(80, '-') << std::endl;
            ftest_file << "Simpler Model: " << comp.simpler_model << std::endl;
            ftest_file << "  Chi^2 = " << std::setprecision(2) << comp.simpler_chi2
                       << ", NDF = " << comp.simpler_ndf << ", NPar = " << comp.simpler_npar << std::endl;
            ftest_file << "Complex Model: " << comp.complex_model << std::endl;
            ftest_file << "  Chi^2 = " << std::setprecision(2) << comp.complex_chi2
                       << ", NDF = " << comp.complex_ndf << ", NPar = " << comp.complex_npar << std::endl;
            ftest_file << "\nF-Test Results:" << std::endl;
            ftest_file << "  F-statistic = " << std::setprecision(4) << comp.F_statistic << std::endl;
            ftest_file << "  F-critical (Alpha=0.05) = " << std::setprecision(4) << comp.F_critical
                       << std::endl;
            ftest_file << "  df1 = " << comp.df1 << ", df2 = " << comp.df2 << std::endl;
            ftest_file << "  p-value = " << std::setprecision(6) << comp.p_value << std::endl;
            ftest_file << "  Significant (p < 0.05): " << (comp.significant ? "YES" : "NO") << std::endl;
            ftest_file << "  Decision: " << comp.decision << std::endl;
        }
        ftest_file << "\n\nFAMILY FINALISTS (MC-template context-reference; the final choice\n"
                   << "is made by k-fold cross-validation, see the next section):\n"
                   << std::endl;
        ftest_file << std::string(80, '=') << std::endl;
        ftest_file << std::setw(15) << "Type" << std::setw(20) << "Function" << std::setw(12)
                   << "Data Chi^2/ndf" << std::setw(12) << "MC Chi^2/ndf" << std::setw(15) << "Combined"
                   << std::setw(15) << "Spurious" << std::endl;
        ftest_file << std::string(89, '-') << std::endl;

        for (const auto& fit : selectedFits) {
            // Find the index of this fit
            size_t idx = std::find(selectedFits.begin(), selectedFits.end(), fit) - selectedFits.begin();
            double mc_chi2_ndf = mc_chi2_ndf_values[idx];
            double combined = (fit->chi2_ndf + mc_chi2_ndf) / 2.0;

            ftest_file << std::setw(15) << fit->type << std::setw(20) << fit->name << std::setw(12)
                       << std::setprecision(4) << fit->chi2_ndf << std::setw(12) << std::setprecision(4)
                       << mc_chi2_ndf << std::setw(15) << std::setprecision(4) << combined;

            // Mark the selected one
            if (fit == bestFit) {
                ftest_file << " *SELECTED*";
            }
            ftest_file << std::endl;
        }

        // ---- k-fold cross-validation: the actual selection criterion ----
        ftest_file << "\n\nk-FOLD CROSS-VALIDATION (k=5) ON MC BACKGROUND TEMPLATE (SIDEBAND REGION):\n"
                   << std::endl;
        ftest_file << "Out-of-sample test that penalises overfitting. This, with the 1-SE\n";
        ftest_file << "parsimony rule, is the criterion that selects the optimal function.\n";
        ftest_file << "An overfit model has a low in-sample chi2/ndf but a WORSE CV chi2/bin.\n";
        ftest_file << "  * = family finalist (passed the within-family F-test)   <<< = selected\n"
                   << std::endl;
        ftest_file << std::setw(20) << "Function" << std::setw(8) << "NPar" << std::setw(20)
                   << "in-sample chi2/ndf" << std::setw(16) << "CV chi2/bin" << std::setw(12) << "CV SE"
                   << std::setw(8) << "final" << std::setw(8) << "opt" << std::endl;
        ftest_file << std::string(92, '-') << std::endl;
        for (const auto& f : allFits) {
            bool isFinalist = false;
            for (auto* sf : selectedFits)
                if (sf == &f) isFinalist = true;
            ftest_file << std::setw(20) << f.name << std::setw(8) << f.npar;
            if (f.chi2_ndf < 999)
                ftest_file << std::setw(20) << std::setprecision(4) << f.chi2_ndf;
            else
                ftest_file << std::setw(20) << "FAILED";
            if (f.cv_valid)
                ftest_file << std::setw(16) << std::setprecision(4) << f.cv_score << std::setw(12)
                           << std::setprecision(3) << f.cv_sem;
            else
                ftest_file << std::setw(16) << "n/a" << std::setw(12) << "-";
            ftest_file << std::setw(8) << (isFinalist ? "*" : "") << std::setw(8)
                       << (&f == bestFit ? "<<<" : "") << std::endl;
        }
        ftest_file << "\n1-SE band: [" << std::setprecision(4) << cvMin << ", " << cvThreshold << "]"
                   << std::endl;
        ftest_file << "Decision: simplest function whose CV score is within 1 SE of the best CV;"
                   << std::endl;
        ftest_file << "          equal-parameter ties broken by the better combined data+MC chi2/ndf."
                   << std::endl;
        ftest_file << "Plots: see cross_validation/cv_score_vs_complexity.pdf and cv_folds_*.pdf"
                   << std::endl;

        ftest_file << "\n\nOPTIMAL FUNCTION (selected by "
                   << (bkgSelMode == "dof" ? "smallest-DoF among F-test finalists" : "k-fold CV + 1-SE rule")
                   << "):\n"
                   << std::endl;
        ftest_file << std::string(80, '=') << std::endl;
        ftest_file << "Name:                      " << bestFit->name << std::endl;
        ftest_file << "Type:                      " << bestFit->type << std::endl;
        ftest_file << "Order:                     " << bestFit->order << std::endl;
        ftest_file << "N Parameters:              " << bestFit->npar << std::endl;
        ftest_file << "---" << std::endl;
        ftest_file << "Data Sidebands Chi^2:         " << std::setprecision(2) << bestFit->chi2 << std::endl;
        ftest_file << "Data Sidebands NDF:        " << bestFit->ndf << std::endl;
        ftest_file << "Data Sidebands Chi^2/NDF:     " << std::setprecision(4) << bestFit->chi2_ndf
                   << std::endl;
        ftest_file << "MC Background Chi^2/NDF:      " << std::setprecision(4) << mc_chi2_ndf_best
                   << std::endl;
        ftest_file << "Combined Score:            " << std::setprecision(4) << best_combined_score
                   << (bkgSelMode == "dof" ? "" : " (reference only)") << std::endl;
        if (bestFit->cv_valid)
            ftest_file << "CV Chi^2/bin:                 " << std::setprecision(4) << bestFit->cv_score
                       << " #pm " << std::setprecision(3) << bestFit->cv_sem
                       << (bkgSelMode == "dof" ? "   (reference only)" : "   (SELECTION CRITERION)")
                       << std::endl;
        ftest_file << "---" << std::endl;
        ftest_file << "Spurious Signal:           " << std::setprecision(2) << bestFit->spurious_signal
                   << " events (reference)" << std::endl;
        if (bkgSelMode == "dof") {
            ftest_file << "\nSelection Criterion: smallest DoF (fewest parameters) among the" << std::endl;
            ftest_file << "  F-test finalists; ties broken by the better combined data+MC Chi^2/ndf."
                       << std::endl;
            ftest_file << "  k-fold CV is computed and shown above as a cross-check only." << std::endl;
            ftest_file << "  (Set EASYBJETS_BKG_SELECT=cv to use the CV + 1-SE rule instead.)" << std::endl;
        } else {
            ftest_file
                << "\nSelection Criterion: k-fold cross-validation on the MC background template (sideband region)"
                << std::endl;
            ftest_file << "  + 1-standard-error parsimony rule (simplest function within 1 SE of the"
                       << std::endl;
            ftest_file << "  best out-of-sample CV score). The within-family F-tests above reduce each"
                       << std::endl;
            ftest_file << "  family to one finalist; CV then chooses between the finalists and resists"
                       << std::endl;
            ftest_file << "  overfitting. Combined Chi^2/ndf is reported for reference only; it is NO"
                       << std::endl;
            ftest_file << "  longer the selection criterion. (Set EASYBJETS_BKG_SELECT=dof for the"
                       << std::endl;
            ftest_file << "  smallest-DoF rule; falls back to lowest combined Chi^2/ndf if CV unavailable.)"
                       << std::endl;
        }

        ftest_file << "\n\nFIT PARAMETERS (OPTIMAL FUNCTION):\n" << std::endl;
        ftest_file << std::string(80, '=') << std::endl;
        for (int i = 0; i < bestFit->npar; ++i) {
            ftest_file << "  [" << i << "] = " << std::scientific << std::setprecision(6)
                       << bestFit->func->GetParameter(i) << " #pm " << bestFit->func->GetParError(i)
                       << std::endl;
        }

        ftest_file.close();
        std::cout << "F-test results saved to: " << outputPath << "ftest_results.txt" << std::endl;
    }


    // Save all fit parameters
    std::ofstream params_file((outputPath + "all_fit_parameters.txt").c_str());
    if (params_file.is_open()) {
        params_file << "ALL FIT PARAMETERS COMPARISON" << std::endl;
        params_file << "============================\n" << std::endl;

        for (const auto& fit : allFits) {
            params_file << "\n" << fit.name << " (" << fit.type << " Order " << fit.order << ")" << std::endl;
            params_file << std::string(60, '-') << std::endl;
            params_file << "Chi^2 = " << std::setprecision(4) << fit.chi2 << ", NDF = " << fit.ndf
                        << ", Chi^2/NDF = " << std::setprecision(4) << fit.chi2_ndf << std::endl;
            params_file << "Parameters:" << std::endl;
            for (int j = 0; j < fit.npar; ++j) {
                params_file << "  [" << j << "] = " << std::scientific << std::setprecision(6)
                            << fit.func->GetParameter(j) << " #pm " << fit.func->GetParError(j) << std::endl;
            }
        }
        params_file.close();
        std::cout << "All fit parameters saved to: " << outputPath << "all_fit_parameters.txt" << std::endl;
    }

    // Save comprehensive F-test results for ALL functions "ADDITIONAL"
    std::ofstream ftest_all((outputPath + "ftest_all_results.txt").c_str());
    if (ftest_all.is_open()) {
        ftest_all << "COMPREHENSIVE F-TEST RESULTS FOR ALL FUNCTIONS" << std::endl;
        ftest_all << "===============================================\n" << std::endl;

        ftest_all << "ALL F-TEST COMPARISONS:\n" << std::endl;
        ftest_all << std::setw(25) << "Comparison" << std::setw(15) << "F-statistic" << std::setw(12)
                  << "p-value" << std::setw(15) << "Significant" << std::setw(20) << "Decision" << std::endl;
        ftest_all << std::string(87, '-') << std::endl;

        for (const auto& comp : ftest_results) {
            std::string comparison_name = comp.simpler_model + " vs " + comp.complex_model;
            ftest_all << std::setw(25) << comparison_name << std::setw(15) << std::fixed
                      << std::setprecision(6) << comp.F_statistic << std::setw(12) << std::setprecision(6)
                      << comp.p_value << std::setw(15) << (comp.significant ? "YES" : "NO") << std::setw(20)
                      << (comp.significant ? "Use complex" : "Use simpler") << std::endl;
        }

        ftest_all << "\n\nDETAILED COMPARISON TABLE:\n" << std::endl;
        ftest_all << std::setw(20) << "Simpler Model" << std::setw(20) << "Complex Model" << std::setw(12)
                  << "DeltaChi^2" << std::setw(10) << "Deltadf" << std::setw(15) << "F-stat" << std::setw(12)
                  << "F-crit" << std::setw(12) << "p-value" << std::endl;
        ftest_all << std::string(101, '-') << std::endl;

        for (const auto& comp : ftest_results) {
            ftest_all << std::setw(20) << comp.simpler_model << std::setw(20) << comp.complex_model
                      << std::setw(12) << std::fixed << std::setprecision(2)
                      << (comp.simpler_chi2 - comp.complex_chi2) << std::setw(10) << comp.df1 << std::setw(15)
                      << std::setprecision(6) << comp.F_statistic << std::setw(12) << std::setprecision(4)
                      << comp.F_critical << std::setw(12) << std::setprecision(6) << comp.p_value
                      << std::endl;
        }

        ftest_all.close();
        std::cout << "Comprehensive F-test results saved to: " << outputPath << "ftest_all_results.txt"
                  << std::endl;
    }

    // Save optimal fit parameters
    std::ofstream best_params_file((outputPath + "optimal_fit_parameters.txt").c_str());
    if (best_params_file.is_open()) {
        best_params_file << "OPTIMAL BACKGROUND FIT PARAMETERS" << std::endl;
        best_params_file << "=================================\n" << std::endl;
        best_params_file << "Function Name: " << bestFit->name << std::endl;
        best_params_file << "Function Type: " << bestFit->type << std::endl;
        best_params_file << "Polynomial Order: " << bestFit->order << std::endl;
        best_params_file << "Number of Parameters: " << bestFit->npar << "\n" << std::endl;

        best_params_file << "Fit Quality:" << std::endl;
        best_params_file << "  Chi^2 = " << std::fixed << std::setprecision(2) << bestFit->chi2 << std::endl;
        best_params_file << "  NDF = " << bestFit->ndf << std::endl;
        best_params_file << "  Chi^2/NDF = " << std::setprecision(4) << bestFit->chi2_ndf << std::endl;
        best_params_file << "  Spurious Signal = " << std::setprecision(2) << bestFit->spurious_signal
                         << " events\n"
                         << std::endl;

        best_params_file << "Fit Parameters:" << std::endl;
        best_params_file << std::string(60, '-') << std::endl;
        for (int k = 0; k < bestFit->npar; ++k) {
            best_params_file << "  Parameter [" << k << "]:" << std::endl;
            best_params_file << "    Value = " << std::scientific << std::setprecision(8)
                             << bestFit->func->GetParameter(k) << std::endl;
            best_params_file << "    Error = " << std::scientific << std::setprecision(8)
                             << bestFit->func->GetParError(k) << std::endl;
            best_params_file << "    Relative Error = " << std::fixed << std::setprecision(2)
                             << (bestFit->func->GetParameter(k) != 0
                                     ? 100.0 * bestFit->func->GetParError(k)
                                           / fabs(bestFit->func->GetParameter(k))
                                     : 0)
                             << " %\n"
                             << std::endl;
        }

        best_params_file << "\nFunctional Form:" << std::endl;
        best_params_file << std::string(60, '-') << std::endl;
        if (bestFit->type == "PolyExp") {
            best_params_file << "f(x) = ([0]";
            for (int k = 1; k < bestFit->order; ++k) {
                if (k == 1)
                    best_params_file << " + [1]*x";
                else
                    best_params_file << " + [" << k << "]*x^" << k;
            }
            best_params_file << ") * exp(-[" << bestFit->order << "]*x)" << std::endl;
            best_params_file << "\nWith fitted values:" << std::endl;
            best_params_file << "f(x) = (" << std::scientific << std::setprecision(4)
                             << bestFit->func->GetParameter(0);
            for (int k = 1; k < bestFit->order; ++k) {
                if (k == 1)
                    best_params_file << " + " << bestFit->func->GetParameter(1) << "*x";
                else
                    best_params_file << " + " << bestFit->func->GetParameter(k) << "*x^" << k;
            }
            best_params_file << ") * exp(-" << bestFit->func->GetParameter(bestFit->order) << "*x)"
                             << std::endl;
        } else if (bestFit->type == "Poly") {
            best_params_file << "f(x) = [0]";
            for (int m = 1; m <= bestFit->order; ++m) {
                best_params_file << " + [" << m << "]*x^" << m;
            }
            best_params_file << std::endl;
            best_params_file << "\nWith fitted values:" << std::endl;
            best_params_file << "f(x) = " << std::scientific << std::setprecision(4)
                             << bestFit->func->GetParameter(0);
            for (int m = 1; m <= bestFit->order; ++m) {
                best_params_file << " + " << bestFit->func->GetParameter(m) << "*x^" << m;
            }
            best_params_file << std::endl;
        } else if (bestFit->type == "Exp") {
            best_params_file << "f(x) = exp([0] + [1]*x";
            for (int m = 2; m <= bestFit->order; ++m) {
                best_params_file << " + [" << m << "]*x^" << m;
            }
            best_params_file << ")" << std::endl;
            best_params_file << "\nWith fitted values:" << std::endl;
            best_params_file << "f(x) = exp(" << std::scientific << std::setprecision(4)
                             << bestFit->func->GetParameter(0) << " + " << bestFit->func->GetParameter(1)
                             << "*x";
            for (int m = 2; m <= bestFit->order; ++m) {
                best_params_file << " + " << bestFit->func->GetParameter(m) << "*x^" << m;
            }
            best_params_file << ")" << std::endl;
        }

        best_params_file.close();
        std::cout << "Optimal fit parameters saved to: " << outputPath << "optimal_fit_parameters.txt"
                  << std::endl;
    }

    // Create comprehensive comparison plot
    TCanvas* c_compare = new TCanvas("c_ftest_comparison", "F-Test Function Comparison", 1600, 1200);
    c_compare->Divide(4, 3);

    int pad_idx = 1;
    for (const auto& fit : allFits) {
        c_compare->cd(pad_idx);
        h_sideband->Draw("EP");
        h_sideband->SetTitle((fit.name + ";M^{J} [GeV];Events").c_str());

        fit.func->SetLineColor(kRed);
        fit.func->SetLineWidth(2);
        fit.func->Draw("same");

        TLatex latex;
        latex.SetNDC();
        latex.SetTextSize(0.05);
        latex.DrawLatex(0.65, 0.85, fit.name.c_str());
        latex.SetTextSize(0.04);
        latex.DrawLatex(0.65, 0.78, Form("#chi^{2}/ndf = %.2f", fit.chi2_ndf));

        bool is_optimal = (fit.name == bestFit->name);
        bool is_selected = false;
        for (const auto& sel : selectedFits) {
            if (fit.name == sel->name) {
                is_selected = true;
                break;
            }
        }

        if (is_optimal) {
            latex.SetTextColor(kRed);
            latex.SetTextSize(0.045);
            latex.DrawLatex(0.65, 0.71, "OPTIMAL");
        } else if (is_selected) {
            latex.SetTextColor(kBlue);
            latex.SetTextSize(0.04);
            latex.DrawLatex(0.65, 0.71, "Selected");
        }
        latex.SetTextColor(kBlack);

        pad_idx++;
    }

    SaveCanvasPdf(c_compare, outputPath + "ftest_all_functions_comparison.pdf");
    delete c_compare;

    // Create selected functions comparison
    TCanvas* c_selected = new TCanvas("c_selected_comparison", "Selected Functions Comparison", 1200, 400);
    c_selected->Divide(3, 1);

    for (size_t i = 0; i < selectedFits.size(); ++i) {
        c_selected->cd(i + 1);
        h_sideband->Draw("EP");
        h_sideband->SetTitle((selectedFits[i]->name + ";M^{J} [GeV];Events").c_str());

        selectedFits[i]->func->SetLineColor(kRed);
        selectedFits[i]->func->SetLineWidth(2);
        selectedFits[i]->func->Draw("same");

        TLatex latex;
        latex.SetNDC();
        latex.SetTextSize(0.045);
        latex.DrawLatex(0.5, 0.85, selectedFits[i]->name.c_str());
        latex.SetTextSize(0.04);
        latex.DrawLatex(0.5, 0.80, Form("#chi^{2}/ndf = %.3f", selectedFits[i]->chi2_ndf));
        if (ShowSpuriousSignal()) {
            latex.DrawLatex(0.5, 0.75, Form("Spurious signal = %.1f", selectedFits[i]->spurious_signal));
        }

        if (selectedFits[i] == bestFit) {
            latex.SetTextColor(kRed);
            latex.SetTextSize(0.05);
            latex.DrawLatex(0.5, 0.68, "OPTIMAL");
        }
    }

    SaveCanvasPdf(c_selected, outputPath + "ftest_selected_comparison.pdf");
    delete c_selected;

    std::cout << "\nAll comparison plots saved successfully!\n" << std::endl;

    // Clone the best fit to ensure it persists after temporary fits are deleted
    FitFunction* persistentBestFit = new FitFunction(*bestFit);
    persistentBestFit->func = (TF1*)bestFit->func->Clone((bestFit->name + "_persistent").c_str());

    // Clean up temporary fits (except the cloned one)
    for (auto& fit : allFits) {
        if (fit.func) {
            delete fit.func;
            fit.func = nullptr;
        }
    }

    return persistentBestFit;

    //return bestFit;
}

#endif // MASSFIT_SELECTION_H
