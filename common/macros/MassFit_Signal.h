/**
 * ============================================================================
 *  MassFit_Signal.h - fit of the Z peak
 * ============================================================================
 * ===========================================
 * @author:    Thiziri AMEZZA
 * @email:     thiziri.amezza@cern.ch
 * @date:      18 April 2026
 * ===========================================
 *
 * PerformOptimizedGaussianFit() fits the Z peak (the MC template, or data minus
 * background) in two stages:
 *   1. a first Gaussian fit in a fixed window, [80, 115] GeV unless given;
 *   2. iterative (default): refit in [mu - 2 sigma, mu + 2 sigma] and repeat
 *      with the new mu, sigma until sigma changes by less than 2 % (at most 5
 *      iterations). If sigma keeps flipping between two values, the iterate with
 *      the best chi2/ndf is taken; if no iteration gives a valid fit, one fit in
 *      the 2 sigma window of stage 1 is used.
 *      windowshop: fit in the 1, 2, 2.5 and 3 sigma windows of stage 1 and keep
 *      the best chi2/ndf.
 * With EASYBJETS_FIT_MODEL set, a Bukin, Crystal Ball or double-sided Crystal
 * Ball is also fitted over the window of the best Gaussian. In "all" mode the
 * best chi2/ndf wins; in a single-model mode that model is used whenever its fit
 * succeeded. Downstream code reads the peak as par[1] and the width as par[2].
 *
 * Part of the mass fit: it is not run on its own but through Mass_Fit.cxx
 * (how to run and all options: Mass_Fit.h).
 *
 * Options read here:
 *   EASYBJETS_GAUSS_FIT   iterative (default) | windowshop
 *   EASYBJETS_FIT_MODEL   gauss (default) | bukin | cb | dscb | all
 *                         (parsed in CommonUtils.h; dgauss is not implemented in
 *                         the mass fit and falls back to the Gaussian)
 *
 * ==========================================================================
 * @version:  1.0
 * ==========================================================================
 */

#ifndef MASSFIT_SIGNAL_H
#define MASSFIT_SIGNAL_H

#include "MassFit_Types.h"

std::pair<double, double> FindOptimalFitRange(TH1F* hist, double nSigma = 1.5) {
    int maxBin = hist->GetMaximumBin();
    double peakPosition = hist->GetBinCenter(maxBin);
    double maxContent = hist->GetBinContent(maxBin);
    double halfMax = maxContent / 2.0;

    int leftBin = maxBin;
    while (leftBin > 1 && hist->GetBinContent(leftBin) > halfMax) leftBin--;
    double leftFWHM = hist->GetBinCenter(leftBin);

    int rightBin = maxBin;
    while (rightBin < hist->GetNbinsX() && hist->GetBinContent(rightBin) > halfMax) rightBin++;
    double rightFWHM = hist->GetBinCenter(rightBin);

    const int minBins = 5;
    if ((rightBin - leftBin) < minBins) {
        int binsToAdd = minBins - (rightBin - leftBin);
        leftBin = std::max(1, leftBin - binsToAdd / 2);
        rightBin = std::min(hist->GetNbinsX(), rightBin + binsToAdd / 2 + binsToAdd % 2);
        leftFWHM = hist->GetBinCenter(leftBin);
        rightFWHM = hist->GetBinCenter(rightBin);
    }

    double estimatedSigma = (rightFWHM - leftFWHM) / 2.355;
    double lowerBound = std::max(peakPosition - nSigma * estimatedSigma, hist->GetXaxis()->GetXmin());
    double upperBound = std::min(peakPosition + nSigma * estimatedSigma, hist->GetXaxis()->GetXmax());

    return std::make_pair(lowerBound, upperBound);
}

inline TF1* PerformOptimizedGaussianFit(TH1F* hist, const std::string& fitName,
                                        std::pair<double, double> initialRange = {0, 0},
                                        bool createDiagnosticPlot = false,
                                        const std::string& outputPath = "") {
    if (!hist || hist->GetEntries() < 10) {
        std::cout << "Warning: Insufficient entries for Gaussian fit, using dummy" << std::endl;
        TF1* dummy = new TF1(fitName.c_str(), "gaus", 80, 115);
        dummy->SetParameters(1, 91, 10);
        return dummy;
    }

    std::cout << "\n=== Two-Stage Gaussian Fitting Process ===" << std::endl;
    std::cout << "Histogram: " << hist->GetName() << std::endl;

    // ========================================
    // STAGE 1: Wide initial fit to find peak (80-110 GeV)
    // ========================================
    std::cout << "\nSTAGE 1: Initial wide-range fit (80-115 GeV)" << std::endl;

    double stage1_xmin, stage1_xmax;
    if (initialRange.first == 0 && initialRange.second == 0) {
        stage1_xmin = 80;
        stage1_xmax = 115;
    } else {
        stage1_xmin = initialRange.first;
        stage1_xmax = initialRange.second;
    }

    std::cout << "  Fit range: [" << stage1_xmin << ", " << stage1_xmax << "] GeV" << std::endl;

    TF1* fit_stage1 = new TF1((fitName + "_stage1").c_str(), "gaus", stage1_xmin, stage1_xmax);

    // Get initial estimates from histogram
    int maxBin = hist->GetMaximumBin();
    double amplitude = hist->GetBinContent(maxBin);
    double mean_initial = hist->GetBinCenter(maxBin);
    double rms = hist->GetRMS();

    fit_stage1->SetParameter(0, amplitude);
    fit_stage1->SetParameter(1, mean_initial);
    fit_stage1->SetParameter(2, rms * 0.5);
    fit_stage1->SetParLimits(0, 0, amplitude * 2);
    fit_stage1->SetParLimits(1, stage1_xmin, stage1_xmax);
    fit_stage1->SetParLimits(2, 0.1, (stage1_xmax - stage1_xmin) / 2);

    hist->Fit(fit_stage1, "RQ0");

    double mean_stage1 = fit_stage1->GetParameter(1);
    double sigma_stage1 = fit_stage1->GetParameter(2);
    double chi2_stage1 = (fit_stage1->GetNDF() > 0) ? fit_stage1->GetChisquare() / fit_stage1->GetNDF() : 999;

    std::cout << "  Results: mu = " << mean_stage1 << " GeV, sigma = " << sigma_stage1 << " GeV" << std::endl;
    std::cout << "  Chi^2/ndf = " << chi2_stage1 << std::endl;

    // ========================================
    // STAGE 2: Refined fit around peak
    // ========================================
    // Selectable via EASYBJETS_GAUSS_FIT env var:
    //   "iterative" (default) - 2sigma self-consistent window:
    //         fit in [mu_1-2sigma_1, mu_1+2sigma_1], extract sigma_2, refit in [mu_2-2sigma_2, mu_2+2sigma_2],
    //         iterate until |dsigma/sigma| < 2% or MAX_ITERATIONS hit.
    //         Removes the bias toward narrow windows that window-shopping has.
    //
    //   "windowshop" - legacy behavior: fit in four candidate windows
    //         (1sigma, 2sigma, 2.5sigma, 3sigma), pick the one with lowest Chi^2/ndf.
    //
    // If an iterative fit fails to converge or produces non-finite params,
    // falls back to a fixed 2sigma window.
    //
    // Tunable parameters:
    constexpr double kNSigmaWindow = 2.0;   // width of the refit window
    constexpr double kConvTolerance = 0.02; // |dsigma|/sigma threshold (2%)
    constexpr int kMaxIterations = 5;       // max iteration guard

    // Read mode from env var once per call.
    std::string fitMode = "iterative";
    if (const char* env = std::getenv("EASYBJETS_GAUSS_FIT")) {
        std::string s(env);
        for (auto& c : s) c = (char)std::tolower((unsigned char)c);
        if (s == "windowshop" || s == "window" || s == "legacy")
            fitMode = "windowshop";
        else if (s == "iterative" || s == "iter")
            fitMode = "iterative";
        else
            std::cerr << "[GaussFit] WARN: unknown EASYBJETS_GAUSS_FIT='" << env << "', using iterative\n";
    }

    TF1* bestFit = nullptr;
    double bestChi2NDF = 1e9;
    std::vector<TF1*> candidateFits; // kept for the diagnostic canvas

    // ===========================================================
    // BRANCH A: iterative self-consistent 2sigma window
    // ===========================================================
    if (fitMode == "iterative") {
        std::cout << "\nSTAGE 2: iterative " << kNSigmaWindow << "sigma self-consistent refit" << std::endl;

        double mu_cur = mean_stage1;
        double sigma_cur = sigma_stage1;
        TF1* fit_cur = nullptr;
        bool converged = false;

        // Oscillation detection: store each iterate's sigma and Chi^2/ndf so we can
        // identify a flip-flop and pick the best iterate instead of the last.
        struct IterRecord {
            TF1* f;
            double sigma;
            double chi2ndf;
        };
        std::vector<IterRecord> iter_history;
        int oscillationCount = 0;

        for (int iter = 1; iter <= kMaxIterations; ++iter) {
            double lo = std::max(hist->GetXaxis()->GetXmin(), mu_cur - kNSigmaWindow * sigma_cur);
            double hi = std::min(hist->GetXaxis()->GetXmax(), mu_cur + kNSigmaWindow * sigma_cur);
            if (lo >= hi) {
                std::cout << "  iter " << iter << ": window collapsed, stopping iteration" << std::endl;
                break;
            }

            std::string name = fitName + "_iter" + std::to_string(iter);
            TF1* f = new TF1(name.c_str(), "gaus", lo, hi);
            f->SetParameter(0, fit_stage1->GetParameter(0));
            f->SetParameter(1, mu_cur);
            f->SetParameter(2, sigma_cur);
            f->SetParLimits(0, 0, fit_stage1->GetParameter(0) * 3);
            f->SetParLimits(1, lo, hi);
            f->SetParLimits(2, 0.1, (hi - lo) / 2);

            hist->Fit(f, "RQ0");

            double mu_new = f->GetParameter(1);
            double sigma_new = f->GetParameter(2);
            double chi2ndf = (f->GetNDF() > 0) ? f->GetChisquare() / f->GetNDF() : 1e9;

            // Sanity-check the fit result. Non-finite or out-of-physical-range values
            // indicate a failed fit - abandon iteration, keep fit_cur (last good) if any.
            bool ok = std::isfinite(mu_new) && std::isfinite(sigma_new) && sigma_new > 0.1 && sigma_new < 50
                      && mu_new > hist->GetXaxis()->GetXmin() && mu_new < hist->GetXaxis()->GetXmax();
            if (!ok) {
                std::cout << "  iter " << iter << ": invalid result (mu=" << mu_new << ", sigma=" << sigma_new
                          << "), abandoning iteration" << std::endl;
                delete f;
                break;
            }

            candidateFits.push_back(f);
            std::cout << "  iter " << iter << ": window=[" << std::fixed << std::setprecision(2) << lo << ","
                      << hi << "]"
                      << "  mu=" << std::setprecision(3) << mu_new << "  sigma=" << sigma_new
                      << "  Chi^2/ndf=" << std::setprecision(4) << chi2ndf << std::endl;

            iter_history.push_back({f, sigma_new, chi2ndf});
            if (iter_history.size() >= 2) {
                const auto& prev = iter_history[iter_history.size() - 2];
                double rel = (prev.sigma > 1e-12) ? std::abs((sigma_new - prev.sigma) / prev.sigma)
                                                  : std::abs(sigma_new - prev.sigma);
                if (rel > 0.05)
                    ++oscillationCount;
                else
                    oscillationCount = 0;
            }

            // Convergence test on sigma (sigma drives the window).
            double rel_change = std::fabs(sigma_new - sigma_cur) / std::max(sigma_cur, 1e-6);
            if (iter > 1 && rel_change < kConvTolerance) {
                std::cout << "  converged (|dsigma|/sigma = " << std::setprecision(4) << rel_change << " < "
                          << kConvTolerance << ")" << std::endl;
                fit_cur = f;
                mu_cur = mu_new;
                sigma_cur = sigma_new;
                converged = true;
                break;
            }

            // Test for oscillation: if sigma changes by more than 5% in two consecutive
            //  iterations, we're probably bouncing between two minima. Take the best
            // iterate by Chi^2/ndf instead of continuing to the last one.
            if (iter >= 3 && oscillationCount >= 2) {
                std::cout << "  oscillation detected; selecting best Chi^2/ndf iterate" << std::endl;
                double bestIterChi2ndf = 1e9;
                TF1* bestIter = nullptr;
                for (const auto& rec : iter_history) {
                    if (rec.chi2ndf < bestIterChi2ndf) {
                        bestIterChi2ndf = rec.chi2ndf;
                        bestIter = rec.f;
                    }
                }
                if (bestIter) {
                    fit_cur = bestIter;
                    mu_cur = bestIter->GetParameter(1);
                    sigma_cur = bestIter->GetParameter(2);
                    converged = true;
                    std::cout << "  selected iterate with Chi^2/ndf=" << std::setprecision(4)
                              << bestIterChi2ndf << std::endl;
                }
                break;
            }

            fit_cur = f;
            mu_cur = mu_new;
            sigma_cur = sigma_new;
        }

        if (converged && fit_cur) {
            bestFit = fit_cur;
            bestChi2NDF = (bestFit->GetNDF() > 0) ? bestFit->GetChisquare() / bestFit->GetNDF() : 1e9;
        } else if (fit_cur) {
            // Iteration exhausted MAX_ITERATIONS without converging.
            // take the last iterate, but flag non-convergence in the log.
            std::cout << "  did not converge in " << kMaxIterations
                      << " iterations, using last iterate (may indicate instability)" << std::endl;
            bestFit = fit_cur;
            bestChi2NDF = (bestFit->GetNDF() > 0) ? bestFit->GetChisquare() / bestFit->GetNDF() : 1e9;
        } else {
            // No iteration produced a valid fit. Fall back to a single fixed-window
            // fit using Stage 1's (mu,sigma) - this is the safest landing.
            std::cout << "  iteration produced no valid fit; falling back to fixed " << kNSigmaWindow
                      << "sigma single-shot" << std::endl;
            double lo = std::max(hist->GetXaxis()->GetXmin(), mean_stage1 - kNSigmaWindow * sigma_stage1);
            double hi = std::min(hist->GetXaxis()->GetXmax(), mean_stage1 + kNSigmaWindow * sigma_stage1);
            TF1* fb = new TF1((fitName + "_fallback").c_str(), "gaus", lo, hi);
            fb->SetParameters(fit_stage1->GetParameters());
            fb->SetParLimits(0, 0, fit_stage1->GetParameter(0) * 3);
            fb->SetParLimits(1, lo, hi);
            fb->SetParLimits(2, 0.1, (hi - lo) / 2);
            hist->Fit(fb, "RQ0");
            candidateFits.push_back(fb);
            bestFit = fb;
            bestChi2NDF = (fb->GetNDF() > 0) ? fb->GetChisquare() / fb->GetNDF() : 1e9;
        }
    }
    // ===========================================================
    // BRANCH B: legacy window-shopping (kept for A/B comparison)
    // ===========================================================
    else {
        std::cout
            << "\nSTAGE 2: legacy window-shop (1sigma, 2sigma, 2.5sigma, 3sigma - pick lowest Chi^2/ndf)"
            << std::endl;

        std::vector<double> nSigmaOptions = {1.0, 2.0, 2.5, 3.0};
        for (double nSigma : nSigmaOptions) {
            double stage2_xmin = std::max(hist->GetXaxis()->GetXmin(), mean_stage1 - nSigma * sigma_stage1);
            double stage2_xmax = std::min(hist->GetXaxis()->GetXmax(), mean_stage1 + nSigma * sigma_stage1);
            if (stage2_xmin >= stage2_xmax) continue;

            std::cout << "  Testing " << nSigma << "sigma range: [" << std::fixed << std::setprecision(2)
                      << stage2_xmin << ", " << stage2_xmax << "] GeV" << std::endl;

            TF1* fit_candidate = new TF1((fitName + "_" + std::to_string(nSigma) + "sigma").c_str(), "gaus",
                                         stage2_xmin, stage2_xmax);
            fit_candidate->SetParameters(fit_stage1->GetParameters());
            fit_candidate->SetParLimits(0, 0, fit_stage1->GetParameter(0) * 3);
            fit_candidate->SetParLimits(1, stage2_xmin, stage2_xmax);
            fit_candidate->SetParLimits(2, 0.1, (stage2_xmax - stage2_xmin) / 2);

            hist->Fit(fit_candidate, "RQ0");

            double chi2ndf =
                (fit_candidate->GetNDF() > 0) ? fit_candidate->GetChisquare() / fit_candidate->GetNDF() : 1e9;

            std::cout << "    -> mu = " << std::setprecision(3) << fit_candidate->GetParameter(1) << " +/- "
                      << fit_candidate->GetParError(1) << " GeV" << std::endl;
            std::cout << "    -> sigma = " << std::setprecision(3) << fit_candidate->GetParameter(2)
                      << " +/- " << fit_candidate->GetParError(2) << " GeV" << std::endl;
            std::cout << "    -> Chi^2/ndf = " << std::setprecision(4) << chi2ndf << std::endl;

            candidateFits.push_back(fit_candidate);

            if (chi2ndf < bestChi2NDF) {
                bestChi2NDF = chi2ndf;
                bestFit = fit_candidate;
            }
        }
    }

    // ===========================================================
    // STAGE 2c: optional alternative signal candidates (EASYBJETS_FIT_MODEL)
    // ===========================================================
    // When EASYBJETS_FIT_MODEL is set, try alternative fit functions:
    //  Gaussian, Bukin, Crystal Ball and double-sided CB, the best chi2/ndf wins
    //   bukin / cb / dscb : only that shape is reported, even if a Gaussian
    //   double-gauss (dgauss) is NOT supported here (only in the response fit).
    FitUtils::FitModelSet fit_mode = FitUtils::CurrentFitModelSet();
    const bool do_alt_models =
        (fit_mode == FitUtils::FitModelSet::AllModels || fit_mode == FitUtils::FitModelSet::BukinOnly
         || fit_mode == FitUtils::FitModelSet::CrystalBallOnly
         || fit_mode == FitUtils::FitModelSet::DoubleCBOnly);
    const bool want_bukin =
        (fit_mode == FitUtils::FitModelSet::AllModels || fit_mode == FitUtils::FitModelSet::BukinOnly);
    const bool want_cb =
        (fit_mode == FitUtils::FitModelSet::AllModels || fit_mode == FitUtils::FitModelSet::CrystalBallOnly);
    const bool want_dscb =
        (fit_mode == FitUtils::FitModelSet::AllModels || fit_mode == FitUtils::FitModelSet::DoubleCBOnly);

    TF1* alt_winner = nullptr;

    if (bestFit && do_alt_models) {
        double seed_lo = bestFit->GetXmin();
        double seed_hi = bestFit->GetXmax();
        double seed_mu = bestFit->GetParameter(1);
        double seed_sig = bestFit->GetParameter(2);
        double seed_amp = bestFit->GetParameter(0);

        // Bukin candidate
        if (want_bukin) {
            TF1* bukin = new TF1((fitName + "_bukin").c_str(), FitUtils::BukinFunction, seed_lo, seed_hi, 4);
            bukin->SetParameters(seed_amp, seed_mu, seed_sig, 0.0);
            bukin->SetParLimits(2, seed_sig * 0.3, seed_sig * 3.0);
            bukin->SetParLimits(3, -0.5, 0.5);
            int status = hist->Fit(bukin, "RQ0");
            double bu_mu = bukin->GetParameter(1);
            double bu_sig = bukin->GetParameter(2);
            double bu_chi2ndf = (bukin->GetNDF() > 0) ? bukin->GetChisquare() / bukin->GetNDF() : 1e9;
            bool bu_ok = (status == 0) && std::isfinite(bu_mu) && std::isfinite(bu_sig) && bu_sig > 0.1
                         && bu_sig < 50 && std::isfinite(bu_chi2ndf);
            std::cout << "  bukin:  mu=" << std::setprecision(3) << bu_mu << "  sigma=" << bu_sig
                      << "  Chi^2/ndf=" << std::setprecision(4) << bu_chi2ndf << (bu_ok ? "" : " [INVALID]")
                      << std::endl;
            if (bu_ok) {
                candidateFits.push_back(bukin);
                alt_winner = bukin;
                if (bu_chi2ndf < bestChi2NDF) {
                    bestFit = bukin;
                    bestChi2NDF = bu_chi2ndf;
                    std::cout << "  -> bukin wins selection" << std::endl;
                }
            } else {
                delete bukin;
            }
        }

        // Crystal Ball candidate
        if (want_cb) {
            TF1* cb =
                new TF1((fitName + "_cball").c_str(), FitUtils::CrystalBallFunction, seed_lo, seed_hi, 5);
            cb->SetParameters(seed_amp, seed_mu, seed_sig, 1.5, 3.0);
            cb->SetParLimits(2, seed_sig * 0.3, seed_sig * 3.0);
            cb->SetParLimits(3, 0.5, 5.0);
            cb->SetParLimits(4, 1.1, 20.0);
            int status = hist->Fit(cb, "RQ0");
            double cb_mu = cb->GetParameter(1);
            double cb_sig = cb->GetParameter(2);
            double cb_chi2ndf = (cb->GetNDF() > 0) ? cb->GetChisquare() / cb->GetNDF() : 1e9;
            bool cb_ok = (status == 0) && std::isfinite(cb_mu) && std::isfinite(cb_sig) && cb_sig > 0.1
                         && cb_sig < 50 && std::isfinite(cb_chi2ndf);
            std::cout << "  cball:  mu=" << std::setprecision(3) << cb_mu << "  sigma=" << cb_sig
                      << "  Chi^2/ndf=" << std::setprecision(4) << cb_chi2ndf << (cb_ok ? "" : " [INVALID]")
                      << std::endl;
            if (cb_ok) {
                candidateFits.push_back(cb);
                alt_winner = cb; // overrides bukin if both succeeded
                if (cb_chi2ndf < bestChi2NDF) {
                    bestFit = cb;
                    bestChi2NDF = cb_chi2ndf;
                    std::cout << "  -> crystal ball wins selection" << std::endl;
                }
            } else {
                delete cb;
            }
        }

        // Double-sided Crystal Ball candidate (Gaussian core + two power-law tails)
        if (want_dscb) {
            TF1* dscb = new TF1((fitName + "_dscb").c_str(), FitUtils::DoubleCrystalBallFunction, seed_lo,
                                seed_hi, 7);
            dscb->SetParameters(seed_amp, seed_mu, seed_sig, 1.5, 3.0, 1.5, 3.0);
            dscb->SetParLimits(2, seed_sig * 0.3, seed_sig * 3.0);
            dscb->SetParLimits(3, 0.3, 5.0);  // alphaL
            dscb->SetParLimits(4, 1.1, 50.0); // nL
            dscb->SetParLimits(5, 0.3, 5.0);  // alphaR
            dscb->SetParLimits(6, 1.1, 50.0); // nR
            int status = hist->Fit(dscb, "RQ0");
            double d_mu = dscb->GetParameter(1);
            double d_sig = dscb->GetParameter(2);
            double d_chi2ndf = (dscb->GetNDF() > 0) ? dscb->GetChisquare() / dscb->GetNDF() : 1e9;
            bool d_ok = (status == 0) && std::isfinite(d_mu) && std::isfinite(d_sig) && d_sig > 0.1
                        && d_sig < 50 && std::isfinite(d_chi2ndf);
            std::cout << "  dscb:   mu=" << std::setprecision(3) << d_mu << "  sigma=" << d_sig
                      << "  Chi^2/ndf=" << std::setprecision(4) << d_chi2ndf << (d_ok ? "" : " [INVALID]")
                      << std::endl;
            if (d_ok) {
                candidateFits.push_back(dscb);
                alt_winner = dscb; // overrides bukin/cb if it also succeeded
                if (d_chi2ndf < bestChi2NDF) {
                    bestFit = dscb;
                    bestChi2NDF = d_chi2ndf;
                    std::cout << "  -> double-sided CB wins selection" << std::endl;
                }
            } else {
                delete dscb;
            }
        }

        // For BukinOnly / CrystalBallOnly / DoubleCBOnly: force the alt-model
        // result to be the final reported fit even if a Gaussian candidate has a
        // better Chi^2/ndf.
        if (fit_mode == FitUtils::FitModelSet::BukinOnly || fit_mode == FitUtils::FitModelSet::CrystalBallOnly
            || fit_mode == FitUtils::FitModelSet::DoubleCBOnly) {
            if (alt_winner) {
                bestFit = alt_winner;
                bestChi2NDF =
                    (alt_winner->GetNDF() > 0) ? alt_winner->GetChisquare() / alt_winner->GetNDF() : 1e9;
                const char* fmName = (fit_mode == FitUtils::FitModelSet::BukinOnly) ? "bukin"
                                     : (fit_mode == FitUtils::FitModelSet::CrystalBallOnly)
                                         ? "crystal ball"
                                         : "double-sided CB";
                const char* fmTag = (fit_mode == FitUtils::FitModelSet::BukinOnly)         ? "bukin"
                                    : (fit_mode == FitUtils::FitModelSet::CrystalBallOnly) ? "cb"
                                                                                           : "dscb";
                std::cout << "  -> forced selection: " << fmName << " (FIT_MODEL=" << fmTag << ")"
                          << std::endl;
            } else {
                std::cerr << "[FitModel] WARN: requested alt-model fit failed; "
                             "falling back to best Gaussian\n";
            }
        }
    }

    if (!bestFit) {
        std::cout << "Warning: No valid refined fit found, using Stage 1 fit" << std::endl;
        bestFit = fit_stage1;
        bestChi2NDF = chi2_stage1;
    }

    std::cout << "\n*** FINAL FIT SELECTED ***" << std::endl;
    std::cout << "  Best window: " << bestFit->GetName() << std::endl;
    std::cout << "  Range: [" << std::fixed << std::setprecision(2) << bestFit->GetXmin() << ", "
              << bestFit->GetXmax() << "] GeV" << std::endl;
    std::cout << "  mu = " << std::setprecision(3) << bestFit->GetParameter(1) << " +/- "
              << bestFit->GetParError(1) << " GeV" << std::endl;
    std::cout << "  sigma = " << std::setprecision(3) << bestFit->GetParameter(2) << " +/- "
              << bestFit->GetParError(2) << " GeV" << std::endl;
    std::cout << "  Chi^2/ndf = " << std::setprecision(4) << bestChi2NDF << std::endl;
    std::cout << "=========================================\n" << std::endl;

    // Validate parameters
    double final_mean = bestFit->GetParameter(1);
    double final_sigma = bestFit->GetParameter(2);
    double final_amp = bestFit->GetParameter(0);

    if (!std::isfinite(final_mean) || !std::isfinite(final_sigma) || !std::isfinite(final_amp)
        || final_mean <= 0 || final_sigma <= 0 || final_amp <= 0 || final_mean > 200 || final_sigma > 50) {
        std::cout << "Warning: Fit produced invalid parameters, using safe defaults" << std::endl;
        bestFit->SetParameter(0, hist->GetMaximum());
        bestFit->SetParameter(1, hist->GetMean());
        bestFit->SetParameter(2, hist->GetRMS() * 0.5);
    }


    // ============================================================
    // Diagnostic canvas: stage 1, each candidate, overlay, residuals
    // ============================================================
    if (createDiagnosticPlot && outputPath != "") {
        TCanvas* c_diag =
            new TCanvas(("c_" + fitName + "_diag").c_str(), "Two-Stage Gaussian Fit Diagnostic", 1800, 1200);
        c_diag->Divide(3, 2);

        // Panel 1: Stage 1 fit (70-130 GeV)
        c_diag->cd(1);
        gPad->SetLeftMargin(0.15);
        hist->Draw("EP");
        hist->SetTitle("Stage 1: Wide Range (70-130 GeV)");
        hist->GetYaxis()->SetTitleOffset(1.4);
        fit_stage1->SetLineColor(kBlue);
        fit_stage1->SetLineWidth(2);
        fit_stage1->Draw("same");

        TBox* box1 = new TBox(fit_stage1->GetXmin(), 0, fit_stage1->GetXmax(), hist->GetMaximum() * 0.3);
        box1->SetFillColorAlpha(kBlue, 0.15);
        box1->Draw("same");
        hist->Draw("EP same");
        fit_stage1->Draw("same");

        TLatex lat1;
        lat1.SetNDC();
        lat1.SetTextSize(0.04);
        lat1.SetTextColor(kBlue);
        lat1.DrawLatex(0.18, 0.87, "INITIAL FIT RANGE");
        lat1.SetTextColor(kBlack);
        lat1.SetTextSize(0.032);
        lat1.DrawLatex(0.18, 0.82, Form("Mass Range: [%.1f, %.1f] GeV", stage1_xmin, stage1_xmax));
        lat1.DrawLatex(0.18, 0.77, Form("#mu = %.2f GeV", mean_stage1));
        lat1.DrawLatex(0.18, 0.72, Form("#sigma = %.2f GeV", sigma_stage1));
        lat1.DrawLatex(0.18, 0.67, Form("#chi^{2}/ndf = %.2f", chi2_stage1));

        // Panels 2-4: Test all sigma windows (2sigma, 2.5sigma, 3sigma)
        int panel = 2;
        for (size_t i = 0; i < candidateFits.size(); ++i) {
            c_diag->cd(panel++);
            gPad->SetLeftMargin(0.15);

            TH1F* h_copy = (TH1F*)hist->Clone(Form("h_temp_%d", panel));
            h_copy->SetDirectory(0);
            h_copy->Draw("EP");
            h_copy->SetTitle(Form("Stage 2: Test %s", candidateFits[i]->GetName()));
            h_copy->GetYaxis()->SetTitleOffset(1.4);

            candidateFits[i]->SetLineColor(kMagenta);
            candidateFits[i]->SetLineWidth(2);
            candidateFits[i]->Draw("same");

            TBox* box = new TBox(candidateFits[i]->GetXmin(), 0, candidateFits[i]->GetXmax(),
                                 hist->GetMaximum() * 0.3);
            box->SetFillColorAlpha(kMagenta, 0.15);
            box->Draw("same");
            h_copy->Draw("EP same");
            candidateFits[i]->Draw("same");

            // Highlight the best fit
            if (candidateFits[i] == bestFit) {
                TLatex bestLabel;
                bestLabel.SetNDC();
                bestLabel.SetTextSize(0.04);
                bestLabel.SetTextColor(kRed);
                bestLabel.DrawLatex(0.18, 0.87, "*** SELECTED ***");
            }

            TLatex lat;
            lat.SetNDC();
            lat.SetTextSize(0.032);
            lat.DrawLatex(0.18, 0.82,
                          Form("Mass Range: [%.1f, %.1f] GeV", candidateFits[i]->GetXmin(),
                               candidateFits[i]->GetXmax()));
            lat.DrawLatex(0.18, 0.77,
                          Form("#mu = %.2f #pm %.2f", candidateFits[i]->GetParameter(1),
                               candidateFits[i]->GetParError(1)));
            lat.DrawLatex(0.18, 0.72,
                          Form("#sigma = %.2f #pm %.2f", candidateFits[i]->GetParameter(2),
                               candidateFits[i]->GetParError(2)));
            double chi2ndf = (candidateFits[i]->GetNDF() > 0)
                                 ? candidateFits[i]->GetChisquare() / candidateFits[i]->GetNDF()
                                 : 0;
            lat.DrawLatex(0.18, 0.67, Form("#chi^{2}/ndf = %.3f", chi2ndf));
        }

        // Panel 5: Comparison of all fits
        c_diag->cd(5);
        gPad->SetLeftMargin(0.15);
        TH1F* h_compare = (TH1F*)hist->Clone("h_compare");
        h_compare->SetDirectory(0);
        h_compare->Draw("EP");
        h_compare->SetTitle("Comparison: All Fit Windows");
        h_compare->GetYaxis()->SetTitleOffset(1.4);

        fit_stage1->SetLineColor(kBlue);
        fit_stage1->SetLineStyle(2);
        fit_stage1->SetLineWidth(2);
        fit_stage1->Draw("same");

        int color_idx = 0;
        std::vector<int> colors = {kGreen + 2, kOrange + 1, kViolet};
        for (auto* fit : candidateFits) {
            fit->SetLineColor(colors[color_idx % colors.size()]);
            fit->SetLineWidth(2);
            if (fit == bestFit) {
                fit->SetLineWidth(3);
                fit->SetLineStyle(1);
            } else {
                fit->SetLineStyle(1);
            }
            fit->Draw("same");
            color_idx++;
        }

        TLegend* leg = new TLegend(0.18, 0.65, 0.88, 0.9);
        leg->SetBorderSize(0);
        leg->SetTextSize(0.028);
        leg->AddEntry(hist, "Data", "lep");
        //Stage 1
        leg->AddEntry(fit_stage1, Form("Initial Fit Range (70-130 GeV): #chi^{2}/ndf=%.2f", chi2_stage1),
                      "l");
        color_idx = 0;

        for (size_t i = 0; i < candidateFits.size(); ++i) {
            auto* fit = candidateFits[i];
            double chi2ndf = (fit->GetNDF() > 0) ? fit->GetChisquare() / fit->GetNDF() : 0;

            // Extract the iteration number from the fit name
            std::string fullName = fit->GetName();
            std::string sigmaLabel;
            size_t underscore = fullName.rfind('_');
            std::string tail = (underscore != std::string::npos) ? fullName.substr(underscore + 1) : fullName;
            if (tail.rfind("iter", 0) == 0) {
                // "iter3" -> "iter 3"
                sigmaLabel = "iter " + tail.substr(4);
            } else if (tail.find("sigma") != std::string::npos) {
                std::string num = tail.substr(0, tail.find("sigma"));
                while (!num.empty() && num.back() == '0') num.pop_back();
                if (!num.empty() && num.back() == '.') num.pop_back();
                sigmaLabel = num + "#sigma";
            } else if (tail == "fallback") {
                sigmaLabel = "fallback";
            } else {
                sigmaLabel = tail;
            }

            std::string label = Form("%s: #chi^{2}/ndf=%.2f", sigmaLabel.c_str(), chi2ndf);
            TLegendEntry* entry = leg->AddEntry(fit, label.c_str(), "l");
            if (fit == bestFit) {
                entry->SetLabel(Form("%s #color[2]{*** SELECTED ***}", label.c_str()));
            }
            color_idx++;
        }

        leg->Draw();

        // Panel 6: Residuals for best fit
        c_diag->cd(6);
        gPad->SetLeftMargin(0.15);
        TH1F* h_residuals = (TH1F*)hist->Clone((fitName + "_residuals").c_str());
        h_residuals->SetDirectory(0);
        for (int i = 1; i <= h_residuals->GetNbinsX(); ++i) {
            double x = h_residuals->GetBinCenter(i);
            if (x >= bestFit->GetXmin() && x <= bestFit->GetXmax()) {
                double data = h_residuals->GetBinContent(i);
                double fit_val = bestFit->Eval(x);
                double residual = 0;
                if (fit_val > 0 && data >= 0) {
                    residual = (data - fit_val) / sqrt(fit_val);
                } else if (fit_val > 0) {
                    residual = -sqrt(fit_val);
                }
                h_residuals->SetBinContent(i, residual);
            } else {
                h_residuals->SetBinContent(i, 0);
            }
        }
        h_residuals->SetTitle("Normalized Residuals (Best Fit);M^{J} [GeV];(Data - Fit) / #sqrt{Fit}");
        h_residuals->GetYaxis()->SetTitleOffset(1.4);
        h_residuals->SetLineColor(kBlack);
        h_residuals->SetMinimum(-5);
        h_residuals->SetMaximum(5);
        h_residuals->Draw("EP");

        TLine* zero_line =
            new TLine(h_residuals->GetXaxis()->GetXmin(), 0, h_residuals->GetXaxis()->GetXmax(), 0);
        zero_line->SetLineStyle(2);
        zero_line->SetLineColor(kRed);
        zero_line->SetLineWidth(2);
        zero_line->Draw("same");

        TLine* upper_2sigma =
            new TLine(h_residuals->GetXaxis()->GetXmin(), 2, h_residuals->GetXaxis()->GetXmax(), 2);
        upper_2sigma->SetLineStyle(3);
        upper_2sigma->SetLineColor(kBlue);
        upper_2sigma->Draw("same");

        TLine* lower_2sigma =
            new TLine(h_residuals->GetXaxis()->GetXmin(), -2, h_residuals->GetXaxis()->GetXmax(), -2);
        lower_2sigma->SetLineStyle(3);
        lower_2sigma->SetLineColor(kBlue);
        lower_2sigma->Draw("same");

        SaveCanvasPdf(c_diag, outputPath + fitName + "_two_stage_diagnostic.pdf");
        delete c_diag;
        delete h_residuals;
        delete h_compare;
    }

    // Clean up non-selected fits
    for (auto* fit : candidateFits) {
        if (fit != bestFit) {
            delete fit;
        }
    }
    if (fit_stage1 != bestFit) {
        delete fit_stage1;
    }

    return bestFit;
}

#endif // MASSFIT_SIGNAL_H
