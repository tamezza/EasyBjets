/**
 * ============================================================================
 *  MassFit_Types.h - common includes and types of the mass fit
 * ============================================================================
 * ===========================================
 * @author:    Thiziri AMEZZA
 * @email:     thiziri.amezza@cern.ch
 * @date:      18 April 2026
 * ===========================================
 *
 * Included by every mass-fit module. Holds the shared structs (FitFunction,
 * FTestComparison, CVResult), the per-run state (current channel and variant
 * label) and the flags used across the modules.
 *
 * The mass fit is split over these headers:
 *   MassFit_Types.h      this file: shared types, state and common includes
 *   MassFit_Models.h     background functions (poly, exp, poly*exp)
 *   MassFit_Signal.h     fit of the background-subtracted peak (Gaussian, Bukin,
 *                        Crystal Ball, double-sided Crystal Ball)
 *   MassFit_Selection.h  F-test, spurious signal, k-fold CV, choice of the background
 *   MassFit_Plots.h      per-function plots and CV diagnostics
 *   Mass_SB_Plot.h       signal + background overlay with a Data/(S+B) pad
 *   Mass_Fit.h           driver (MassIsolationExtractSig_Generic, Mass_Analysis_Run)
 *
 * Part of the mass fit: it is not run on its own but through Mass_Fit.cxx
 * (how to run and all options: Mass_Fit.h).
 *
 * Option read here:
 *   EASYBJETS_SHOW_SPURIOUS   1 / on to print the spurious-signal test on the
 *                             plots and in the log (default: off)
 *
 * ==========================================================================
 * @version:  1.0
 * ==========================================================================
 */
#ifndef MASSFIT_TYPES_H
#define MASSFIT_TYPES_H

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

// Set once by Mass_Analysis_Run so the helpers do not need a channel argument.
inline std::string& CurrentChannel() {
    static std::string ch = "DB";
    return ch;
}

// Set per variant by the driver (cfg.label), drawn on the fit plots.
inline std::string& CurrentVariantLabel() {
    static std::string lab = "";
    return lab;
}

// EASYBJETS_SHOW_SPURIOUS=1 prints the spurious signal on plots and in the log.
// Defaults to off: (spurious signal is computed and written to the output file but not shown).
// Kept from old studies in case we want to revist this again in the future.
inline bool ShowSpuriousSignal() {
    static bool cached = []() {
        const char* env = std::getenv("EASYBJETS_SHOW_SPURIOUS");
        if (!env || !*env) return false;
        std::string s(env);
        for (auto& c : s) c = (char)std::tolower((unsigned char)c);
        return !(s == "0" || s == "false" || s == "no" || s == "off");
    }();
    return cached;
}

// One candidate background function with its fit quality and cross-validation (CV) score.
struct FitFunction {
    TF1* func;
    std::string name; // e.g. "polyexp_3"
    std::string type; // "PolyExp", "Poly", "Exp"
    int order;
    double chi2;
    int ndf;
    double chi2_ndf;
    int npar;
    double spurious_signal;
    // k-fold cross-validation results (out-of-sample, on the MC template).
    double cv_score; // held-out chi2 per bin
    double cv_sem;   // standard error over the folds
    bool cv_valid;   // false if there were too few bins for CV

    FitFunction(TF1* f, std::string n, std::string t, int o)
        : func(f), name(n), type(t), order(o), chi2(0), ndf(0), chi2_ndf(999), npar(0), spurious_signal(999),
          cv_score(1e9), cv_sem(0.0), cv_valid(false) {}
};

// One nested-model F-test comparison (simpler vs more complex order).
struct FTestComparison {
    std::string simpler_model;
    std::string complex_model;
    int simpler_npar;
    int complex_npar;
    double simpler_chi2;
    double complex_chi2;
    int simpler_ndf;
    int complex_ndf;
    double F_statistic;
    int df1;
    int df2;
    double F_critical;
    double p_value;
    bool significant;
    std::string decision;
};

// Result of a k-fold cross-validation evaluation for one candidate.
struct CVResult {
    double score = 1e9; // total held-out chi2 / total held-out bins (per-bin)
    double sem = 0.0;   // standard error of the per-fold mean (for the 1-SE rule)
    int nTest = 0;      // total held-out bins used
    bool valid = false;
};

#endif // MASSFIT_TYPES_H
