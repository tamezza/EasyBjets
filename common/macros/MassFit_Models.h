/**
 * ============================================================================
 *  MassFit_Models.h - background functions of the sideband fit
 * ============================================================================
 * ===========================================
 * @author:    Thiziri AMEZZA
 * @email:     thiziri.amezza@cern.ch
 * @date:      18 April 2026
 * ===========================================
 *
 * Three families, orders 2 to 5 each:
 *   poly_N     a_0 + a_1 m + ... + a_N m^N                     N+1 parameters
 *   exp_N      exp(a_0 + a_1 m + ... + a_N m^N)                N+1 parameters
 *   polyexp_N  (a_0 + ... + a_(N-1) m^(N-1)) exp(-lambda m)    N+1 parameters
 *
 * Within a family consecutive orders differ by one parameter, which is what the
 * F-test in MassFit_Selection.h needs. polyexp cannot be rewritten as one of the
 * other two (its prefactor can vanish or change sign), so the three families are
 * really different functional forms.
 *
 * Part of the mass fit: it is not run on its own but through Mass_Fit.cxx
 * (how to run and all options: Mass_Fit.h). No option is read here.
 *
 * ==========================================================================
 * @version:  1.0
 * ==========================================================================
 */

#ifndef MASSFIT_MODELS_H
#define MASSFIT_MODELS_H

#include "MassFit_Types.h"

// Create a polynnomial-exponential fit function
inline TF1* CreatePolyExpFit(int order, double xmin, double xmax, const std::string& name) {
    // Create polynomial part of formula
    std::string formula = "(";
    for (int k = 0; k < order; ++k) {
        if (k > 0) formula += " + ";
        formula += "[" + std::to_string(k) + "]";
        if (k == 1)
            formula += "*x";
        else if (k >= 2)
            formula += "*pow(x," + std::to_string(k) + ")";
    }
    // Add exponential part
    formula += ")*exp(-[" + std::to_string(order) + "]*x)";

    TF1* fit = new TF1(name.c_str(), formula.c_str(), xmin, xmax);

    fit->SetParameter(0, 100.0);
    if (order >= 2) fit->SetParameter(1, -0.5);
    for (int k = 2; k < order; ++k) {
        fit->SetParameter(k, TMath::Power(-1, k + 1) * 0.001 / TMath::Power(10, k - 2));
    }
    fit->SetParameter(order, 0.02);
    fit->SetParLimits(order, 0.0, 1.0);

    return fit;
}

// Create a polynomial fit function
inline TF1* CreatePolyFit(int order, double xmin, double xmax, const std::string& name) {
    std::string formula = "[0]";
    for (int i = 1; i <= order; ++i) {
        formula += " + [" + std::to_string(i) + "]*pow(x," + std::to_string(i) + ")";
    }

    TF1* fit = new TF1(name.c_str(), formula.c_str(), xmin, xmax);
    fit->SetParameter(0, 100);
    fit->SetParameter(1, -0.5);
    for (int i = 2; i <= order; ++i) {
        fit->SetParameter(i, TMath::Power(-1, i + 1) * 0.001 / TMath::Power(10, i - 2));
    }
    return fit;
}

// Create an exponential fit function
inline TF1* CreateExpFit(int order, double xmin, double xmax, const std::string& name) {
    std::string formula = "exp([0]";
    for (int i = 1; i <= order; ++i) {
        formula += " + [" + std::to_string(i) + "]*pow(x," + std::to_string(i) + ")";
    }
    formula += ")";

    TF1* fit = new TF1(name.c_str(), formula.c_str(), xmin, xmax);
    fit->SetParameter(0, 4.5);
    fit->SetParameter(1, -0.01);
    for (int i = 2; i <= order; ++i) {
        fit->SetParameter(i, TMath::Power(-1, i) * 1e-6 / TMath::Power(10, i - 2));
    }
    return fit;
}

#endif // MASSFIT_MODELS_H
