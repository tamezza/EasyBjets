/**
 * ============================================================================
 *  Response_Fit_Core.h - functions shared by the response fits
 * ============================================================================
 * ===========================================
 * @author:    Thiziri AMEZZA
 * @email:     thiziri.amezza@cern.ch
 * @date:      10 April 2026
 * ===========================================
 *
 * Functions used by both response fits (Response_Fit_Inclusive.h and
 * Response_Fit_pTbins_Mcuts.h):
 * 
 *   DataDrivenSeed             seed window of the iterative Gaussian fit
 *   ComputePerSampleJER        in-situ resolution and closure from the fitted
 *                              widths, for one sample (signal, bkg or combined)
 *   FillResponseBinningHists   MC pass that fills the histograms used to choose
 *                              the response binning
 *                              EASYBJETS_RESPONSE_BINNING_MODE=dynamic; the truth
 *                              probe is filled only when EASYBJETS_BINNING_TRUTH
 *                              is auto or algorithm (see BinningConfig.h)
 * 
 * Options that act here:
 * 
 *   EASYBJETS_MASS_CUT                 applied in FillResponseBinningHists exactly
 *                                      as in the fit loops (MassCutConfig.h)
 *   EASYBJETS_RESPONSE_BINNING_MODE    static (default) | dynamic: whether the
 *                                      binning probe runs at all
 *   EASYBJETS_BINNING_TRUTH            match (default) | auto | algorithm: whether
 *                                      the truth probe is filled
 *
 * JER closure, with sigma the fitted width of the balance R = pT_jet / pT_ref:
 * 
 *   sigma_insitu    = sqrt(sigma_data^2 - sigma_truth^2)
 *   sigma^MC_insitu = sqrt(sigma_RDB^MC^2 - sigma_truth^2)
 *   sigma^MC        = width of pT_reco / pT_truth
 *   Delta_closure   = (sigma^MC - sigma^MC_insitu) / sigma_insitu
 * 
 * The errors are propagated as independent; sigma_truth enters both square
 * roots, so its correlation between them is not taken into account.
 *
 * Not run on its own; see the two fit headers for how to run them.
 *
 * ==========================================================================
 * @version:  1.0
 * ==========================================================================
 */

#ifndef RESPONSE_FIT_CORE_H
#define RESPONSE_FIT_CORE_H

#include "Response_Fit_Common.h" // RDBFitResults + channel/pT-bin accessors
#include "CommonUtils.h"         // Log::, Form, FitUtils::, ...
#include "VariantConfig.h"       // VariantConfig (branch names, channel, maxTruthDR)
#include "MassCutConfig.h"       // MassCutConfig::Reject (optional mass window)
#include "TTree.h"
#include "TH1F.h"
#include <vector>
#include <string>
#include <iostream>
#include <cmath>
#include <utility>

/* ----------------------------------------------------------------------------
 * Fit window
 * -------------------------------------------------------------------------- */

// Starting window of the iterative Gaussian fit, taken from the histogram:
//   [mean - nSig * RMS, mean + nSig * RMS]
// with mean and RMS recomputed once inside +/-3 RMS (to reduce the effect of the
// tails), then limited to the histogram range. The same window is used for the
// reco, truth and reco/truth fits. The truth peaks are narrow, so the fixed
// default window [0.8, 1.2] would leave too few bins for the fit scan.
// Returns [0.8, 1.2] for an empty or degenerate histogram.
inline std::pair<double, double> DataDrivenSeed(TH1* h, double nSig = 2.0) {
    if (!h) return {0.8, 1.2};
    double m = h->GetMean(), s = h->GetRMS();
    if (s > 0) {
        const int first = h->GetXaxis()->GetFirst();
        const int last = h->GetXaxis()->GetLast();
        h->GetXaxis()->SetRangeUser(m - 3.0 * s, m + 3.0 * s);
        const double m2 = h->GetMean(), s2 = h->GetRMS();
        h->GetXaxis()->SetRange(first, last); // restore original range
        if (s2 > 0) {
            m = m2;
            s = s2;
        }
    }
    if (s <= 0) return {0.8, 1.2};
    double lo = m - nSig * s, hi = m + nSig * s;
    const double xlo = h->GetXaxis()->GetXmin();
    const double xhi = h->GetXaxis()->GetXmax();
    if (lo < xlo) lo = xlo;
    if (hi > xhi) hi = xhi;
    return {lo, hi};
}

// In-situ resolution of one sample. From the fitted widths of the data balance
// (sigma_RDB), the MC balance (sigma_RDB^MC), the truth balance (sigma_truth)
// and the MC reco/truth response (sigma_MC):
//   sigma_insitu    = sqrt(sigma_RDB^2    - sigma_truth^2)   data in-situ resolution
//   sigma_insitu^MC = sqrt(sigma_RDB^MC^2 - sigma_truth^2)   MC in-situ resolution
//   Delta_closure   = (sigma_MC - sigma_insitu^MC) / sigma_insitu
// with their errors (widths treated as independent). If a difference of squares
// is negative, the result is left at 0 and a warning is printed.
struct PerSampFit {
    double sig_RDB = 0, sig_RDB_err = 0, sig_truth = 0, sig_truth_err = 0;
    double sig_MC = 0, sig_MC_err = 0, sig_insitu = 0, sig_insitu_err = 0;
    double sig_MCi = 0, sig_MCi_err = 0, dC = 0, dC_err = 0;
    double sig_RDB_MC = 0, sig_RDB_MC_err = 0;
    bool ok = false;
};

inline PerSampFit ComputePerSampleJER(double sRDB_d, double sRDB_d_err, double sRDB_mc, double sRDB_mc_err,
                                      double sTr, double sTr_err, double sMC, double sMC_err) {
    PerSampFit p;
    p.sig_RDB = sRDB_d;
    p.sig_RDB_err = sRDB_d_err;
    p.sig_RDB_MC = sRDB_mc;
    p.sig_RDB_MC_err = sRDB_mc_err;
    p.sig_truth = sTr;
    p.sig_truth_err = sTr_err;
    p.sig_MC = sMC;
    p.sig_MC_err = sMC_err;
    if (sTr <= 0 || sMC <= 0) return p;
    double dd = sRDB_d * sRDB_d - sTr * sTr, dm = sRDB_mc * sRDB_mc - sTr * sTr;

    if (dd > 0) {
        p.sig_insitu = std::sqrt(dd);
    } else {
        Log::warning("FitResponse::ComputePerSampleJER",
                     Form("Negative variance for sig_insitu: sigma_data^2=%.5f sigma_truth^2=%.5f",
                          sRDB_d * sRDB_d, sTr * sTr));
    }
    if (dm > 0) {
        p.sig_MCi = std::sqrt(dm);
    } else {
        Log::warning("FitResponse::ComputePerSampleJER",
                     Form("Negative variance for sig^MC_insitu: sigma_MC^2=%.5f sigma_truth^2=%.5f",
                          sRDB_mc * sRDB_mc, sTr * sTr));
    }

    if (p.sig_insitu > 0)
        p.sig_insitu_err = std::sqrt(std::pow(sRDB_d * sRDB_d_err / p.sig_insitu, 2)
                                     + std::pow(sTr * sTr_err / p.sig_insitu, 2));
    if (p.sig_MCi > 0)
        p.sig_MCi_err = std::sqrt(std::pow(sRDB_mc * sRDB_mc_err / p.sig_MCi, 2)
                                  + std::pow(sTr * sTr_err / p.sig_MCi, 2));
    if (p.sig_insitu > 0) {
        p.dC = (p.sig_MC - p.sig_MCi) / p.sig_insitu;
        double num = p.sig_MC - p.sig_MCi;
        p.dC_err =
            std::sqrt(std::pow(p.sig_MC_err / p.sig_insitu, 2) + std::pow(p.sig_MCi_err / p.sig_insitu, 2)
                      + std::pow(num * p.sig_insitu_err / (p.sig_insitu * p.sig_insitu), 2));
    }
    p.ok = true;
    return p;
}

// One pass over the MC to fill the histograms used to choose the response
// binning (dynamic binning mode):
//   h_reco       reco balance R (always filled)
//   h_recoTruth  reco/truth response R / R_truth (only with need_truth_probe and
//                if the truth branch exists)
// The EASYBJETS_MASS_CUT selection is applied as in the fit loops, so the
// histograms describe the same events. The branch addresses are reset at the end.
// Returns whether the truth branch R_truth_<channel> was found; if not, the
// caller falls back to the "match" truth binning. warnPrefix is the tag of
// the warning printed when the branch is missing.
inline bool FillResponseBinningHists(TTree* tree_mc, const VariantConfig& cfg, bool need_truth_probe,
                                     TH1F* h_reco, TH1F* h_recoTruth, const char* warnPrefix) {
    std::vector<float>* t_m = nullptr;
    std::vector<float>* t_pt = nullptr;
    std::vector<float>* t_dR = nullptr;
    double t_w = 1.0, t_R = 0.0, t_R_truth = -1.0;
    int t_idx = -1;
    tree_mc->SetBranchAddress(cfg.branch_mass.c_str(), &t_m);
    tree_mc->SetBranchAddress(cfg.branch_pt.c_str(), &t_pt);
    tree_mc->SetBranchAddress(cfg.branch_R.c_str(), &t_R);
    tree_mc->SetBranchAddress("weight", &t_w);
    tree_mc->SetBranchAddress("selected_ljet_index", &t_idx);

    bool have_truth_branch = false, have_dR_truth = false;
    if (need_truth_probe) {
        std::string truth_R_br = "R_truth_" + cfg.channel;
        if (tree_mc->GetBranch(truth_R_br.c_str())) {
            tree_mc->SetBranchAddress(truth_R_br.c_str(), &t_R_truth);
            have_truth_branch = true;
        }
        if (tree_mc->GetBranch("matched_tjet_dR_sig")) {
            tree_mc->SetBranchAddress("matched_tjet_dR_sig", &t_dR);
            have_dR_truth = true;
        }
        if (!have_truth_branch)
            std::cerr << warnPrefix << " WARN: truth probe requested but " << truth_R_br
                      << " not found in MC tree; truth_bin will fall back to match.\n";
    }

    Long64_t n = tree_mc->GetEntries();
    for (Long64_t i = 0; i < n; ++i) {
        tree_mc->GetEntry(i);
        if (t_idx < 0 || !t_m || t_idx >= (int)t_m->size()) continue;
        if (MassCutConfig::Reject(t_m, t_idx)) continue; // EASYBJETS_MASS_CUT (no-op if off)
        h_reco->Fill(t_R, t_w);

        if (need_truth_probe && h_recoTruth && have_truth_branch && t_R_truth > 0) {
            double dR = -1.0;
            if (have_dR_truth && t_dR && t_idx < (int)t_dR->size()) dR = t_dR->at(t_idx);
            bool truthOK = (dR >= 0 && dR < cfg.maxTruthDR);
            if (!have_dR_truth) truthOK = true; // no dR branch: accept on R_truth > 0 alone
            if (truthOK) {
                double r_rt = t_R / t_R_truth; // pT_reco / pT_truth
                if (r_rt > 0 && r_rt < 5) h_recoTruth->Fill(r_rt, t_w);
            }
        }
    }
    tree_mc->ResetBranchAddresses();
    return have_truth_branch;
}

#endif // RESPONSE_FIT_CORE_H
