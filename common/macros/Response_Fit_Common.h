/**
 * ============================================================================
 *  Response_Fit_Common.h - types and pT bins shared by the response fits
 * ============================================================================
 * ===========================================
 * @author:    Thiziri AMEZZA
 * @email:     thiziri.amezza@cern.ch
 * @date:      10 April 2026
 * ===========================================
 *
 * Small header included by both response fits (and by Response_Fit_Core.h):
 * 
 *   Response_Fit_Inclusive.h      fit of R in one inclusive bin
 *   Response_Fit_pTbins_Mcuts.h   fit of R in pT bins x mass cuts, and in pT bins
 *                                 without mass cut (Response_Fit_pTbins.cxx)
 * 
 * It holds the fit result (RDBFitResults), the current channel, and the pT bins
 * of that channel, taken from BinConfig in CommonUtils.h:
 * 
 *   DB   200-300, 300-450 GeV              (BinConfig::ptBins_DB, 2 bins)
 *   MJB  450-600, 600-900, 900-1200 GeV    (BinConfig::ptBins_3, 3 bins)
 *
 * Two versions of the bin edges:
 * 
 *   CurrentPtBins()          for the selection: the last bin has no upper limit
 *   CurrentPtBinsDisplay()   for axes, labels and file names: the last edge is
 *                            the channel's upper edge (450 or 1200 GeV)
 *
 * Not run on its own and reads no option; how to run is in the two fit headers.
 *
 * ==========================================================================
 * @version:  1.0
 * ==========================================================================
 */

#ifndef RESPONSE_FIT_COMMON_H
#define RESPONSE_FIT_COMMON_H

#include "CommonUtils.h" 
#include <string>

// mean, sigma and chi2/ndf of one fit of an R (balance) distribution.
struct RDBFitResults {
    double mean = 0;
    double meanErr = 0;
    double sigma = 0;
    double sigmaErr = 0;
    double chi2ndf = 0;
};

// Largest number of pT bins of a channel (DB 2, MJB 3), used as array size.
constexpr int MAX_PT_BINS = 3;

// Current channel ("DB" or "MJB"), set once at the start of a run, so the
// helpers below do not need a channel argument. The pT bins come from BinConfig:
// BinConfig::ptBins_DB (2 bins) for DB, BinConfig::ptBins_3 (3 bins) for MJB.
inline std::string& CurrentChannel() {
    static std::string ch = "DB";
    return ch;
}
inline int CurrentNPtBins() {
    return (CurrentChannel() == "MJB") ? BinConfig::nPtBins_3 : BinConfig::nPtBins_DB;
}
inline const double* CurrentPtBins() {
    return (CurrentChannel() == "MJB") ? BinConfig::ptBins_3 : BinConfig::ptBins_DB;
}
// Bin edges for display (axes, bin centres, labels, file names): the same as
// CurrentPtBins(), but the last edge is 450 or 1200 GeV instead of "no limit".
// The selection uses CurrentPtBins(), so jets above that edge still go into
// the last bin.
inline const double* CurrentPtBinsDisplay() {
    return (CurrentChannel() == "MJB") ? BinConfig::ptBins_3_display : BinConfig::ptBins_DB_display;
}
inline const char* CurrentPtBinLabel(int i) {
    return (CurrentChannel() == "MJB") ? BinConfig::PtBinLabel_3(i) : BinConfig::PtBinLabel_DB(i);
}

#endif // RESPONSE_FIT_COMMON_H
