/**
 * Mass_Ratios_PeakWidth.cxx - DB (Z(bb)+gamma) entry point of the mass ratios 
 *
 * Reads the fit_summary.txt of every variant and computes the data/MC
 * ratios b-JMDP = mu_data/mu_MC and b-JMDR = sigma_data/sigma_MC, and their
 * double ratios variant/nominal, with three error treatments. The code is
 * shared with MJB and lives in common/macros/Mass_Ratios_PeakWidth.h.
 *
 * Usage (from ZbbGamma/macros, after Mass_Fit.cxx with the same EASYBJETS_RUN,
 * EASYBJETS_BINNING, EASYBJETS_FIT_MODEL and pT mode):
 *   export EASYBJETS_RUN=run3
 *   root 'Mass_Ratios_PeakWidth.cxx()'              // one pT bin (inclusive)
 *   root 'Mass_Ratios_PeakWidth.cxx("ptbins")'      // pT bins of the channel
 *   root 'Mass_Ratios_PeakWidth.cxx("ALL")'         // both
 *
 * Next step: Mass_Summary_Combined.cxx, once both channels are done.
 */
#include "../../common/macros/Mass_Ratios_PeakWidth.h"

void Mass_Ratios_PeakWidth(const std::string& ptMode = "inclusive") {
    Mass_Ratios_Width_Run(ptMode, "DB");
}
