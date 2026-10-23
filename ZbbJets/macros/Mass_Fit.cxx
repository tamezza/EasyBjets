/**
 * Mass_Fit.cxx - MJB (Z(bb)+jets) entry point of the mass peak fit
 *
 * For each variant, WP and pT bin: chooses the background function on the
 * data sidebands (F-test + k-fold CV), subtracts it, fits the Z peak in data
 * and in the signal MC, and writes fit_summary.txt. The code is shared with
 * DB and lives in common/macros/Mass_Fit.h (all options in its header).
 *
 * Usage (from ZbbJets/macros):
 *   export EASYBJETS_RUN=run3
 *   export EASYBJETS_DATA_DIR=$HOME/lpnatlas_data/ZbbJets
 *   root 'Mass_Fit.cxx()'                                // all variants, one pT bin
 *   root 'Mass_Fit.cxx("ALL","ptbins")'                  // all variants, pT bins
 *   root 'Mass_Fit.cxx("MJB_bJER10v00","inclusive")'    // one variant
 *   arguments: variant ("ALL", "MJB", "MJB_bJER10v00", ...), pT mode (inclusive | ptbins)
 *
 * Next step: Mass_Ratios_PeakWidth.cxx with the same options and pT mode.
 */
#include "../../common/macros/Mass_Fit.h"

void Mass_Fit(const std::string& variant = "ALL", const std::string& ptMode = "inclusive") {
    Mass_Analysis_Run(variant, ptMode, "MJB");
}
