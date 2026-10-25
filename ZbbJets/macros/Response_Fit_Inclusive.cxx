/**
 * Response_Fit_Inclusive.cxx - MJB (Z(bb)+jets) for the inclusive response fit
 *
 * Fits the balance R in data and MC in one inclusive bin, for each variant and
 * WP, and derives R_JES and R_JER. The code is shared with DB and lives in
 * common/macros/Response_Fit_Inclusive.h (all options in its header).
 *
 * Usage (from ZbbJets/macros):
 * 
 *   export EASYBJETS_RUN=run3
 *   export EASYBJETS_DATA_DIR=$HOME/lpnatlas_data/ZbbJets
 *   export EASYBJETS_MASS_CUT=on        // optional mass window, 50-150 GeV (default: no cut)
 *   root 'Response_Fit_Inclusive.cxx()'                  // nominal MJB variant
 *   root 'Response_Fit_Inclusive.cxx("MJB_bJER10v00")'   // one MJB variant
 *   root 'Response_Fit_Inclusive.cxx("ALL")'             // all MJB variants
 *
 * Next step: Response_Ratios_Inclusive.cxx with the same options.
 */
#include "../../common/macros/Response_Fit_Inclusive.h"

void Response_Fit_Inclusive(const std::string& variant = "MJB") {
    FitResponse_Run(variant, "MJB");
}
