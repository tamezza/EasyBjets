/**
 * Response_Fit_Inclusive.cxx - DB (Z(bb)+gamma) for the inclusive response fit
 *
 * Fits the balance R in data and MC in one inclusive bin, for each variant and
 * WP, and derives R_JES and R_JER. The code is shared with MJB and lives in
 * common/macros/Response_Fit_Inclusive.h (all options in its header).
 *
 * Usage (from ZbbGamma/macros):
 * 
 *   export EASYBJETS_RUN=run3
 *   export EASYBJETS_DATA_DIR=$HOME/lpnatlas_data/ZbbGamma
 *   export EASYBJETS_MASS_CUT=on        // set mass window, 50-150 GeV (default: no cut)
 *   root 'Response_Fit_Inclusive.cxx()'                  // nominal DB variant
 *   root 'Response_Fit_Inclusive.cxx("DB_bJER10v00")'   // one DB variant
 *   root 'Response_Fit_Inclusive.cxx("ALL")'             // all DB variants
 *
 * Next step: Response_Ratios_Inclusive.cxx with the same options.
 */
#include "../../common/macros/Response_Fit_Inclusive.h"

void Response_Fit_Inclusive(const std::string& variant = "DB") {
    FitResponse_Run(variant, "DB");
}
