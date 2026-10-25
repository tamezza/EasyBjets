/**
 * Response_Ratios_Inclusive.cxx - MJB (Z(bb)+jets) for the inclusive response ratios
 *
 * From the inclusive fits, computes the data/MC ratios R_JES and R_JER and the
 * double ratios correction / nominal, with their error treatments, and writes
 * the inclusive CSV used by the combined summary. The code is shared with DB
 * and lives in common/macros/Response_Ratios_Inclusive.h.
 *
 * Usage (from ZbbJets/macros, after Response_Fit_Inclusive.cxx for all variants):
 * 
 *   export EASYBJETS_RUN=run3          // same as for the fit
 *   root 'Response_Ratios_Inclusive.cxx()'                  // all corrections
 *   root 'Response_Ratios_Inclusive.cxx("MJB_bJER10v00")'   // one correction
 *
 * Next step: Response_Summary_Combined_Inclusive.cxx, once both channels are done.
 */
#include "../../common/macros/Response_Ratios_Inclusive.h"

void Response_Ratios_Inclusive(const std::string& variant = "ALL") {
    Ratio_Response_Inclusive_Run(variant, "MJB");
}
