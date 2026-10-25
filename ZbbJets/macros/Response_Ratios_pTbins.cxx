/**
 * Response_Ratios_pTbins.cxx - MJB (Z(bb)+jets) for the response ratios
 * in pT bins without mass cut 
 *
 * Same as Response_Ratios_pTbins_Mcuts.cxx, for the fits without mass cut
 * (Response_Fit_pTbins.cxx); output in the R2pTbin_MassInclusive folders. The
 * code is shared with DB and lives in common/macros/Response_Ratios_pTbins_Mcuts.h.
 *
 * Usage (from ZbbJets/macros, after Response_Fit_pTbins.cxx for all variants):
 * 
 *   export EASYBJETS_RUN=run3          // same as for the fit
 *   root 'Response_Ratios_pTbins.cxx()'                  // all corrections
 *   root 'Response_Ratios_pTbins.cxx("MJB_bJER10v00")'   // one correction
 *
 * Next step: Response_Summary_Combined_pTbins.cxx, once both channels are done.
 */
#include "../../common/macros/Response_Ratios_pTbins_Mcuts.h"

void Response_Ratios_pTbins(const std::string& variant = "ALL") {
    Ratio_Response_Run(variant, "MJB", "inclusive");
}
