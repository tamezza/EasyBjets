/**
 * Response_Fit_pTbins.cxx - MJB (Z(bb)+jets) for the response fit in pT bins without mass cut
 *
 * Same as Response_Fit_pTbins_Mcuts.cxx, but with one mass bin ("MassInclusive")
 * instead of the three mass cuts. The code is shared with DB and lives in
 * common/macros/Response_Fit_pTbins_Mcuts.h (all options in its header).
 *
 * Usage (from ZbbJets/macros):
 * 
 *   export EASYBJETS_RUN=run3
 *   export EASYBJETS_DATA_DIR=$HOME/lpnatlas_data/ZbbJets
 *   export EASYBJETS_MASS_CUT=on        // set mass window, 50-150 GeV (default: no cut)
 *   root 'Response_Fit_pTbins.cxx()'                  // nominal MJB variant
 *   root 'Response_Fit_pTbins.cxx("MJB_bJER10v00")'   // one MJB variant
 *   root 'Response_Fit_pTbins.cxx("ALL")'             // all MJB variants
 *
 * Next step: Response_Ratios_pTbins.cxx with the same options.
 */
#include "../../common/macros/Response_Fit_pTbins_Mcuts.h"

void Response_Fit_pTbins(const std::string& variant = "MJB") {
    Rdb_pT_bins_MCut_Run(variant, "MJB", "inclusive");
}
