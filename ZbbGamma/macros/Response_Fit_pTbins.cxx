/**
 * Response_Fit_pTbins.cxx - DB (Z(bb)+gamma) for the response fit in pT bins without mass cut
 *
 * Same as Response_Fit_pTbins_Mcuts.cxx, but with one mass bin ("MassInclusive")
 * instead of the three mass cuts. The code is shared with MJB and lives in
 * common/macros/Response_Fit_pTbins_Mcuts.h (all options in its header).
 *
 * Usage (from ZbbGamma/macros):
 * 
 *   export EASYBJETS_RUN=run3
 *   export EASYBJETS_DATA_DIR=$HOME/lpnatlas_data/ZbbGamma
 *   export EASYBJETS_MASS_CUT=on        // optional mass window, 50-150 GeV (default: no cut)
 *   root 'Response_Fit_pTbins.cxx()'                  // nominal DB variant
 *   root 'Response_Fit_pTbins.cxx("DB_bJER10v00")'   // one DB variant
 *   root 'Response_Fit_pTbins.cxx("ALL")'             // all DB variants
 *
 * Next step: Response_Ratios_pTbins.cxx with the same options.
 */
#include "../../common/macros/Response_Fit_pTbins_Mcuts.h"

void Response_Fit_pTbins(const std::string& variant = "DB") {
    Rdb_pT_bins_MCut_Run(variant, "DB", "inclusive");
}
