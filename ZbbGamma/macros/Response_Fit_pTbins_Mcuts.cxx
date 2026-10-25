/**
 * Response_Fit_pTbins_Mcuts.cxx - DB (Z(bb)+gamma) for response fit in pT bins x mass cuts
 *
 * Fits the balance R in each pT bin of the channel and in three jet-mass cuts (50-80, 80-110, >= 110 GeV), 
 * for each variant and WP. The code is shared with MJB and lives in common/macros/Response_Fit_pTbins_Mcuts.h 
 * (all options in its header).
 *
 * Usage (from ZbbGamma/macros):
 * 
 *   export EASYBJETS_RUN=run3
 *   export EASYBJETS_DATA_DIR=$HOME/lpnatlas_data/ZbbGamma
 *   root 'Response_Fit_pTbins_Mcuts.cxx()'                  // nominal DB variant
 *   root 'Response_Fit_pTbins_Mcuts.cxx("DB_bJER10v00")'   // one DB variant
 *   root 'Response_Fit_pTbins_Mcuts.cxx("ALL")'             // all DB variants
 *
 * Next step: Response_Ratios_pTbins_Mcuts.cxx with the same options.
 */
#include "../../common/macros/Response_Fit_pTbins_Mcuts.h"

void Response_Fit_pTbins_Mcuts(const std::string& variant = "DB") {
    Rdb_pT_bins_MCut_Run(variant, "DB");
}
