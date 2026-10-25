/**
 * Response_Fit_pTbins_Mcuts.cxx - MJB (Z(bb)+jets) for the response fit in pT bins x mass cuts (thin shim)
 *
 * Fits the balance R in each pT bin of the channel and in three jet-mass cuts
 * (50-80, 80-110, >= 110 GeV), for each variant and WP. The code is shared with
 * DB and lives in common/macros/Response_Fit_pTbins_Mcuts.h (all options in
 * its header).
 *
 * Usage (from ZbbJets/macros):
 * 
 *   export EASYBJETS_RUN=run3
 *   export EASYBJETS_DATA_DIR=$HOME/lpnatlas_data/ZbbJets
 *   root 'Response_Fit_pTbins_Mcuts.cxx()'                  // nominal MJB variant
 *   root 'Response_Fit_pTbins_Mcuts.cxx("MJB_bJER10v00")'   // one MJB variant
 *   root 'Response_Fit_pTbins_Mcuts.cxx("ALL")'             // all MJB variants
 *
 * Next step: Response_Ratios_pTbins_Mcuts.cxx with the same options.
 */
#include "../../common/macros/Response_Fit_pTbins_Mcuts.h"

void Response_Fit_pTbins_Mcuts(const std::string& variant = "MJB") {
    Rdb_pT_bins_MCut_Run(variant, "MJB");
}
