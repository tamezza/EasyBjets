/**
 * Mass_Distributions_pTbins.cxx - DB (Z(bb)+gamma)
 *
 * pT-binned version of Mass_Distributions: the same stacked Data/MC plots
 * (mass, asymmetry, response, flavour split) in each pT bin of the channel,
 * with the sideband SF refitted per bin. The code is shared with MJB and
 * lives in common/macros/Mass_Distributions_pTbins.h (bins and options in
 * its header).
 *
 * Usage (from ZbbGamma/macros):
 *   export EASYBJETS_RUN=run3 // or run2
 *   export EASYBJETS_DATA_DIR=$HOME/lpnatlas_data/ZbbGamma
 *   export EASYBJETS_MASS_CUT=on        // 50-150 GeV, as used in this analysis
 *   root 'Mass_Distributions_pTbins.cxx()'                   // nominal DB variant
 *   root 'Mass_Distributions_pTbins.cxx("DB_bJER10v00")'     // one DB variant
 *   root 'Mass_Distributions_pTbins.cxx("ALL")'              // all DB variants
 */
#include "../../common/macros/Mass_Distributions_pTbins.h"

void Mass_Distributions_pTbins(const std::string& variant = "DB") {
    Mass_Response_pTbins_Run(variant, "DB");
}
