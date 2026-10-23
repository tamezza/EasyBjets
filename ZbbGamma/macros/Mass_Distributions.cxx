/**
 * Mass_Distributions.cxx - DB (Z(bb)+gamma)
 *
 * Data/MC stack plots of the jet mass, the asymmetry and the response, the
 * sideband scale factor of the background MC, and the background flavour
 * composition. The code is shared with MJB and lives in
 * common/macros/Mass_Distributions.h (all options in its header).
 *
 * Usage (from ZbbGamma/macros):
 *   export EASYBJETS_RUN=run3 // or run2
 *   export EASYBJETS_DATA_DIR=$HOME/lpnatlas_data/ZbbGamma
 *   export EASYBJETS_MASS_CUT=on        // 50-150 GeV, as used in this analysis
 *   root 'Mass_Distributions.cxx()'                   // nominal DB variant
 *   root 'Mass_Distributions.cxx("DB_bJER10v00")'     // one DB variant
 *   root 'Mass_Distributions.cxx("ALL")'              // all DB variants
 */
#include "../../common/macros/Mass_Distributions.h"

void Mass_Distributions(const std::string& variant = "DB") {
    Mass_Response_Run(variant, "DB");
}
