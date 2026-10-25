/**
 * CutFlow.cxx - DB (Z(bb)+gamma) - cutflow plots and tables
 *
 * Stacked cutflows, efficiencies, S/B and data/MC vs working point, from the
 * h_cutflow histograms of the slimmed ntuples. The code is shared with MJB
 * and lives in common/macros/CutFlow.h (outputs and details in its header).
 *
 * Usage (from ZbbGamma/macros):
 * 
 *   export EASYBJETS_RUN=run3
 *   export EASYBJETS_DATA_DIR=$HOME/lpnatlas_data/ZbbGamma
 *   root 'CutFlow.cxx()'                 // nominal DB (default)
 *   root 'CutFlow.cxx("DB_bJER10v00")'   // one DB variant
 *   root 'CutFlow.cxx("ALL")'            // all DB variants
 */
#include "../../common/macros/CutFlow.h"

void CutFlow(const std::string& variant = "DB") {
    CutFlow_Run(variant, "DB");
}
