/**
 * CutFlow.cxx - MJB (Z(bb)+jets) - cutflow plots and tables
 *
 * Stacked cutflows, efficiencies, S/B and data/MC vs working point, from the
 * h_cutflow histograms of the slimmed ntuples. The code is shared with DB
 * and lives in common/macros/CutFlow.h (outputs and details in its header).
 *
 * Usage (from ZbbJets/macros):
 * 
 *   export EASYBJETS_RUN=run3
 *   export EASYBJETS_DATA_DIR=$HOME/lpnatlas_data/ZbbJets
 *   root 'CutFlow.cxx()'                  // nominal MJB (default)
 *   root 'CutFlow.cxx("MJB_bJER10v00")'   // one MJB variant
 *   root 'CutFlow.cxx("ALL")'             // all MJB variants
 */
#include "../../common/macros/CutFlow.h"

void CutFlow(const std::string& variant = "MJB") {
    CutFlow_Run(variant, "MJB");
}
