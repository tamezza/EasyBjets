/**
 * Efficiency.cxx - MJB (Z(bb)+jets) - selection efficiency
 *
 * Efficiency after each cut for data, signal and background, per working point
 * and per sample, from the h_cutflow histograms of the slimmed ntuples. The code
 * is shared with DB and lives in common/macros/Efficiency.h.
 *
 * Usage (from ZbbJets/macros):
 * 
 *   export EASYBJETS_RUN=run3
 *   export EASYBJETS_DATA_DIR=$HOME/lpnatlas_data/ZbbJets
 *   root 'Efficiency.cxx()'                  // nominal MJB (default)
 *   root 'Efficiency.cxx("MJB_bJER10v00")'   // one MJB variant
 *   root 'Efficiency.cxx("ALL")'             // all MJB variants
 */
#include "../../common/macros/Efficiency.h"

void Efficiency(const std::string& variant = "MJB") {
    Efficiency_Run(variant, "MJB");
}
