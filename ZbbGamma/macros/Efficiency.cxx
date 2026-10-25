/**
 * Efficiency.cxx - DB (Z(bb)+gamma) - selection efficiency
 *
 * Efficiency after each cut for data, signal and background, per working point
 * and per sample, from the h_cutflow histograms of the slimmed ntuples. The code
 * is shared with MJB and lives in common/macros/Efficiency.h.
 *
 * Usage (from ZbbGamma/macros):
 * 
 *   export EASYBJETS_RUN=run3
 *   export EASYBJETS_DATA_DIR=$HOME/lpnatlas_data/ZbbGamma
 *   root 'Efficiency.cxx()'                 // nominal DB (default)
 *   root 'Efficiency.cxx("DB_bJER10v00")'   // one DB variant
 *   root 'Efficiency.cxx("ALL")'            // all DB variants
 */
#include "../../common/macros/Efficiency.h"

void Efficiency(const std::string& variant = "DB") {
    Efficiency_Run(variant, "DB");
}
