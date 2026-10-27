/**
 * Bootstrap_Correlations.cxx - DB (Z(bb)+gamma) - the bootstrap correlations
 *
 * Measures the correlations (rho) between each bJR10 correction and the nominal
 * for every quantity (b-JMDP, b-JMDR, b-JES, b-JER) and every binning
 * (inclusive, pT bins, pT x mass bins), in one run. They are used by the
 * "full" errors of the ratio macros. The code is shared with MJB and lives in
 * common/macros/Bootstrap_Correlations.h (method and outputs in its header).
 *
 * Usage (from ZbbGamma/macros):
 * 
 *   export EASYBJETS_RUN=run3
 *   export EASYBJETS_DATA_DIR=$HOME/lpnatlas_data/ZbbGamma
 *   export EASYBJETS_MASS_CUT=on               // same setting as the response fits
 *   export EASYBJETS_BOOT_N=200                 // number of replicas (default 200)
 *   root 'Bootstrap_Correlations.cxx()'                 // every DB bJER correction
 *   root 'Bootstrap_Correlations.cxx("DB_bJER10v01")'   // one correction
 *
 * Output: ../Run3Plots/Bootstrap/correlations.txt (../Run2Plots for run2).
 * Then run the ratio macros with EASYBJETS_CORR_METHOD=file and the same
 * EASYBJETS_RUN. DR uncertainties and plots: Bootstrap_Analyze.cxx().
 */
#include "../../common/macros/Bootstrap_Correlations.h"

void Bootstrap_Correlations(const std::string& variant = "ALL") {
    Bootstrap_Correlations_Run(variant, "DB");
}
