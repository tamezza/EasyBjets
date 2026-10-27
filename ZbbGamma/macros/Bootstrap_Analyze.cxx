/**
 * Bootstrap_Analyze.cxx - DB (Z(bb)+gamma) - the bootstrap analysis
 *
 * From the replicas written by Bootstrap_Correlations.cxx, computes the double
 * ratio of every observable (b-JMDP, b-JMDR, b-JES, sigma_insitu, b-JER) and its
 * uncertainty: bootstrap spread, and the full, corr-aware and uncorrelated
 * linear errors; plus the diagnostic plots. The code is shared with MJB and
 * lives in common/macros/Bootstrap_Analyze.h (outputs in its header).
 *
 * Usage (from ZbbGamma/macros, after Bootstrap_Correlations.cxx):
 * 
 *   export EASYBJETS_RUN=run3                      // same era as the bootstrap
 *   root 'Bootstrap_Analyze.cxx()'                 // all DB corrections
 *   root 'Bootstrap_Analyze.cxx("DB_bJER10v01")'   // one correction
 *
 * Output: ../Run3Plots/Bootstrap/analysis/ (../Run2Plots for run2).
 */
#include "../../common/macros/Bootstrap_Analyze.h"

void Bootstrap_Analyze(const std::string& variant = "ALL") {
    Bootstrap_Analyze_Run(variant, "DB");
}
