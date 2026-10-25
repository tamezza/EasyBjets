/**
 * Response_Summary_Combined.cxx - DB + MJB summary of R_JES and R_JER in pT bins
 * x mass cuts 
 *
 * Single ratio, double ratio and b_JES / b_JER response figures of both channels
 * on one pT axis. Same file in both channels; the plots go to the Run{2,3}Plots
 * of the folder it is run from. The code lives in
 * common/macros/Response_Summary_Combined.h.
 *
 * Usage (from ZbbGamma/macros, after Response_Ratios_pTbins_Mcuts.cxx in both channels):
 *   export EASYBJETS_RUN=run3
 *   root 'Response_Summary_Combined.cxx()'                       // all corrections, DB + MJB
 *   root 'Response_Summary_Combined.cxx("bJER10v01,MuonCorr")'   // a subset
 *   root 'Response_Summary_Combined.cxx("ALL","DB")'             // one channel only
 * 
 *   arguments: corrections (no channel prefix), topology (combined | DB | MJB)
 */
#include "../../common/macros/Response_Summary_Combined.h"

void Response_Summary_Combined(const std::string& corrections = "ALL",
                               const std::string& topology = "combined") {
    Ratio_Response_Combined_Summary_Run(corrections, topology);
}
