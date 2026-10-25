/**
 * Response_Summary_Combined_Inclusive.cxx - DB + MJB summary of the inclusive
 * R_JES and R_JER
 *
 * One point per channel, read from the inclusive CSVs. Same file in both
 * channels; the plots go to the Run{2,3}Plots of the folder it is run from.
 * The code lives in common/macros/Response_Summary_Combined_Inclusive.h.
 *
 * Usage (from ZbbJets/macros, after Response_Ratios_Inclusive.cxx in both channels):
 * 
 *   export EASYBJETS_RUN=run3
 *   root 'Response_Summary_Combined_Inclusive.cxx()'                       // all corrections, DB + MJB
 *   root 'Response_Summary_Combined_Inclusive.cxx("bJER10v01,MuonCorr")'   // a subset
 *   root 'Response_Summary_Combined_Inclusive.cxx("ALL","DB")'             // one channel only
 * 
 *   arguments: corrections (no channel prefix), topology (combined | DB | MJB)
 */
#include "../../common/macros/Response_Summary_Combined_Inclusive.h"

void Response_Summary_Combined_Inclusive(const std::string& corrections = "ALL",
                                         const std::string& topology = "combined") {
    Ratio_Response_Combined_Summary_Inclusive_Run(corrections, topology);
}
