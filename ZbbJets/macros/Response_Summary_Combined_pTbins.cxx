/**
 * Response_Summary_Combined_pTbins.cxx - DB + MJB summary of R_JES and R_JER in
 * pT bins without mass cut
 *
 * One point per pT bin, read from the R2pTbin_MassInclusive CSVs. Same file in
 * both channels; the plots go to the Run{2,3}Plots of the folder it is run from.
 * The code lives in common/macros/Response_Summary_Combined_Inclusive.h.
 *
 * Usage (from ZbbJets/macros, after Response_Ratios_pTbins.cxx in both channels):
 * 
 *   export EASYBJETS_RUN=run3
 *   root 'Response_Summary_Combined_pTbins.cxx()'                       // all corrections, DB + MJB
 *   root 'Response_Summary_Combined_pTbins.cxx("bJER10v01,MuonCorr")'   // a subset
 *   root 'Response_Summary_Combined_pTbins.cxx("ALL","DB")'             // one channel only
 * 
 *   arguments: corrections (no channel prefix), topology (combined | DB | MJB)
 */
#include "../../common/macros/Response_Summary_Combined_Inclusive.h"

void Response_Summary_Combined_pTbins(const std::string& corrections = "ALL",
                                      const std::string& topology = "combined") {
    Ratio_Response_Combined_Summary_Inclusive_Run(corrections, topology, "ptbins");
}
