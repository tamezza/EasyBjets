/**
 * Mass_Summary_Combined.cxx - combined DB + MJB mass summary 
 *
 * Puts the b-JMDP and b-JMDR results of DB and MJB on one pT axis. Same file
 * in both channels: it reads the CSVs of DB and MJB, and the plots go to the
 * Run{2,3}Plots of the folder it is run from. The code lives in
 * common/macros/Mass_Summary_Combined.h.
 *
 * Usage (from ZbbGamma/macros, after Mass_Ratios_PeakWidth.cxx in both channels):
 *   export EASYBJETS_RUN=run3
 *   root 'Mass_Summary_Combined.cxx()'      // all corrections, DB + MJB, inclusive
 *   root 'Mass_Summary_Combined.cxx("ALL","combined","inclusive","ptbins")'
 *   root 'Mass_Summary_Combined.cxx("bJER10v01,MuonCorr","DB","inclusive")'
 *   arguments: corrections, topology (DB | MJB | combined), DB pT mode, MJB pT mode
 */
#include "../../common/macros/Mass_Summary_Combined.h"

void Mass_Summary_Combined(const std::string& corrections = "ALL",
                           const std::string& topology = "combined",
                           const std::string& dbPtMode = "inclusive",
                           const std::string& mjbPtMode = "inclusive") {
    Mass_Ratios_Combined_Summary_Run(corrections, topology, dbPtMode, mjbPtMode);
}
