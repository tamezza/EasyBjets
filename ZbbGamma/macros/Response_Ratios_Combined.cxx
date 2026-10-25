/**
 * Response_Ratios_Combined.cxx - DB + MJB double ratios on one pT axis 
 *
 * Reads the pT x mass-cut CSVs of both channels and draws the double ratios
 * correction / nominal from 200 to 1200 GeV. Same file in both channels; the
 * plots go to the Run{2,3}Plots of the folder it is run from. The code lives in
 * common/macros/Response_Ratios_Combined.h.
 *
 * Usage (from ZbbGamma/macros, after Response_Ratios_pTbins_Mcuts.cxx in both channels):
 *   export EASYBJETS_RUN=run3
 *   root 'Response_Ratios_Combined.cxx()'              // every correction
 *   root 'Response_Ratios_Combined.cxx("bJER10v01")'   // one correction (no channel prefix)
 */
#include "../../common/macros/Response_Ratios_Combined.h"

void Response_Ratios_Combined(const std::string& variant = "ALL") {
    Ratio_Response_Combined_Run(variant);
}
