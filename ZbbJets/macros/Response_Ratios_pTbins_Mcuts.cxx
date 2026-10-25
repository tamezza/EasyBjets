/**
 * Response_Ratios_pTbins_Mcuts.cxx - MJB (Z(bb)+jets) for the response ratios
 * in pT bins x mass cuts
 *
 * From the pT x mass-cut fits, computes R_JES, R_JER and the double ratios
 * correction / nominal, with their error treatments, and writes the CSV used
 * by the combined plots. The code is shared with DB and lives in
 * common/macros/Response_Ratios_pTbins_Mcuts.h.
 *
 * Usage (from ZbbJets/macros, after Response_Fit_pTbins_Mcuts.cxx for all variants):
 * 
 *   export EASYBJETS_RUN=run3          // same as for the fit
 *   root 'Response_Ratios_pTbins_Mcuts.cxx()'                  // all corrections
 *   root 'Response_Ratios_pTbins_Mcuts.cxx("MJB_bJER10v00")'   // one correction
 *
 * Next steps, once both channels are done: Response_Ratios_Combined.cxx and
 * Response_Summary_Combined.cxx.
 */
#include "../../common/macros/Response_Ratios_pTbins_Mcuts.h"

void Response_Ratios_pTbins_Mcuts(const std::string& variant = "ALL") {
    Ratio_Response_Run(variant, "MJB");
}
