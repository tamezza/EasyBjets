/**
 * ============================================================================
 *  Mass_Distributions_pTbins.h - Mass_Distributions in pT bins
 * ============================================================================
 * ===========================================
 * @author:    Thiziri AMEZZA
 * @email:     thiziri.amezza@cern.ch
 * @date:      18 April 2026
 * ===========================================
 *
 * The same plots as Mass_Distributions.h, once per pT bin of the channel. The
 * sideband SF is refitted in each bin and the bin is written on every plot.
 *
 * The pT bins are BinConfig in CommonUtils.h (not the PtBinConfig used by
 * Mass_Fit). The last bin has no upper limit, although its label shows the
 * channel's upper edge:
 *   DB   200-300, 300-450 GeV
 *   MJB  450-600, 600-900, 900-1200 GeV
 *
 * How to run (from ZbbGamma/macros for DB, from ZbbJets/macros for MJB):
 *   export EASYBJETS_RUN=run3
 *   export EASYBJETS_DATA_DIR=$HOME/lpnatlas_data/ZbbGamma    (ZbbJets for MJB)
 *   export EASYBJETS_MASS_CUT=on       for this analysis: jets with 50-150 GeV
 *   root 'Mass_Distributions_pTbins.cxx()'                nominal variant
 *   root 'Mass_Distributions_pTbins.cxx("DB_bJER10v00")'  one variant (MJB_... in ZbbJets)
 *   root 'Mass_Distributions_pTbins.cxx("ALL")'           all four variants
 *
 * Output: <PlotDir>/<variant>/Mass_Distributions_pTbins/<wp>/<ptbin>/
 * Options: the same as Mass_Distributions.h.
 *
 * ==========================================================================
 * @version:  1.0
 * ==========================================================================
 */
#ifndef MASS_DISTRIBUTIONS_PTBINS_H
#define MASS_DISTRIBUTIONS_PTBINS_H

#include "Mass_Distributions.h"

// Called by Mass_Distributions_pTbins.cxx in each channel's
// macros/ folder. "ALL" runs every variant of the channel.
inline void Mass_Response_pTbins_Run(const std::string& variant, const std::string& channel) {
    AtlasUtils::LoadAtlasStyle();
    MassCutConfig::PrintBanner("Mass_Response_pTbins");

    if (variant == "ALL") {
        for (const auto& key : GetPrimaryVariantKeys(channel)) {
            Log::info("Mass_Distributions_pTbins", "Running for: " + key);
            Mass_Response_pTbins_Run(key, channel);
        }
        return;
    }

    VariantConfig cfg = GetVariant(variant);

    const bool isMJB = (cfg.channel == "MJB");
    const int nPt = isMJB ? BinConfig::nPtBins_3 : BinConfig::nPtBins_DB;
    const double* edges = isMJB ? BinConfig::ptBins_3 : BinConfig::ptBins_DB; // last bin has no upper limit
    auto fileLab = [&](int i) {
        return isMJB ? BinConfig::PtBinFileLabel_3(i) : BinConfig::PtBinFileLabel_DB(i);
    };
    auto dispLab = [&](int i) {
        return isMJB ? BinConfig::PtBinLabel_3(i) : BinConfig::PtBinLabel_DB(i);
    };

    auto workingPoints = AnalysisWPConfig::StandardWPs();
    for (const auto& wp : workingPoints) {
        for (int b = 0; b < nPt; ++b) {
            Log::info("Mass_Distributions_pTbins",
                      Form("Processing WP %s, pT bin %s", wp.c_str(), fileLab(b)));
            std::string outDir = PathBuilder::OutputDirWP(cfg, "Mass_Distributions_pTbins", wp)
                                 + std::string(fileLab(b)) + "/";
            PathBuilder::MakeDir(outDir);
            std::vector<std::string> inputFiles = SampleConfig::AllSamplePaths(cfg, wp);

            HistogramManager histManager(wp, outDir, cfg);
            histManager.SetPtWindow(edges[b], edges[b + 1], dispLab(b));
            histManager.ProcessFiles(inputFiles);
            histManager.SaveHistograms();

            Log::info("Mass_Distributions_pTbins",
                      Form("Completed WP %s, pT bin %s", wp.c_str(), fileLab(b)));
        }
    }
}

#endif // MASS_DISTRIBUTIONS_PTBINS_H
