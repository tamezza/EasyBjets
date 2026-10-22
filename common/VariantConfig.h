/**
 * =============================================================================================
 *  VariantConfig.h - configuration shared by the DB and MJB analyses
 * =============================================================================================
 * ===========================================
 * @author:    Thiziri AMEZZA
 * @email:     thiziri.amezza@cern.ch
 * @date:      10 March 2026 
 * ===========================================
  *
 * One file holds everything the macros need to know about a "variant", i.e. a
 * jet calibration applied to one channel:
 *   - which ntuple branches to read (nominal vs regressed pT/mass)
 *   - where the slim ntuples live and where plots go
 *   - the mass window, sidebands and histogram binning used by the fits
 *
 * Contents:
 *   VariantConfig        the per-variant settings
 *   registry             the 8 variants (4 per channel), see detail::VariantConfigsStorage()
 *   PathBuilder          ntuple / plot path construction
 *   AtlasUtils           ATLAS style + label
 *   SampleConfig         sample indices, filenames, legend labels, colours
 *   DatasetConfig        subset of VariantConfig used by the mass fit
 *   PtBinConfig          pT bins per channel
 *   WPConfig             which working points to run
 *
 * To add a variant: add one entry to the registry below and, if needed, its
 * binning overrides in the block right after it. The macros pick it up through
 * GetVariant(key) and GetPrimaryVariantKeys(channel).
 * ==========================================================================
 * @version:  1.0
 * ==========================================================================
 */

#ifndef VARIANT_CONFIG_H
#define VARIANT_CONFIG_H

#include "RunEra.h"             // RunConfig::NtupleSubdir(), PlotDir(), LumiLabel()
#include "WorkingPointConfig.h" // WorkingPoint:: registry - WP labels / branch names

#include <string>
#include <map>
#include <vector>
#include <iostream>
#include <TSystem.h>
#include <stdexcept>
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <TROOT.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TH1F.h>
#include <THStack.h>
#include <TLegend.h>
#include <TGraphAsymmErrors.h>

// ============================================================================
// VARIANT CONFIGURATION STRUCTURE
// ============================================================================
struct VariantConfig {
    // --- Identity ---
    std::string key;            // "DB", "DB_bJER10v00", "MJB", "MJB_bJER10v00", ... 
    std::string label;          //  what goes in legends ("Nominal", "bJR10v00", "#mu-in-jet")
    std::string channel = "DB"; // "DB" or "MJB" 

    // slim ntuple location, see PathBuilder::SamplePath:
    //   <DataBase>/NtupleSlim_Run{2,3}/<channel>_<path_tag>_<wp>wp_slim/<sample>_<file_tag>_<wp>wp_slim.root
    std::string path_tag;       // "nominal", "bJER10v00", "bJER10v01", "MuonCorr"
    std::string file_tag;       // "nominal", "bJER10v00", "bJER10v01", "MuonCorr" 

    // Branch names - eta/phi/e are the same across all variants, only pT and mass change.
    std::string branch_pt;      
    std::string branch_mass;    
    std::string branch_R = "R_DB";  // Response branch — "R_DB" for DB, "R_MJB" for MJB.

    // --- Tree name: Slimmer writes "nominal" for the pre-mass-cut tree ---
    std::string tree_name = "nominal";

    // --- Output directory ---
    std::string plot_subdir;    

    // --- Response histogram binning - Default ---
    // Regressed variants need finer binning for the Gaussian fit to resolve the peak properly.
    int nBinsR     = 30;        
    int nBinsTruth = 70;        // truth/ref and reco/truth histograms

    // Optional overrides. An empty vector means "use the scalar above".
    //   *_perMassCut : index 0 = 50-80, 1 = 80-110, 2 = >=110 GeV
    //   *_perPtBin   : index = pT bin of the channel (PtBinConfig), wins over perMassCut

    std::vector<int> nBinsR_perMassCut;    
    std::vector<int> nBinsTruth_perMassCut;
    std::vector<int> nBinsR_perPtBin;       
    std::vector<int> nBinsTruth_perPtBin;   
 
    int GetNBinsR(int massCutIdx) const {
        if (massCutIdx >= 0 && massCutIdx < (int)nBinsR_perMassCut.size())
            return nBinsR_perMassCut[massCutIdx];
        return nBinsR;
    }
    int GetNBinsTruth(int massCutIdx) const {
        if (massCutIdx >= 0 && massCutIdx < (int)nBinsTruth_perMassCut.size())
            return nBinsTruth_perMassCut[massCutIdx];
        return nBinsTruth;
    }
    int GetNBinsR_pt(int ptBinIdx, int massCutIdx) const {
        if (ptBinIdx >= 0 && ptBinIdx < (int)nBinsR_perPtBin.size())
            return nBinsR_perPtBin[ptBinIdx];
        return GetNBinsR(massCutIdx);
    }
    int GetNBinsTruth_pt(int ptBinIdx, int massCutIdx) const {
        if (ptBinIdx >= 0 && ptBinIdx < (int)nBinsTruth_perPtBin.size())
            return nBinsTruth_perPtBin[ptBinIdx];
        return GetNBinsTruth(massCutIdx);
    }

    // --- Mass analysis parameters (per-channel defaults; override per-variant in the registry) ---
    double mass_xmin      = 50;
    double mass_xmax      = 150;    // 200 in nominal for DB, 150 for MJB 
    int    sideband2_max  = 150;    // upper edge of the high-mass sideband
    int    min_bins       = 10;
    int    max_bins       = 20;     
    double sideband1_low  = 50;
    double sideband1_high = 75;     
    double sideband2_low  = 110;

    // --- OPTIONAL: per-pT-bin Mass-Analysis binning overrides of min_bins / max_bins ---
    std::vector<int> min_bins_perPtBin;   
    std::vector<int> max_bins_perPtBin;   

    int GetMinBins(int ptBinIdx) const {
        if (ptBinIdx >= 0 && ptBinIdx < (int)min_bins_perPtBin.size())
            return min_bins_perPtBin[ptBinIdx];
        return min_bins;
    }
    int GetMaxBins(int ptBinIdx) const {
        if (ptBinIdx >= 0 && ptBinIdx < (int)max_bins_perPtBin.size())
            return max_bins_perPtBin[ptBinIdx];
        return max_bins;
    }

    // Max dR(reco, truth) for an event to enter the JER closure histograms.
    double maxTruthDR = 0.75;

    // seed window of the iterative Gaussian fit on R
    double fit_seed_lo = 0.8;
    double fit_seed_hi = 1.2;

    // Positional constructor, 18 arguments. The order here must match the
    // registry entries below, so be careful when adding a field.
    VariantConfig(
        const std::string& key_            = "",
        const std::string& label_          = "",
        const std::string& channel_        = "DB",           
        const std::string& path_tag_       = "",
        const std::string& file_tag_       = "",
        const std::string& branch_pt_      = "ljet_pt",
        const std::string& branch_mass_    = "ljet_m",
        const std::string& plot_subdir_    = "DB",
        double xmin  = 50,  double xmax  = 200, int sb2max = 200,
        int    minb  = 10,  int    maxb  = 20,
        double sb1lo = 50,  double sb1hi = 75,  double sb2lo = 110,
        const std::string& tree_name_     = "nominal",
        int nBinsR_     = 30,
        int nBinsTruth_ = 70,
        const std::string& branch_R_      = "R_DB",
        double maxTruthDR_  = 0.75,
        double fit_seed_lo_ = 0.8,
        double fit_seed_hi_ = 1.2
    ) : key(key_), label(label_), channel(channel_),
        path_tag(path_tag_), file_tag(file_tag_),
        branch_pt(branch_pt_), branch_mass(branch_mass_),
        branch_R(branch_R_),
        tree_name(tree_name_),
        plot_subdir(plot_subdir_),
        nBinsR(nBinsR_),
        nBinsTruth(nBinsTruth_),
        mass_xmin(xmin), mass_xmax(xmax), sideband2_max(sb2max),
        min_bins(minb), max_bins(maxb),
        sideband1_low(sb1lo), sideband1_high(sb1hi), sideband2_low(sb2lo),
        maxTruthDR(maxTruthDR_),
        fit_seed_lo(fit_seed_lo_), fit_seed_hi(fit_seed_hi_) {}
        

    // axis labels that depend on the channel. Use these instead of writing
    // "R_{DB}" in the macros, otherwise MJB plots end up with the wrong label.
    inline std::string RAxisLabel() const {
        return (channel == "MJB") ? std::string("R_{MJB}") : std::string("R_{DB}");
    }
    inline std::string RMeanLabel() const {
        return std::string("#LT ") + RAxisLabel() + std::string("#GT");
    }
    inline std::string RSigmaLabel() const {
        return std::string("#sigma (") + RAxisLabel() + std::string(")");
    }
    inline std::string RFileTag() const {
        return (channel == "MJB") ? std::string("R_MJB") : std::string("R_DB");
    }
};


// ============================================================================
// VARIANT REGISTRY — ADD NEW VARIANTS HERE
// ============================================================================
// Sidebands and bin counts below were tuned on the Run samples. 
// Run 2 and Run 3 needed slightly different values
//
// Read access goes through GetVariantConfigs() / GetVariant(). The mutable
// accessor exists for binning scans only.
// ============================================================================
namespace detail {

    inline std::map<std::string, VariantConfig>& VariantConfigsStorage() {
        static std::map<std::string, VariantConfig> configs = {

            // ==================== DB : Z(->bb)+gamma - Direct Balance ====================

            { "DB", VariantConfig(
                "DB",                       // key
                "Nominal",                  // label
                "DB",                       // channel
                "nominal",                  // path_tag 
                "nominal",                  // file_tag  
                "ljet_pt",                  // branch_pt
                "ljet_m",                   // branch_mass
                "DB",                       // plot_subdir
                50, 160, 160,               // mass_xmin, mass_xmax, sideband2_max
                10, 15,                     // min_bins, max_bins
                50, 75, 110,                // sideband1_low, sideband1_high, sideband2_low 
                "nominal",                  // tree_name
                40,                         // nBinsR 
                30,                         // nBinsTruth
                "R_DB"                      // R_branch
            ) },

            { "DB_bJER10v00", VariantConfig(
                "DB_bJER10v00", "bJR10v00", "DB",
                "bJER10v00", "bJER10v00",
                "ljet_bJR10v00_pt", "ljet_bJR10v00_mass",
                "DB_bJER10v00",
                50, 150, 150,               
                10, 25,                     
                50, 75, 110,
                "nominal",
                30,                         
                30,                         
                "R_DB"
            ) },

            { "DB_bJER10v01", VariantConfig(
                "DB_bJER10v01", "bJR10v01", "DB",
                "bJER10v01", "bJER10v01",
                "ljet_bJR10v01_pt", "ljet_bJR10v01_mass",
                "DB_bJER10v01",
                50, 150, 150,               
                10, 30,                     
                50, 80, 110,
                "nominal",
                40,                         
                30,                         
                "R_DB"
            ) },

            { "DB_MuonCorr", VariantConfig(
                "DB_MuonCorr", "#mu-in-jet", "DB",
                "MuonCorr", "MuonCorr",
                "ljet_pt", "ljet_m",
                "DB_MuonCorr",
                50, 160, 160,               
                10, 15, 
                50, 75, 110,
                "nominal",
                40,                         
                30,                         
                "R_DB"
            ) },

            // ==================== MJB : Z(->bb)+jets - Multijet Balance) ====================

            { "MJB", VariantConfig(
                "MJB", "Nominal", "MJB",
                "nominal", "nominal",
                "ljet_pt", "ljet_m",
                "MJB",
                50, 150, 150,
                10, 50, 
                50, 65, 110,
                "nominal",
                50, 90, 
                "R_MJB"
            ) },

            { "MJB_bJER10v00", VariantConfig(
                "MJB_bJER10v00", "bJR10v00", "MJB",
                "bJER10v00", "bJER10v00",
                "ljet_bJR10v00_pt", "ljet_bJR10v00_mass",
                "MJB_bJER10v00",
                50, 150, 150,
                0, 50, 
                50, 65, 110, 
                "nominal",
                50, 90, 
                "R_MJB"
            ) },

            // The multijet SF window for this variant (bJER10v01 regression model) starts 
            // at 60 GeV instead of 50, hence mass_xmin = 60 and sideband1 = [60, 75].
            { "MJB_bJER10v01", VariantConfig(
                "MJB_bJER10v01", "bJR10v01", "MJB",
                "bJER10v01", "bJER10v01",
                "ljet_bJR10v01_pt", "ljet_bJR10v01_mass",
                "MJB_bJER10v01",
                60, 150, 150,
                10, 45,
                60, 75, 110,
                "nominal",
                50, 90, 
                "R_MJB"
            ) },

            { "MJB_MuonCorr", VariantConfig(
                "MJB_MuonCorr", "#mu-in-jet", "MJB",
                "MuonCorr", "MuonCorr",
                "ljet_pt", "ljet_m",
                "MJB_MuonCorr",
                50, 150, 150,
                10, 50,  
                50, 65, 110,
                "nominal",
                50, 90,
                "R_MJB"
            ) },
        };

        // Binning overrides. Applied once; without the flag every call would
        // reset vectors that a scan may have modified at runtime.
        static bool _overrides_applied = false;
        if (!_overrides_applied) {
            _overrides_applied = true;

            // per mass cut: {50-80, 80-110, >=110}
            configs["DB"].nBinsR_perMassCut      = {25, 30, 15}; 
            configs["DB"].nBinsTruth_perMassCut  = {30, 55, 20}; 
            configs["DB_bJER10v00"].nBinsR_perMassCut      = {25, 25, 20};
            configs["DB_bJER10v00"].nBinsTruth_perMassCut  = {30, 55, 30};
            configs["DB_bJER10v01"].nBinsR_perMassCut      = { 15, 15, 20}; 
            configs["DB_bJER10v01"].nBinsTruth_perMassCut  = {30, 30, 55};  
            configs["DB_MuonCorr"].nBinsR_perMassCut       = {30, 35, 15}; 
            configs["DB_MuonCorr"].nBinsTruth_perMassCut   = {30, 35, 20}; 

            // MJB keeps the scalar defaults except for bJER10v01, whose peak is much narrower
            configs["MJB_bJER10v01"].nBinsR_perMassCut     = {60, 75, 30}; //{20, 20, 15};
            configs["MJB_bJER10v01"].nBinsTruth_perMassCut = {70, 95, 50}; //{40, 40, 25};

            // per pT bin, mass fit. The last MJB bin has the lowest statistics
            // and does better with a coarser binning.
            configs["MJB"].max_bins_perPtBin           = {50, 20};
            configs["MJB_MuonCorr"].max_bins_perPtBin           = {50, 20};

            // per pT bin, response fit (one count per pT region, all mass cuts)
            configs["DB"].nBinsR_perPtBin       = {40, 45};
            configs["DB"].nBinsTruth_perPtBin       = {30, 30};
            configs["DB_bJER10v00"].nBinsR_perPtBin       = {30, 25}; 
            configs["DB_bJER10v00"].nBinsTruth_perPtBin       = {30, 30};
            configs["DB_bJER10v01"].nBinsR_perPtBin       = {20, 20}; 
            configs["DB_bJER10v01"].nBinsTruth_perPtBin       = {68, 65};
            configs["DB_MuonCorr"].nBinsR_perPtBin       = {40, 20};
            configs["DB_MuonCorr"].nBinsTruth_perPtBin       = {60, 30};

            configs["MJB"].nBinsR_perPtBin       = {60, 40, 45};
            configs["MJB"].nBinsTruth_perPtBin       = {100, 100, 100};
            configs["MJB_bJER10v00"].nBinsR_perPtBin       = {40, 45, 35}; 
            configs["MJB_bJER10v00"].nBinsTruth_perPtBin       = {90, 90, 90};
            configs["MJB_bJER10v01"].nBinsR_perPtBin       = {75, 60, 35}; 
            configs["MJB_bJER10v01"].nBinsTruth_perPtBin       = {95, 95, 95};
            configs["MJB_MuonCorr"].nBinsR_perPtBin       = {60, 40, 45};
            configs["MJB_MuonCorr"].nBinsTruth_perPtBin       = {100, 100, 100};

        }

        // ---- Run sideband overrides ----
        // Global sideband / mass-window override from the environment, applied
        // to every variant, e.g.; 
        //   export EASYBJETS_SIDEBAND1_HIGH=80
        //   export EASYBJETS_SIDEBAND2_LOW=115

        static bool _era_overrides_applied = false;
        if (!_era_overrides_applied) {
            _era_overrides_applied = true;

            //if (RunConfig::IsRun3()) {
            //}

            auto get_env_double = [](const char* name) -> std::pair<bool, double> {
                if (const char* env = std::getenv(name)) {
                    if (env[0] != '\0') {
                        try { return {true, std::stod(env)}; } catch (...) {}
                    }
                }
                return {false, 0.0};
            };
            auto [has_sb1lo, env_sb1lo]  = get_env_double("EASYBJETS_SIDEBAND1_LOW");
            auto [has_sb1hi, env_sb1hi]  = get_env_double("EASYBJETS_SIDEBAND1_HIGH");
            auto [has_sb2lo, env_sb2lo]  = get_env_double("EASYBJETS_SIDEBAND2_LOW");
            auto [has_sb2max, env_sb2max] = get_env_double("EASYBJETS_SIDEBAND2_MAX");
            auto [has_xmin, env_xmin]    = get_env_double("EASYBJETS_MASS_XMIN");
            auto [has_xmax, env_xmax]    = get_env_double("EASYBJETS_MASS_XMAX");

            if (has_sb1lo || has_sb1hi || has_sb2lo || has_sb2max || has_xmin || has_xmax) {
                std::cout << "[VariantConfig] Applying env-var sideband overrides:";
                if (has_sb1lo)  std::cout << " sb1_lo="  << env_sb1lo;
                if (has_sb1hi)  std::cout << " sb1_hi="  << env_sb1hi;
                if (has_sb2lo)  std::cout << " sb2_lo="  << env_sb2lo;
                if (has_sb2max) std::cout << " sb2_max=" << env_sb2max;
                if (has_xmin)   std::cout << " xmin="    << env_xmin;
                if (has_xmax)   std::cout << " xmax="    << env_xmax;
                std::cout << "\n";
                for (auto& [k, v] : configs) {
                    if (has_sb1lo)  v.sideband1_low  = env_sb1lo;
                    if (has_sb1hi)  v.sideband1_high = env_sb1hi;
                    if (has_sb2lo)  v.sideband2_low  = env_sb2lo;
                    if (has_sb2max) v.sideband2_max  = (int)env_sb2max;
                    if (has_xmin)   v.mass_xmin      = env_xmin;
                    if (has_xmax)   v.mass_xmax      = env_xmax;
                }
            }
        }

        return configs;
    }

}  // namespace detail


inline const std::map<std::string, VariantConfig>& GetVariantConfigs() {
    return detail::VariantConfigsStorage();
}

// Get variant keys for a given channel, in the order they are usually plotted.
// An unknown channel returns every key.
inline std::vector<std::string> GetPrimaryVariantKeys(const std::string& channel = "DB") {
    if (channel == "DB")  return {"DB", "DB_bJER10v00", "DB_bJER10v01", "DB_MuonCorr"};
    if (channel == "MJB") return {"MJB", "MJB_bJER10v00", "MJB_bJER10v01", "MJB_MuonCorr"};
    std::vector<std::string> all;
    for (const auto& [k, v] : GetVariantConfigs()) all.push_back(k);
    return all;
}

inline const VariantConfig& GetVariant(const std::string& key) {
    const auto& configs = GetVariantConfigs();
    auto it = configs.find(key);
    if (it == configs.end()) {
        std::cerr << "ERROR: Unknown variant '" << key << "'. Available variants:\n";
        for (const auto& [k, v] : configs)
            std::cerr << "  - " << k << " (" << v.label << ")\n";
        throw std::runtime_error("Unknown variant: " + key);
    }
    return it->second;
}

// Normal analysis code should use the const version GetVariantConfigs().
//For binning scans. Changes persist for the rest of the ROOT session, so
// prefer copying one VariantConfig and editing the copy when possible.
inline std::map<std::string, VariantConfig>& GetVariantConfigsMutable() {
    return detail::VariantConfigsStorage();
}


// ============================================================================
// PATH BUILDER UTILITIES
// ============================================================================
namespace PathBuilder {

    inline std::string GetEosPath() {
        // $EASYBJETS_DATA_DIR if set (e.g. an eos path), otherwise fall back to "../"
        // for local development, e.g.
        //   export EASYBJETS_DATA_DIR=root://eosuser.cern.ch//eos/user/t/tamezza/Documents/EasyBjets/Zbby/
        if (const char* env = std::getenv("EASYBJETS_DATA_DIR")) {
            std::string p(env);
            if (!p.empty() && p.back() != '/') p.push_back('/');
            return p;
        }
        return "../";
    }

    // "0p3" / "30" / "0p30" -> "0p30". Unknown WPs are guessed from the digits.
    inline std::string WPtoFileWP(const std::string& wp) {
        if (WorkingPoint::IsValid(wp)) return WorkingPoint::FolderLabel(wp);
        if (wp.find('p') != std::string::npos) return wp;
        if (wp == "125") return "1p25";
        if (wp == "155") return "1p55";
        return "0p" + wp;
    }

    // Display label for a working point, e.g. "0p30" -> "0.30 % QCD Eff. WP".
    inline std::string WPtoLabel(const std::string& wp) {
        if (WorkingPoint::IsValid(wp)) return WorkingPoint::DisplayLabel(wp);
        if (wp.find('p') != std::string::npos) {
            std::string tmp = wp;
            std::replace(tmp.begin(), tmp.end(), 'p', '.');
            return tmp + " % QCD Eff. WP";
        }
        double val = std::stod(wp) / 100.0;
        char buf[64];
        snprintf(buf, sizeof(buf), "%.2f%% QCD Eff. WP", val);
        return std::string(buf);
    }

    // Build ntuple directory name. Prefix is driven by cfg.channel, e.g.
    //   DB nominal:   "DB_nominal_0p46wp_slim"
    //   MJB bJER10v01:"MJB_bJER10v01_0p46wp_slim"
    inline std::string NtupleDir(const VariantConfig& cfg, const std::string& fileWP) {
        return cfg.channel + "_" + cfg.path_tag + "_" + fileWP + "wp_slim";
    }

    // Build filename for a given sample, e.g.
    //   DB nominal:   "data_nominal_0p46wp_slim.root"
    //   MJB bJER10v01:"data_bJER10v01_0p46wp_slim.root"
    inline std::string NtupleFile(const VariantConfig& /*cfg*/,
                                   const std::string& sample,
                                   const std::string& /*fileWP*/,
                                   const std::string& file_tag,
                                   const std::string& fileWP_for_name) {
        return sample + "_" + file_tag + "_" + fileWP_for_name + "wp_slim.root";
    }

    // Overload to match the old DB/MJB API signature.
    inline std::string NtupleFile(const VariantConfig& cfg,
                                   const std::string& sample,
                                   const std::string& fileWP) {
        return sample + "_" + cfg.file_tag + "_" + fileWP + "wp_slim.root";
    }

    // Topology directory under DataSets/, derived from the channel:
    //   DB  -> ZbbGamma,  MJB -> ZbbJets.
    inline std::string TopologyDir(const VariantConfig& cfg) {
        return (cfg.channel == "MJB") ? "ZbbJets" : "ZbbGamma";
    }

    // Base directory that contains the NtupleSlim_{Run2|Run3} subdirectories.
    // $EASYBJETS_DATA_DIR set it is taken as the per-topology directory;
    // otherwise we assume the macro runs from <topology>/macros/.
    inline std::string DataBase(const VariantConfig& cfg) {
        if (const char* env = std::getenv("EASYBJETS_DATA_DIR")) {
            std::string p(env);
            if (!p.empty() && p.back() != '/') p.push_back('/');
            return p;
        }
        return std::string("../../DataSets/") + TopologyDir(cfg) + "/";
    }

    inline std::string SamplePath(const VariantConfig& cfg,
                                   const std::string& sample,
                                   const std::string& wp) {
        std::string fileWP = WPtoFileWP(wp);
        return DataBase(cfg) + RunConfig::NtupleSubdir() + "/"
             + NtupleDir(cfg, fileWP) + "/"
             + NtupleFile(cfg, sample, fileWP);
    }

    inline std::string DataPath(const VariantConfig& cfg, const std::string& wp)
        { return SamplePath(cfg, "data", wp); }

    // NOTE: SignalPath and BkgPath now read their filenames from SampleConfig
    // using cfg.channel. They are defined AFTER SampleConfig below to avoid
    // a forward-declaration tangle.

    // <PlotDir>/<variant>/<analysis>/ and <PlotDir>/<variant>/<analysis>/<wp>/
    inline std::string OutputDir(const VariantConfig& cfg, const std::string& analysis) {
        return (RunConfig::PlotDir() + "/") + cfg.plot_subdir + "/" + analysis + "/";
    }

    inline std::string OutputDirWP(const VariantConfig& cfg,
                                    const std::string& analysis,
                                    const std::string& wp) {
        return (RunConfig::PlotDir() + "/") + cfg.plot_subdir + "/" + analysis + "/"
             + WPtoFileWP(wp) + "/";
    }

    inline void MakeDir(const std::string& path) {
        gSystem->Exec(("mkdir -p " + path).c_str());
    }

}  // namespace PathBuilder


// ============================================================================
// ATLAS STYLE UTILITIES
// ============================================================================
namespace AtlasUtils {

    inline void LoadAtlasStyle() {
        static bool loaded = false;
        if (loaded) return;
        if (gROOT->GetGlobalFunction("SetAtlasStyle")) {
            gROOT->ProcessLine("SetAtlasStyle()");
        } else {
            int ret = gROOT->LoadMacro("../../common/Utils/AtlasStyle.C");
            if (ret == 0) {
                gROOT->ProcessLine("SetAtlasStyle()");
            } else {
                std::cerr << "Warning: AtlasStyle.C not found. Using default style.\n";
                gStyle->SetOptStat(0);
                gStyle->SetPadTickX(1);
                gStyle->SetPadTickY(1);
                gStyle->SetLegendBorderSize(0);
            }
        }
        loaded = true;
    }

    inline void DrawAtlasLabel(double x1 = 0.15,  double y1 = 0.85,
                                double x2_offset = 0.1, double y2_offset = -0.04) {
        TLatex latex;
        latex.SetNDC();
        latex.SetTextFont(72);
        latex.SetTextColor(kBlack);
        latex.SetTextSize(0.04);
        latex.DrawLatex(x1, y1, "ATLAS");
        latex.SetTextFont(42);
        latex.DrawLatex(x1 + x2_offset, y1, "Internal");
        latex.DrawLatex(x1, y1 + y2_offset, RunConfig::LumiLabel().c_str());
    }

}  // namespace AtlasUtils


// ============================================================================
// WORKING POINT CONFIGURATION
// ============================================================================
namespace WPConfig {

    inline std::vector<std::string> DefaultWPs() {
        return {"46"};
        //return {"25", "30", "37", "46", "58", "74", "94", "125", "155"};
    }

    inline std::vector<std::string> AllWPs() {
        return {"25", "30", "37", "46", "58", "74", "94", "125", "155"};
    }

    inline std::vector<int> DefaultColors() {
        return {kRed, kMagenta, kOrange+7};
    }
    inline std::vector<int> AllColors() {
        return {kRed, kBlue, kGreen+2, kMagenta, kCyan+2,
                kOrange+7, kViolet+3, kPink+2, kSpring+4};
    }

}  // namespace WPConfig


// ============================================================================
// SAMPLE CONFIGURATION
// ============================================================================
// Three logical samples for both channels:
//   kData      = real data
//   kMultijets = main non-resonant background
//                  DB:  gamma+jets        (file: gamma_jets_sherpa)
//                  MJB: multijets         (file: multijets_pythia)
//   kZbb       = Z(->bb)+X MC signal
//                  DB:  Z(->bb)+gamma     (file: Zbby)
//                  MJB: Z(->bb)+jets      (file: ZbbJets)
//
// Sample index convention (channel-independent)
// --------------------------------------------- 
// Sample indices are the same in both channels so the stack / cutflow helpers
// do not need to know which one they are in, 
// The core three samples are ALWAYS present for every variant:
//
// ============================================================================
namespace SampleConfig {

    // ------------------------------------------------------------------------
    // Sample index convention (channel-independent)
    // ------------------------------------------------------------------------
    // The core three samples are ALWAYS present for every variant:
    //   kData       = 0  (data)
    //   kMultijets  = 1  main MC background - gets the sideband scale factor
    //                    (gamma+jets for DB, multijets for MJB; the name kept as
    //                    "kMultijets" for channel-agnostic code)
    //   kZbb        = 2  Z(->bb) MC signal 
    //
    // Optional: the minor MC backgrounds live has indices [3..kNSamples]. 
    // (ttbar, W(qq), Z(qq)) are optional and off by default.
    //
    // TWO INDEPENDENT FLAGS — one per channel:
    //   EASYBJETS_DB_MINOR_BKGS  -> DB:  tty / Wqqy / Zqqy
    //   EASYBJETS_MJB_MINOR_BKGS -> MJB: top / WqqJets / ZqqJets
    //
    // To enable them per channel with env vars:
    //   root EASYBJETS_DB_MINOR_BKGS=1 ...
    //   root EASYBJETS_MJB_MINOR_BKGS=1 ...
    // Or edit the defaults below.
    //
    // When only one channel enables them, the other gets "__unused_*" filenames at
    // indices 3..5 so the vectors keep the same length; the macros skip files that
    // do not exist.
    // ------------------------------------------------------------------------
    #ifndef EASYBJETS_DB_MINOR_BKGS
    #define EASYBJETS_DB_MINOR_BKGS  0   // set to 1 to enable DB: tty/Wqqy/Zqqy
    #endif
    #ifndef EASYBJETS_MJB_MINOR_BKGS
    #define EASYBJETS_MJB_MINOR_BKGS 0   // set to 1 to enable MJB: top/WqqJets/ZqqJets
    #endif
    #define EASYBJETS_ANY_MINOR_BKGS (EASYBJETS_DB_MINOR_BKGS || EASYBJETS_MJB_MINOR_BKGS)

    enum SampleIdx {
        kData       = 0,
        kMultijets  = 1,  // main background (gets sideband SF)
        kZbb        = 2,  // signal
    #if EASYBJETS_ANY_MINOR_BKGS
        kTop        = 3,   // tty (DB) / top (MJB)
        kWqq        = 4,   // Wqqy     / WqqJets
        kZqq        = 5,   // Zqqy     / ZqqJets
        kNSamples   = 6
    #else
        kNSamples   = 3
    #endif
    };

    // ------------------------------------------------------------------------
    // PER-SAMPLE FILENAME TABLE (edit here when a sample's filename changes)
    // ------------------------------------------------------------------------
    struct SampleName {
        std::string run2;
        std::string run3;
        SampleName(const std::string& s) : run2(s), run3(s) {}
        SampleName(const char* s)        : run2(s), run3(s) {}
        SampleName(const std::string& r2, const std::string& r3) : run2(r2), run3(r3) {}

        std::string resolve() const {
            return RunConfig::IsRun3() ? run3 : run2;
        }
    };

    // Edit here when a sample gets renamed; nothing else needs to change. 
    namespace FileNamesDB {
        inline const SampleName kDataFile      = {"data"};
        inline const SampleName kMainBkgFile   = {"gamma_jets_sherpa", "gamma_jets_pythia"};
        inline const SampleName kSignalFile    = {"Zbby"};
        inline const SampleName kTopFile       = {"tty"};  
        inline const SampleName kWqqFile       = {"Wqqy"}; 
        inline const SampleName kZqqFile       = {"Zqqy"};  
    }
    namespace FileNamesMJB {
        inline const SampleName kDataFile      = {"data"};
        inline const SampleName kMainBkgFile   = {"multijets_pythia", "multijets"};
        inline const SampleName kSignalFile    = {"ZbbJets"};
        inline const SampleName kTopFile       = {"top"};
        inline const SampleName kWqqFile       = {"WqqJets"};
        inline const SampleName kZqqFile       = {"ZqqJets"};
    }

    inline std::vector<std::string> FileNames(const std::string& channel = "DB") {
        if (channel == "MJB")
            return {FileNamesMJB::kDataFile.resolve(),
                    FileNamesMJB::kMainBkgFile.resolve(),
                    FileNamesMJB::kSignalFile.resolve()
    #if EASYBJETS_ANY_MINOR_BKGS
        #if EASYBJETS_MJB_MINOR_BKGS
                   , FileNamesMJB::kTopFile.resolve()
                   , FileNamesMJB::kWqqFile.resolve()
                   , FileNamesMJB::kZqqFile.resolve()
        #else
                   , "__unused_top", "__unused_Wqq", "__unused_Zqq"
        #endif
    #endif
            };
        return {FileNamesDB::kDataFile.resolve(),
                FileNamesDB::kMainBkgFile.resolve(),
                FileNamesDB::kSignalFile.resolve()
    #if EASYBJETS_ANY_MINOR_BKGS
        #if EASYBJETS_DB_MINOR_BKGS
                   , FileNamesDB::kTopFile.resolve()
                   , FileNamesDB::kWqqFile.resolve()
                   , FileNamesDB::kZqqFile.resolve()
        #else
                   , "__unused_tty", "__unused_Wqqy", "__unused_Zqqy"
        #endif
    #endif
        };
    }
    inline std::vector<std::string> FileNames(const VariantConfig& cfg) {
        return FileNames(cfg.channel);
    }

    // Legend display names
    inline std::vector<std::string> DisplayNames(const std::string& channel = "DB") {
        if (channel == "MJB")
            return {"Data", "Multijets", "Z(#rightarrow b#bar{b})+jets"
    #if EASYBJETS_ANY_MINOR_BKGS
        #if EASYBJETS_MJB_MINOR_BKGS
                   , "t#bar{t}", "W(qq)+jets", "Z(qq)+jets"
        #else
                   , "", "", ""
        #endif
    #endif
            };
        return {"Data", "#gamma+jets", "Z(#rightarrow b#bar{b})+#gamma"
    #if EASYBJETS_ANY_MINOR_BKGS
        #if EASYBJETS_DB_MINOR_BKGS
                   , "t#bar{t}+#gamma", "W(qq)+#gamma", "Z(qq)+#gamma"
        #else
                   , "", "", ""
        #endif
    #endif
        };
    }
    inline std::vector<std::string> DisplayNames(const VariantConfig& cfg) {
        return DisplayNames(cfg.channel);
    }

    // Histogram name suffixes 
    inline std::vector<std::string> HistNames(const std::string& channel = "DB") {
        if (channel == "MJB")
            return {"data", "MultiJets", "Zbbjets"
    #if EASYBJETS_ANY_MINOR_BKGS
        #if EASYBJETS_MJB_MINOR_BKGS
                   , "top", "WqqJets", "ZqqJets"
        #else
                   , "unused_top", "unused_WqqJ", "unused_ZqqJ"
        #endif
    #endif
            };
        return {"data", "GammaJets", "Zbby"
    #if EASYBJETS_ANY_MINOR_BKGS
        #if EASYBJETS_DB_MINOR_BKGS
                   , "tty", "Wqqy", "Zqqy"
        #else
                   , "unused_tty", "unused_Wqqy", "unused_Zqqy"
        #endif
    #endif
        };
    }
    inline std::vector<std::string> HistNames(const VariantConfig& cfg) {
        return HistNames(cfg.channel);
    }

    inline std::vector<int> Colors() {
        return {kBlack, kBlue, kRed
    #if EASYBJETS_ANY_MINOR_BKGS
               , kOrange + 7   // minor bkg 1: kTop
               , kGreen  + 2   // minor bkg 2: kWqq
               , kViolet + 2   // minor bkg 3: kZqq
    #endif
        };
    }
    inline std::vector<int> FillStyles() {
        return {0, 1001, 1001
    #if EASYBJETS_ANY_MINOR_BKGS
               , 1001, 1001, 1001
    #endif
        };
    }

    // Core samples (0..2) are required; anything above that is optional and
    // handled by the tree_extra[] path (missing-file tolerant).
    inline bool IsOptional(int idx) { return idx >= 3; }

    inline std::vector<int> MCIndices() {
        std::vector<int> mc;
        mc.push_back(kMultijets);
        mc.push_back(kZbb);
        for (int i = 3; i < kNSamples; ++i) mc.push_back(i);
        return mc;
    }

    // Build all input file paths for a given variant + working point
    inline std::vector<std::string> AllSamplePaths(const VariantConfig& cfg,
                                                    const std::string& wp) {
        auto names = FileNames(cfg);
        std::vector<std::string> paths;
        for (const auto& n : names)
            paths.push_back(PathBuilder::SamplePath(cfg, n, wp));
        return paths;
    }

    // Stack order: minor backgrounds (smallest contributions at the bottom)
    // then main background (dominant) then Z(bb) signal on top.
    // This is the order in which the histograms are added to the stack.
    inline std::vector<int> StackOrder() {
        std::vector<int> order;
        for (int i = 3; i < kNSamples; ++i) order.push_back(i);  
        order.push_back(kMultijets);
        order.push_back(kZbb);
        return order;
    }

    inline void BuildStack(THStack* hs, TH1F** histograms, const bool sampleLoaded[]) {
        for (int idx : StackOrder())
            if (sampleLoaded[idx] && histograms[idx]) hs->Add(histograms[idx]);
    }

    inline void AddLegendEntries(TLegend* legend, TH1F** h,
                                  TGraphAsymmErrors* err_band,
                                  const bool sampleLoaded[],
                                  const VariantConfig& cfg) {
        auto dn = DisplayNames(cfg);
        if (sampleLoaded[kData]      && h[kData])      legend->AddEntry(h[kData],      "Data",                "lep");
        if (sampleLoaded[kZbb]       && h[kZbb])       legend->AddEntry(h[kZbb],       dn[kZbb].c_str(),      "f");
        if (sampleLoaded[kMultijets] && h[kMultijets]) legend->AddEntry(h[kMultijets], dn[kMultijets].c_str(),"f");
        // Minor backgrounds (indices >= 3) - only present when kNSamples > 3.
        for (int i = 3; i < kNSamples && i < (int)dn.size(); ++i) {
            if (sampleLoaded[i] && h[i])
                legend->AddEntry(h[i], dn[i].c_str(), "f");
        }
        if (err_band) legend->AddEntry(err_band, "MC Stat. Unc.", "f");
    }

    inline void AddLegendEntries(TLegend* legend, TH1F** h,
                                  TGraphAsymmErrors* err_band,
                                  const bool sampleLoaded[]) {
        auto dn = DisplayNames("DB");
        if (sampleLoaded[kData]      && h[kData])      legend->AddEntry(h[kData],      "Data",                "lep");
        if (sampleLoaded[kZbb]       && h[kZbb])       legend->AddEntry(h[kZbb],       dn[kZbb].c_str(),      "f");
        if (sampleLoaded[kMultijets] && h[kMultijets]) legend->AddEntry(h[kMultijets], dn[kMultijets].c_str(),"f");
        for (int i = 3; i < kNSamples && i < (int)dn.size(); ++i) {
            if (sampleLoaded[i] && h[i])
                legend->AddEntry(h[i], dn[i].c_str(), "f");
        }
        if (err_band) legend->AddEntry(err_band, "MC Stat. Unc.", "f");
    }

    // Total MC yield at a given cutflow bin
    inline double TotalMC(TH1F** h, const bool sampleLoaded[], int bin) {
        double total = 0;
        for (int i = 0; i < kNSamples; ++i)
            if (i != kData && sampleLoaded[i] && h[i])
                total += h[i]->GetBinContent(bin);
        return total;
    }

    // Create a summed MC histogram 
    inline TH1F* MakeTotalMC(TH1F** h, const bool sampleLoaded[],
                               const std::string& name = "h_sumMC") {
        TH1F* total = nullptr;
        for (int i = 0; i < kNSamples; ++i) {
            if (i != kData && sampleLoaded[i] && h[i]) {
                total = static_cast<TH1F*>(h[i]->Clone(name.c_str()));
                total->SetDirectory(nullptr);
                for (int j = i + 1; j < kNSamples; ++j)
                    if (j != kData && sampleLoaded[j] && h[j])
                        total->Add(h[j]);
                break;
            }
        }
        return total;
    }

}  // namespace SampleConfig


// ============================================================================
// PATH BUILDER — Signal/Background convenience wrappers (defined here because
// they depend on SampleConfig filenames).
// ============================================================================
namespace PathBuilder {
    inline std::string SignalPath(const VariantConfig& cfg, const std::string& wp) {
        auto names = SampleConfig::FileNames(cfg);
        return SamplePath(cfg, names[SampleConfig::kZbb], wp);
    }
    inline std::string BkgPath(const VariantConfig& cfg, const std::string& wp) {
        auto names = SampleConfig::FileNames(cfg);
        return SamplePath(cfg, names[SampleConfig::kMultijets], wp);
    }
}


// ============================================================================
// DATASET CONFIG - the part of the VariantConfig that Mass_Analysis uses 
// ============================================================================
struct DatasetConfig {
    std::string name;
    std::string mass_branch;
    std::string pt_branch;
    std::string output_subdir;
    double xmin, xmax;
    int    sideband2_max;
    int    min_bins, max_bins;
    double sideband1_low, sideband1_high, sideband2_low;
    std::vector<int> min_bins_perPtBin;
    std::vector<int> max_bins_perPtBin;

    int GetMinBins(int ptBinIdx) const {
        if (ptBinIdx >= 0 && ptBinIdx < (int)min_bins_perPtBin.size())
            return min_bins_perPtBin[ptBinIdx];
        return min_bins;
    }
    int GetMaxBins(int ptBinIdx) const {
        if (ptBinIdx >= 0 && ptBinIdx < (int)max_bins_perPtBin.size())
            return max_bins_perPtBin[ptBinIdx];
        return max_bins;
    }

    static DatasetConfig FromVariant(const VariantConfig& cfg) {
        DatasetConfig d;
        d.name           = cfg.plot_subdir;
        d.mass_branch    = cfg.branch_mass;
        d.pt_branch      = cfg.branch_pt;
        d.output_subdir  = cfg.plot_subdir;
        d.xmin           = cfg.mass_xmin;
        d.xmax           = cfg.mass_xmax;
        d.sideband2_max  = cfg.sideband2_max;
        d.min_bins       = cfg.min_bins;
        d.max_bins       = cfg.max_bins;
        d.sideband1_low  = cfg.sideband1_low;
        d.sideband1_high = cfg.sideband1_high;
        d.sideband2_low  = cfg.sideband2_low;
        d.min_bins_perPtBin = cfg.min_bins_perPtBin;
        d.max_bins_perPtBin = cfg.max_bins_perPtBin;
        return d;
    }
};

inline std::map<std::string, DatasetConfig> GetDatasetConfigs(const std::string& channel = "DB") {
    std::map<std::string, DatasetConfig> datasets;
    for (const auto& key : GetPrimaryVariantKeys(channel))
        datasets[key] = DatasetConfig::FromVariant(GetVariant(key));
    return datasets;
}


// ============================================================================
// pT BIN CONFIGURATION
// ============================================================================
// DB (photon recoil, lower pT):  2 bins — 200-300, 300-450
// MJB (multijet recoil, higher pT): 3 bins — 450-600, 600-900, 900-1200
//
// Two sets of edges per definition:
//   bins[]         used for the event selection; the last edge is set to a very
//                  large value (1e9) so the highest-pT bin is open and nothing is dropped
//   display_bins[] used for axes, labels and filenames; the last edge is the
//                  physical one (450 for DB, 1200 for MJB)
// They are identical except for that last edge.
// ============================================================================
namespace PtBinConfig {

    struct PtBinDef {
        int    nBins;
        double bins[9];            
        double display_bins[9];    
        const char* labels[8];     
        const char* fileNames[8];  
    };

    constexpr double kPtBinOpenEnd = 1.0e9;

    inline PtBinDef DefaultBins_DB() {
        PtBinDef d;
        d.nBins       = 2;
        d.bins[0]     = 200; d.bins[1] = 300; d.bins[2] = kPtBinOpenEnd;
        d.display_bins[0] = 200; d.display_bins[1] = 300; d.display_bins[2] = 450;
        d.labels[0]   = "200 < p_{T} < 300 GeV";
        d.labels[1]   = "300 < p_{T} < 450 GeV";
        d.fileNames[0]= "pt200to300";
        d.fileNames[1]= "pt300to450";
        return d;
    }

    // Two bins. We tried three (450/600/900/1200) but the last one was too
    // poorly populated for the mass fit.
    inline PtBinDef DefaultBins_MJB() {
        PtBinDef d;
        d.nBins        = 2;
        d.bins[0]      = 450; d.bins[1] = 750; d.bins[2] = kPtBinOpenEnd; //= 900;
        d.display_bins[0] = 450; d.display_bins[1] = 750; d.display_bins[2] = 1200;
        d.labels[0]    = "450 < p_{T} < 750 GeV";
        d.labels[1]    = "750 < p_{T} < 1200 GeV";
        d.fileNames[0] = "pt450to750";
        d.fileNames[1] = "pt750to1200";
        return d;
    }

    // Inclusive (single-bin) captures the entire pT spectrum without losing events
    // above the displayed nominal upper edge.
    inline PtBinDef InclusiveBins_DB() {
        PtBinDef d;
        d.nBins       = 1;
        d.bins[0]     = 200; d.bins[1] = kPtBinOpenEnd;
        d.display_bins[0] = 200; d.display_bins[1] = 450;
        d.labels[0]   = "200 < p_{T} < 450 GeV";
        d.fileNames[0]= "pt200to450";
        return d;
    }

    inline PtBinDef InclusiveBins_MJB() {
        PtBinDef d;
        d.nBins        = 1;
        d.bins[0]      = 450; d.bins[1] = kPtBinOpenEnd;
        d.display_bins[0] = 450; d.display_bins[1] = 1200;
        d.labels[0]    = "450 < p_{T} < 1200 GeV";
        d.fileNames[0] = "pt450to1200";
        return d;
    }

    inline PtBinDef DefaultBins(const std::string& channel = "DB") {
        return (channel == "MJB") ? DefaultBins_MJB() : DefaultBins_DB();
    }
    inline PtBinDef DefaultBins(const VariantConfig& cfg) {
        return DefaultBins(cfg.channel);
    }

    inline PtBinDef InclusiveBins(const std::string& channel = "DB") {
        return (channel == "MJB") ? InclusiveBins_MJB() : InclusiveBins_DB();
    }
    inline PtBinDef InclusiveBins(const VariantConfig& cfg) {
        return InclusiveBins(cfg.channel);
    }

    inline PtBinDef GetBins(const std::string& mode = "ptbins",
                              const std::string& channel = "DB") {
        if (mode == "inclusive" || mode == "1bin") return InclusiveBins(channel);
        return DefaultBins(channel);
    }
    inline PtBinDef GetBins(const VariantConfig& cfg,
                              const std::string& mode = "ptbins") {
        return GetBins(mode, cfg.channel);
    }

}  // namespace PtBinConfig


// ============================================================================
// ANALYSIS-SPECIFIC WORKING POINT SETS
// ============================================================================
namespace AnalysisWPConfig {

    inline std::vector<std::string> StandardWPs() {
        return {"46"};
        //return {"25", "30", "37", "46", "58", "74", "94", "125", "155"};
    }

    inline std::vector<std::string> PrimaryWP() { return {"46"}; }

    inline std::vector<int> StandardColors() {
        return {kRed, kBlue, kGreen+2, kMagenta, kCyan+2,
                kOrange+7, kViolet+3, kPink+2, kSpring+4};
    }

}  // namespace AnalysisWPConfig


#endif // VARIANT_CONFIG_H
