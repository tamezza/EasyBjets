/**
 * ============================================================================
 *  RunEra.h - Configuration for Run 2 / Run 3 switch
 * ============================================================================
 * ===========================================
 * @author:    Thiziri AMEZZA
 * @email:     thiziri.amezza@cern.ch
 * @date:      18 April 2026 
 * ===========================================
 * 
 * Everything that depends on the data-taking period is read from here:
 * the slim-ntuple subdirectory, the plot output directory and the
 * luminosity label on the plots.
 * e.g:
 *    - Input ntuple directory: NtupleSlim_Run2 or NtupleSlim_Run3
 *    - Plot output directory: ../Run2Plots/ or ../Run3Plots/
 *    - ATLAS label luminosity + and centre-of-mass energy labels:
 *
 *  USAGE
 *  -----
 *   Run 2 is the default. To switch eras, you have two options, either
 *
 *    (A) edit the #define below to RUN_ERA_RUN3 and reload ROOT, OR
 *
 *    (B) set the environment variable before starting ROOT:
 *          export EASYBJETS_RUN=run2 (or run3)
 *          root 'Response_Fit_Inclusive.cxx("ALL")'
 * 
 * Option (B) - the env-var override preferred over the compile-time default, 
 * so you can switch between eras without editing this file.
 * If the variable is not set, the compile-time default is used.
 * ==========================================================================
 * @version:  1.0
 * ==========================================================================
 */

#ifndef EASYBJETS_RUN_ERA_H
#define EASYBJETS_RUN_ERA_H

#include <cstdlib>
#include <string>

// Default era; can be overridden with EASYBJETS_RUN.
#define RUN_ERA_RUN2 2
#define RUN_ERA_RUN3 3
#ifndef EASYBJETS_DEFAULT_RUN_ERA
#define EASYBJETS_DEFAULT_RUN_ERA RUN_ERA_RUN3
#endif

namespace RunConfig {

    enum class Era { Run2, Run3 };

    // Read once per session. Accepts "run2"/"2" and "run3"/"3"; anything
    // else falls back to the compile-time default.
    inline Era CurrentEra() {
        static Era cached = [](){
            const char* env = std::getenv("EASYBJETS_RUN");
            if (env && *env) {
                std::string s(env);
                for (auto& ch : s) ch = (char)std::tolower((unsigned char)ch);
                if (s == "run3" || s == "3") return Era::Run3;
                if (s == "run2" || s == "2") return Era::Run2;
            }
            return (EASYBJETS_DEFAULT_RUN_ERA == RUN_ERA_RUN3) ? Era::Run3 : Era::Run2;
        }();
        return cached;
    }

    inline bool IsRun2() { return CurrentEra() == Era::Run2; }
    inline bool IsRun3() { return CurrentEra() == Era::Run3; }

    inline std::string NtupleSubdir() {
        return IsRun3() ? "NtupleSlim_Run3" : "NtupleSlim_Run2";
    }

    // Suffix added to the plot directory when a non-Gaussian signal model is
    // selected with EASYBJETS_FIT_MODEL, so the alternative-model outputs do
    // not overwrite the Gaussian ones. Reads the variable directly rather
    // than going through FitUtils to keep this header free of includes.
    inline std::string FitModelDirTag() {
        const char* env = std::getenv("EASYBJETS_FIT_MODEL");
        if (!env || !*env) return "";
        std::string s(env);
        for (auto& c : s) c = (char)std::tolower((unsigned char)c);
        if (s == "bukin" || s == "bukinonly" || s == "bukin-only")          return "_bukin";
        if (s == "crystalball" || s == "cb" || s == "crystalballonly" || s == "cb-only") return "_cb";
        if (s == "dscb" || s == "doublecb" || s == "double-cb" || s == "dscbonly" ||
            s == "doublesidedcb" || s == "double-sided-cb")                 return "_dscb";
        if (s == "doublegauss" || s == "doublegaussian" || s == "dgauss" ||
            s == "2gauss" || s == "twogauss" || s == "double-gauss")        return "_dgauss";
        if (s == "all" || s == "bukin-cb" || s == "extended")               return "_allmodels";
        return "";   // Default Gaussian output for unknown values.
    }

    // Plot directory, relative to <topology>/macros/
    // The fit-model tag is added when a non-default model is selected.
    inline std::string PlotDir() {
        return (IsRun3() ? "../Run3Plots" : "../Run2Plots") + FitModelDirTag();
    }

    // Luminosity and centre-of-mass energy shown on ATLAS plots.
    inline std::string LumiLabel() {
        return IsRun3()
            ? "#sqrt{s} = 13.6 TeV, 159 fb^{-1}"
            : "#sqrt{s} = 13 TeV, 140 fb^{-1}";
    }

    inline std::string EraTag() {
        return IsRun3() ? "Run3" : "Run2";
    }

}  // namespace RunConfig

#endif  // EASYBJETS_RUN_ERA_H
