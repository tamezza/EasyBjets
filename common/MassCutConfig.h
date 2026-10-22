/**
 * ============================================================================
 *   MassCutConfig.h - optional global cut on the large-R jet mass
 * ============================================================================
 * ===========================================
 * @author:    Thiziri AMEZZA
 * @email:     thiziri.amezza@cern.ch
 * @date:      10 March 2026 
 * ===========================================
 *
 * All event loops call MassCutConfig::Reject() right after validating the
 * selected jet index. With EASYBJETS_MASS_CUT unset nothing is cut and the
 * macros behave exactly as before this was added. Switching it on keeps only
 * events with m_low <= m_J <= m_high, which removes most of the top
 * contamination at high mass without touching the balance observables.
 *
 * Activate the cut via:
 * 
 *   export EASYBJETS_MASS_CUT=on        50-150 GeV (what we set)
 *   export EASYBJETS_MASS_CUT=off       explicit off (default)
 *   export EASYBJETS_MASS_CUT=50:150    explicit window 
 *   export EASYBJETS_MASS_CUT=40:170    custom - any window 
 *
 * If the window is narrower than the mass-fit range, the upper sideband of
 * the mass fit should be narrowed as well, e.g. with
 * EASYBJETS_SIDEBAND2_MAX=150 (see VariantConfig.h), otherwise the SF is
 * fitted on a region that no longer contains events.
 * ==========================================================================
 * @version:  1.0
 * ==========================================================================
 */

#ifndef EASYBJETS_MASS_CUT_CONFIG_H
#define EASYBJETS_MASS_CUT_CONFIG_H

#include <cctype>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

namespace MassCutConfig {

// Default window when EASYBJETS_MASS_CUT=on (no explicit lo:hi)
constexpr double kDefaultLow  = 50.0;
constexpr double kDefaultHigh = 150.0;

struct Settings {
    bool   active = false;
    double low   = kDefaultLow;
    double high  = kDefaultHigh;
};

// Parsed once per process. Unparseable values disable the cut with a warning.
inline const Settings& Current() {
    static const Settings s = []() -> Settings {
        Settings r;
        const char* env = std::getenv("EASYBJETS_MASS_CUT");
        if (!env || !*env) return r;   // unset → disabled

        std::string raw(env);
        std::string lc(raw);
        for (auto& c : lc) c = (char)std::tolower((unsigned char)c);

        if (lc == "off" || lc == "0" || lc == "false" || lc == "none") {
            return r;
        }
        if (lc == "on" || lc == "1" || lc == "true" || lc == "default") {
            r.active = true;
            r.low    = kDefaultLow;
            r.high   = kDefaultHigh;
            return r;
        }

        size_t sep = raw.find(':');
        if (sep == std::string::npos) sep = raw.find(',');
        if (sep != std::string::npos) {
            try {
                double lo = std::stod(raw.substr(0, sep));
                double hi = std::stod(raw.substr(sep + 1));
                if (hi > lo && lo >= 0) {
                    r.active = true;
                    r.low    = lo;
                    r.high   = hi;
                    return r;
                }
            } catch (...) {}
        }
        std::cerr << "[MassCutConfig] WARN: cannot parse EASYBJETS_MASS_CUT='"
                  << env << "' (expected on|off|<lo>:<hi>); cut DISABLED\n";
        return r;
    }();
    return s;
}

inline bool   IsActive() { return Current().active; }
inline double Low()      { return Current().low;    }
inline double High()     { return Current().high;   }

inline std::string Tag() {
    const auto& s = Current();
    if (!s.active) return "off";
    char buf[64];
    std::snprintf(buf, sizeof(buf), "%g:%g", s.low, s.high);
    return std::string(buf);
}

// Print a one-line at the start of each macro so the log says whether the cut was on.
inline void PrintBanner(const std::string& macroName) {
    const auto& s = Current();
    if (s.active) {
        std::cout << "[MassCutConfig] " << macroName
                  << ": event-level mass cut ACTIVE — accepting "
                  << s.low << " ≤ ljet_m ≤ " << s.high << " GeV\n";
    } else {
        std::cout << "[MassCutConfig] " << macroName
                  << ": event-level mass cut OFF (default)\n";
    }
}

// ============================================================================
// Per-event rejection helper — the SINGLE entry point for event loops.
//
// Returns true if the event should be SKIPPED (i.e., its mass falls outside
// the configured window). Returns false if the event should be processed.
// ============================================================================
inline bool Reject(const std::vector<float>* ljet_m, int idx) {
    const auto& s = Current();
    if (!s.active) return false;
    if (!ljet_m || idx < 0 || idx >= (int)ljet_m->size()) return true;
    double m = (*ljet_m)[idx];
    return (m < s.low || m > s.high);
}

inline bool RejectScalar(double mass) {
    const auto& s = Current();
    if (!s.active) return false;
    return (mass < s.low || mass > s.high);
}

}  // namespace MassCutConfig

#endif  // EASYBJETS_MASS_CUT_CONFIG_H
