/**
 * ============================================================================
 *  AtlasLabel.h - the ATLAS label of every plot, set in one place
 * ============================================================================
 * ===========================================
 * @author:    Thiziri AMEZZA
 * @email:     thiziri.amezza@cern.ch
 * @date:      26 October 2026
 * ===========================================
 *
 * All plots take their "ATLAS <status>" label from here, so it is changed for
 * the whole framework without editing any macro.
 *
 *  USAGE
 *  -----
 *   Default: "ATLAS Internal". To change it, either
 *
 *    (A) set the environment variable before starting, e.g.
 *          export EASYBJETS_ATLAS_LABEL=Preliminary
 *          root 'Mass_Summary_Combined.cxx()'
 *
 *    (B) edit kDefaultStatus below.
 *
 *   EASYBJETS_ATLAS_LABEL values:
 *     Internal (default) | Preliminary | "Simulation Internal" | any text
 *     none               no ATLAS label at all (the luminosity line stays)
 *     ""                 "ATLAS" alone
 *
 *  In the code
 *  -----------
 *   AtlasLabel::Title()          "#font[72]{ATLAS} Internal", to draw with your own TLatex
 *   AtlasLabel::Title(0.79)      same, with the status scaled to 79%
 *   AtlasLabel::Draw(x, y)       label + luminosity line below it
 *   AtlasUtils::DrawAtlasLabel   (VariantConfig.h) calls AtlasLabel::Draw
 *
 * ==========================================================================
 * @version:  1.0
 * ==========================================================================
 */

#ifndef EASYBJETS_ATLAS_LABEL_H
#define EASYBJETS_ATLAS_LABEL_H

#include "RunEra.h"  // RunConfig::LumiLabel

#include <TLatex.h>

#include <cstdio>
#include <cstdlib>
#include <string>

namespace AtlasLabel {

// Label used when EASYBJETS_ATLAS_LABEL is not set.
inline const char* kDefaultStatus = "Internal";

// True unless EASYBJETS_ATLAS_LABEL=none (or off).
inline bool Enabled() {
    const char* env = std::getenv("EASYBJETS_ATLAS_LABEL");
    if (!env) return true;
    const std::string s(env);
    return !(s == "none" || s == "None" || s == "NONE" || s == "off");
}

// Text after "ATLAS": Internal, Preliminary, ... Printed once per session.
inline std::string Status() {
    const char* env = std::getenv("EASYBJETS_ATLAS_LABEL");
    const std::string s = env ? std::string(env) : std::string(kDefaultStatus);
    static bool printed = false;
    if (!printed) {
        printed = true;
        if (Enabled())
            std::printf("[AtlasLabel] plots labelled \"ATLAS %s\" (EASYBJETS_ATLAS_LABEL)\n", s.c_str());
        else
            std::printf("[AtlasLabel] no ATLAS label on the plots (EASYBJETS_ATLAS_LABEL=none)\n");
    }
    return Enabled() ? s : std::string();
}

// "ATLAS Internal" or "ATLAS Preliminary" or ...
inline std::string Title(double statusScale = 1.0) {
    const std::string s = Status();
    if (!Enabled()) return "";
    std::string t = "#font[72]{ATLAS}";
    if (s.empty()) return t;
    if (statusScale != 1.0) {
        char buf[32];
        std::snprintf(buf, sizeof(buf), "%.2f", statusScale);
        return t + " #scale[" + buf + "]{" + s + "}";
    }
    return t + " " + s;
}

// Label at (x, y) in NDC, and the luminosity line lumiDy below it.
inline void Draw(double x, double y, double size = 0.04, bool withLumi = true, double lumiDy = -0.04) {
    TLatex l;
    l.SetNDC();
    l.SetTextFont(42);
    l.SetTextColor(kBlack);
    l.SetTextSize(size);
    const std::string t = Title();
    if (!t.empty()) l.DrawLatex(x, y, t.c_str());
    if (withLumi) l.DrawLatex(x, y + lumiDy, RunConfig::LumiLabel().c_str());
}

}  // namespace AtlasLabel

#endif  // EASYBJETS_ATLAS_LABEL_H
