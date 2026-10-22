/**
 * ============================================================================
 *  WorkingPointConfig.h - Single source of truth for (bb)-tagging working points
 * ============================================================================
 * ===========================================
 * @author:    Thiziri AMEZZA
 * @email:     thiziri.amezza@cern.ch
 * @date:      10 March 2026 
 * ===========================================
 * 
 *  A working point appears in THREE different spellings across the framework:
 *
 *    1. key            - what code, scripts, folders and filenames use.
 *                        Standardised: "0p25", "0p30", "0p37", "0p46", ...
 *    2. branch suffix  - fixed by the input ntuple schema. Identical to the
 *                        key for every WP EXCEPT 0.30, whose branch is
 *                        "pass_FlatMassQCDEff_0p3" (note: "0p3", not "0p30").
 *    3. display label  - for plot legends, e.g. "GN2X (0.46% QCD WP)".
 *
 * To add, remove, or rename a working point, edit ONLY WorkingPoint::Table() below.
 * ==========================================================================
 * @version:  1.0
 * ==========================================================================
 */

#ifndef EASYBJETS_WORKING_POINT_CONFIG_H
#define EASYBJETS_WORKING_POINT_CONFIG_H

#include <string>
#include <vector>

namespace WorkingPoint {

    struct Info {
        std::string key;     
        std::string branch;  
        std::string label;   
    };

    // ----------------------------------------------------------------------
    // Edit ONLY this table to change the WP set.
    // ----------------------------------------------------------------------
    inline const std::vector<Info>& Table() {
        static const std::vector<Info> t = {
            {"0p25", "0p25", "GN2X (0.25% QCD WP)"},
            {"0p30", "0p3",  "GN2X (0.30% QCD WP)"},  
            {"0p37", "0p37", "GN2X (0.37% QCD WP)"},
            {"0p46", "0p46", "GN2X (0.46% QCD WP)"},
            {"0p58", "0p58", "GN2X (0.58% QCD WP)"},
            {"0p74", "0p74", "GN2X (0.74% QCD WP)"},
            {"0p94", "0p94", "GN2X (0.94% QCD WP)"},
            {"1p25", "1p25", "GN2X (1.25% QCD WP)"},
            {"1p55", "1p55", "GN2X (1.55% QCD WP)"},
        };
        return t;
    }

    // Accepts the key, the branch suffix or a bare number ("30", "46", "125")
    // and returns the key. Empty string if the WP is not in the table.
    // e.g. "30" -> "0p30"
    inline std::string Canonical(const std::string& wp) {
        std::string cand = wp;
        if (wp.find('p') == std::string::npos) {     // bare number -> XpYY form
            if      (wp == "125") cand = "1p25";
            else if (wp == "155") cand = "1p55";
            else if (wp == "3" || wp == "30") cand = "0p30";
            else if (!wp.empty()) cand = "0p" + wp;  
        }
        for (const auto& i : Table())
            if (cand == i.key || wp == i.branch || wp == i.key) return i.key;
        return "";
    }

    inline const Info* Find(const std::string& wp) {
        const std::string k = Canonical(wp);
        if (k.empty()) return nullptr;
        for (const auto& i : Table())
            if (i.key == k) return &i;
        return nullptr;
    }

    inline bool IsValid(const std::string& wp) { return Find(wp) != nullptr; }

    inline std::string FolderLabel (const std::string& wp) { auto* i = Find(wp); return i ? i->key    : wp; }
    inline std::string BranchSuffix(const std::string& wp) { auto* i = Find(wp); return i ? i->branch : wp; }
    inline std::string DisplayLabel(const std::string& wp) { auto* i = Find(wp); return i ? i->label  : wp; }

    inline std::string BranchName(const std::string& wp) {
        return "pass_FlatMassQCDEff_" + BranchSuffix(wp);
    }

    inline std::vector<std::string> AllKeys() {
        std::vector<std::string> v;
        v.reserve(Table().size());
        for (const auto& i : Table()) v.push_back(i.key);
        return v;
    }

}  // namespace WorkingPoint

#endif  // EASYBJETS_WORKING_POINT_CONFIG_H
