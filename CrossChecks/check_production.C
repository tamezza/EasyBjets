/**
 * ============================================================================
 *  check_production.C - sanity check of the slim ntuple production
 * ============================================================================
  * ===========================================
 * @author:    Thiziri AMEZZA
 * @email:     thiziri.amezza@cern.ch
 * @date:      18 April 2026 
 * ===========================================
 * 
 * The producers do not fail hard: the MJB run scripts ignore ROOT's segfault
 * exit code, so a crash in the middle of a file leaves a truncated .root and
 * the script reports success, and a selection that keeps no event writes an
 * empty tree without complaint. This macro opens every slim ntuple of an era
 * and checks that it can be read and has a non-empty "nominal" tree.
 *
 *   root 'check_production.C("run3","ALL")'
 *   root 'check_production.C("run2","DB")'
 *
 *   era      run2 | run3          (NtupleSlim_Run2 / NtupleSlim_Run3)
 *   channel  DB | MJB | ALL       (ZbbGamma / ZbbJets)
 *
 * Run from EasyBjets/ or any directory from which DataSets/ can be found.
 * Prints one line per file with the entry counts and flags
 *   [DIR MISSING]  [NO .root FILES]  [ZOMBIE]  [NO 'nominal' TREE]  EMPTY
 * then a summary with the number of problems.
 *
 * ==========================================================================
 * @version:  1.0
 * ==========================================================================
 */
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TError.h>
#include <vector>
#include <algorithm>
#include <cstdio>

namespace {
    int g_problems = 0;
    int g_files    = 0;
    long long g_totEntries = 0;

    // Resolve DataSets/ whether we are run from EasyBjets/, a macros/ subdir,
    // or the repo root. Returns "" if none found.
    TString FindDataSetsBase() {
        const char* cands[] = {"DataSets", "../DataSets", "../../DataSets",
                               "EasyBjets/DataSets"};
        for (const char* c : cands)
            if (!gSystem->AccessPathName(c)) return TString(c);
        return TString("");
    }

    void CheckOneFile(const TString& path) {
        ++g_files;
        const char* bn = gSystem->BaseName(path);
        TFile* f = TFile::Open(path, "READ");
        if (!f || f->IsZombie()) {
            printf("    [ZOMBIE/UNOPENABLE] %s\n", bn);
            ++g_problems;
            if (f) delete f;
            return;
        }
        TTree* t = (TTree*)f->Get("nominal");
        if (!t) {
            printf("    [NO 'nominal' TREE]  %s\n", bn);
            ++g_problems;
            f->Close(); delete f;
            return;
        }
        Long64_t n  = t->GetEntries();
        g_totEntries += n;

        // DB has a second tree for the mass window, but not all files do (-1 means "not present" in MJB).
        TTree* tm   = (TTree*)f->Get("CutMassNominal");
        Long64_t nm = tm ? tm->GetEntries() : -1;

        const char* flag = (n == 0) ? "   <== EMPTY 'nominal' tree!" : "";
        if (n == 0) ++g_problems;

        if (nm >= 0)
            printf("    %-52s nominal=%-9lld  CutMassNominal=%-9lld%s\n",
                   bn, (long long)n, (long long)nm, flag);
        else
            printf("    %-52s nominal=%-9lld%s\n", bn, (long long)n, flag);

        f->Close(); delete f;
    }

    void ScanFolder(const TString& dir) {
        void* d = gSystem->OpenDirectory(dir);
        if (!d) { printf("    [DIR MISSING] %s\n", dir.Data()); ++g_problems; return; }
        std::vector<TString> files;
        const char* e;
        while ((e = gSystem->GetDirEntry(d))) {
            TString name(e);
            if (name.EndsWith(".root")) files.push_back(dir + "/" + name);
        }
        gSystem->FreeDirectory(d);
        std::sort(files.begin(), files.end());
        if (files.empty()) { printf("    [NO .root FILES]\n"); ++g_problems; return; }
        for (const auto& fp : files) CheckOneFile(fp);
    }
}

void check_production(const char* era = "run2", const char* channel = "ALL") {
    gErrorIgnoreLevel = kError;  

    TString base = FindDataSetsBase();
    if (base.IsNull()) {
        printf("ERROR: could not locate a DataSets/ directory from here (%s).\n",
               gSystem->pwd());
        printf("Run this from the EasyBjets/ directory.\n");
        return;
    }

    TString eraSub = TString(era).Contains("3") ? "NtupleSlim_Run3"
                                                 : "NtupleSlim_Run2";

    std::vector<std::pair<TString,TString>> chans;  // {topology subdir, label}
    TString c(channel); c.ToUpper();
    if (c == "DB"  || c == "ALL") chans.push_back({"ZbbGamma", "DB"});
    if (c == "MJB" || c == "ALL") chans.push_back({"ZbbJets",  "MJB"});
    if (chans.empty()) {
        printf("ERROR: channel must be DB, MJB, or ALL (got '%s').\n", channel);
        return;
    }

    printf("\n#### Production check  era=%s  channel=%s  base=%s ####\n",
           eraSub.Data(), c.Data(), base.Data());

    for (const auto& ch : chans) {
        TString topoDir = base + "/" + ch.first + "/" + eraSub;
        printf("\n==================== %s   (%s) ====================\n",
               ch.second.Data(), topoDir.Data());
        void* d = gSystem->OpenDirectory(topoDir);
        if (!d) {
            printf("  [TOPOLOGY DIR MISSING] %s\n", topoDir.Data());
            ++g_problems;
            continue;
        }
        std::vector<TString> subs;
        const char* e;
        while ((e = gSystem->GetDirEntry(d))) {
            TString name(e);
            if (name.EndsWith("wp_slim")) subs.push_back(name);
        }
        gSystem->FreeDirectory(d);
        std::sort(subs.begin(), subs.end());
        if (subs.empty()) {
            printf("  [NO *wp_slim FOLDERS] — nothing produced for this era/channel\n");
            ++g_problems;
            continue;
        }
        for (const auto& s : subs) {
            printf("\n  -- %s --\n", s.Data());
            ScanFolder(topoDir + "/" + s);
        }
    }

    printf("\n========================================================\n");
    printf("Checked %d file(s), %lld total 'nominal' entries; %d problem(s).\n",
           g_files, g_totEntries, g_problems);
    if (g_problems > 0)
        printf("RESULT: PROBLEMS DETECTED — inspect the flagged rows above.\n");
    else
        printf("RESULT: OK — every file opens with a non-empty 'nominal' tree.\n");
    printf("========================================================\n");
}
