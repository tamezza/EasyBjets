/**
 * ============================================================================
 *  Insitu_Analysis.h — shared scaffolding for the in-situ calibration analysis
 * ============================================================================
 * ===========================================
 * @author:    Thiziri AMEZZA
 * @email:     thiziri.amezza@cern.ch
 * @date:      10 March 2026 
 * ===========================================
 * 
 * Zbby_Insitu_Analysis.cxx (DB) and ZbbJets_Insitu_Analysis.cxx (MJB) have
 * different selections, branches and histograms, so their event loops stay
 * in their own files. The bits that were copy-pasted between them (jet type
 * parsing, delta phi, working point check, output paths) are here.
 *
 * Output layout, with <era> = NtupleSlim_Run2 or NtupleSlim_Run3 (EASYBJETS_RUN):
 *   ../DataSets/<topology>/<era>/<channel>_<jet|MuonCorr>_<wp>wp_slim/
 *       <input basename>_<jet|MuonCorr>_<wp>wp_slim.root
 *       CutFlow[_<jet>|_MuonCorr]/
 *       BkgComposition/
 *
 *  Each driver includes this header and passes its own topology ("ZbbGamma" /
 *  "ZbbJets") and channel ("DB" / "MJB") into MakeInsituPaths().
 * ==========================================================================
 * @version:  1.0
 * ==========================================================================
 */

#ifndef EASYBJETS_INSITU_ANALYSIS_H
#define EASYBJETS_INSITU_ANALYSIS_H

#include <string>
#include <iostream>
#include <cmath>

#include "Utils/Util.h"          
#include "RunEra.h"              
#include "WorkingPointConfig.h"  

// ============================================================================
// Jet types (shared by both topologies)
// ============================================================================
static const double PI = std::acos(-1.0);

enum class JetType { NOMINAL, BJER10V00, BJER10V01 };

inline JetType parseJetType(const std::string& jetTypeStr) {
    if (jetTypeStr == "bJER10v00") return JetType::BJER10V00;
    if (jetTypeStr == "bJER10v01") return JetType::BJER10V01;
    return JetType::NOMINAL;
}

inline std::string getJetTypeSuffix(JetType jetType) {
    switch (jetType) {
        case JetType::BJER10V00: return "bJER10v00";
        case JetType::BJER10V01: return "bJER10v01";
        default:                 return "nominal";
    }
}

// Signed difference in phi, wrapped to (-pi, pi].
inline double DeltaPhi(double phi1, double phi2) {
    double dphi = phi1 - phi2;
    while (dphi >  PI) dphi -= 2.0 * PI;
    while (dphi <= -PI) dphi += 2.0 * PI;
    return dphi;
}

// ============================================================================
// Working-point normalisation: Validate and canonicalise a WP in place
// (e.g. "0p3" -> "0p30"). Returns false and logs an error if the WP is unknown.
inline bool NormalizeWorkingPoint(std::string& workingPoint) {
    if (!WorkingPoint::IsValid(workingPoint)) {
        std::cerr << "Error: Unknown working point '" << workingPoint << "'\n";
        return false;
    }
    workingPoint = WorkingPoint::Canonical(workingPoint);
    return true;
}

// ============================================================================
// Output paths
// ============================================================================
struct InsituPaths {
    Str baseName;       
    Str outputFolder;   
    Str fn_out;         
    Str cutFlowFolder;  
    Str bkgFolder;      
};

// topology: "ZbbGamma" | "ZbbJets",  channel: "DB" | "MJB"
// The era subdir (NtupleSlim_Run2|Run3) comes from RunConfig - set EASYBJETS_RUN.
inline InsituPaths MakeInsituPaths(const Str& fn,
                                   const std::string& topology,
                                   const std::string& channel,
                                   const std::string& jetTypeSuffix,
                                   const std::string& workingPoint,
                                   bool isMuonCorr) {
    InsituPaths p;
    p.baseName = gSystem->BaseName(fn);

    const std::string folderTag = isMuonCorr ? "MuonCorr" : jetTypeSuffix;
    p.outputFolder = "../DataSets/" + topology + "/" + RunConfig::NtupleSubdir()
                   + "/" + channel + "_" + folderTag + "_" + workingPoint + "wp_slim/";

    const std::string outSuffix = isMuonCorr
        ? ("_MuonCorr_" + workingPoint + "wp_slim.root")
        : ("_" + jetTypeSuffix + "_" + workingPoint + "wp_slim.root");

    const std::string cutflowDirName = isMuonCorr
        ? "CutFlow_MuonCorr"
        : (jetTypeSuffix == "nominal" ? "CutFlow" : "CutFlow_" + jetTypeSuffix);

    p.fn_out        = p.outputFolder + Str(p.baseName).ReplaceAll(".root", outSuffix.c_str());
    p.cutFlowFolder = p.outputFolder + cutflowDirName.c_str() + "/";
    p.bkgFolder     = p.outputFolder + "BkgComposition/";
    return p;
}

inline void CreateInsituDirs(const InsituPaths& p) {
    auto mk = [](const TString& dir) {
        if (gSystem->AccessPathName(dir.Data()))
            gSystem->Exec(("mkdir -p " + dir).Data());
    };
    mk(p.outputFolder);
    mk(p.cutFlowFolder);
    mk(p.bkgFolder);
}

#endif  // EASYBJETS_INSITU_ANALYSIS_H
