/**
 * ============================================================================
 *  Response_Summary_Combined.h
 * ============================================================================
 * ===========================================
 * @author:    Thiziri AMEZZA
 * @email:     thiziri.amezza@cern.ch
 * @date:      10 June 2026 
 * ===========================================
 *
 * Summary plots of the response analysis, DB and MJB on one pT axis, built
 * from the per-channel CSVs written by Response_Ratios_pTbins_Mcuts.h. For
 * each WP, mass cut and correction (bJER10v00, bJER10v01, MuonCorr):
 *
 *   R_JES / R_JER summary   top: data/MC single ratio, nominal and correction
 *                           bottom: double ratio correction / nominal
 *   b_JES / b_JER response  top: data and MC for the nominal calibration
 *                           bottom: data and MC for the correction
 *
 * Each also exists with the three mass cuts overlaid, and the
 * summary with all corrections together when several are selected. The
 * summary figures are drawn once per double-ratio error treatment
 * (uncorrelated, correlation-aware, full-correlation), each in its own folder.
 *
 * Usage (from ZbbGamma/macros or ZbbJets/macros; the plots go to that
 * channel's Run{2,3}Plots), after Response_Ratios_pTbins_Mcuts.cxx("ALL") in
 * both channels with the same EASYBJETS_RUN:
 * 
 *   root 'Response_Summary_Combined.cxx()'                       all corrections, DB + MJB
 *   root 'Response_Summary_Combined.cxx("bJER10v01,MuonCorr")'   a subset
 *   root 'Response_Summary_Combined.cxx("ALL","DB")'             one channel only
 *   arguments: corrections ("ALL" or a comma list), topology (combined | DB | MJB)
 *
 * Input: the CSV written by Response_Ratios_pTbins_Mcuts.cxx, one per
 * channel and correction, e.g. for DB and bJER10v01 at WP 0p46:
 *   ../../ZbbGamma/Run3Plots/Response/Ratios/DB_Ratio_DB_bJER10v01/0p46/R2pTbin/
 *       Summary_Ratios_0p46.csv
 *
 * Output: <PlotDir>/Response/Combined_Ratios/<wp>/
 *   Summary_DataMC/<treatment>/R_JES|R_JER_Summary_<topo>_<corr|ALL>_<mass>_<wp>.pdf
 *   Summary_DataMC/<treatment>/R_JES|R_JER_Summary_<topo>_<corr>_AllMass_<wp>.pdf
 *   Summary_Response/bJES|bJER_Response_<topo>_<corr>_<mass|AllMass>_<wp>.pdf
 *   <treatment> = uncorrelated | correlation_aware | full_correlation
 *
 * Option: EASYBJETS_YRANGE_<JES|JER>_<SR|DR|COMP>="lo,hi" fixes the y range of
 * a panel, e.g. to give Run 2 and Run 3 the same scale (SR: single ratio,
 * DR: double ratio, COMP: response component -- b_JES / b_JER response).
 *
 * Shared helpers. Two other summaries reuse functions from this file (namespace
 * RatioResponseCombined::Summary). If you change one of these functions, re-run
 * those summaries too and check their plots:
 *
 *   Mass_Summary_Combined.h (mass summary) uses:
 *     CorrectionColor, CorrDisp          colour and displayed name of a correction
 *     ParseCorrections                   "ALL" or "a,b,c" -> list of corrections
 *     ExpandRange, ApplySummaryYRange    y range of the panels
 *
 *   Response_Summary_Combined_Inclusive.h (inclusive and pT-binned response
 *   summaries) uses:
 *     CorrectionColor, CorrDisp, ParseCorrections, ExpandRange
 *     QtyAxis, DoubleRatioAxis, RespAxis            axis titles
 *     SingleNom/Var, DoubleRatio, DoubleRatErr      values read from a CSV row
 *     DRErrFlavor, DRErrTag                         double-ratio error treatment
 *     NonClosureHalfWidth 
 *
 * ==========================================================================
 * @version:  1.0
 * ==========================================================================
 */

#ifndef RATIO_RESPONSE_COMBINED_SUMMARY_H
#define RATIO_RESPONSE_COMBINED_SUMMARY_H

#include "Response_Ratios_Combined.h" 
#include <TGraphErrors.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <TPaveText.h>
#include <functional>
#include <algorithm>
#include <sstream>

namespace RatioResponseCombined {
namespace Summary {

// Which quantity a figure shows: the scale (JES) or the resolution (JER).
enum class Qty { JES, JER };

// How the double-ratio error is computed. The CSV holds three versions of it
// and the summary draws a full set of figures for each (one folder per
// treatment). The one in use is set with DRErrFlavor() and read by DoubleRatErr.
//   Uncorr  errors added in quadrature, no correlation (largest error)
//   Corr    variant and nominal fully correlated (rho = 1, same events)
//   Full    rho from "correlations.txt". For R_JER this column equals Corr,
//           since no partial-rho propagation is implemented for the width.
enum class DRErr { Uncorr, Corr, Full };
inline DRErr& DRErrFlavor() {
    static DRErr f = DRErr::Uncorr;
    return f;
}
inline const char* DRErrTag(DRErr f) {
    return f == DRErr::Uncorr ? "uncorrelated"
                              : (f == DRErr::Corr ? "correlation_aware" : "full_correlation");
}
inline const char* DRErrLabel(DRErr f) {
    return f == DRErr::Uncorr ? "Uncorrelated errors"
                              : (f == DRErr::Corr ? "Correlation-aware errors" : "Full-correlation errors");
}

inline std::string QtyTag(Qty q) { return (q == Qty::JES) ? "R_JES" : "R_JER"; }

// Single Data/MC ratio axis title, e.g. R_{JES}(b_{JES}^{Data}/b_{JES}^{MC}).
inline std::string QtyAxis(Qty q) {
    return (q == Qty::JES) ? "R_{JES}(b_{JES}^{Data}/b_{JES}^{MC})" : "R_{JER}(b_{JER}^{Data}/b_{JER}^{MC})";
}

// Name shown on the plots for a correction. The variant keys are called
// "bJER10v.." for historical reasons but the calibration is the b-jet
// regression "bJR10v..", and "MuonCorr" is shown as "#mu-in-jet". Only the
// displayed text changes; keys, folder and branch names keep the original names.
inline std::string CorrDisp(std::string s) {
    std::string::size_type p = 0;
    while ((p = s.find("bJER", p)) != std::string::npos) {
        s.replace(p, 4, "bJR");
        p += 3;
    }
    p = 0;
    while ((p = s.find("MuonCorr", p)) != std::string::npos) {
        s.replace(p, 8, "#mu-in-jet");
        p += 10;
    }
    return s;
}

// Double-ratio axis title, e.g. R_{JES}^{bJR10v01}/R_{JES}^{Nominal}.
// corrLabel = the correction suffix for individual figures, "Corrected" when
// several corrections share one figure.
inline std::string DoubleRatioAxis(Qty q, const std::string& corrLabel) {
    const std::string Q = (q == Qty::JES) ? "R_{JES}" : "R_{JER}";
    return Q + "^{" + CorrDisp(corrLabel) + "}/" + Q + "^{Nominal}";
}

// Response axis title, e.g. b_{JES}^{Nominal} or b_{JER}^{bJR10v01}.
inline std::string RespAxis(Qty q, const std::string& tag) {
    return std::string(q == Qty::JES ? "b_{JES}" : "b_{JER}") + "^{" + CorrDisp(tag) + "}";
}

// Small white box with a text label, in NDC coordinates.
inline void DrawPanelTag(const std::string& text, double x1, double y1, double x2, double y2) {
    TPaveText* p = new TPaveText(x1, y1, x2, y2, "NDC");
    p->SetFillColor(kWhite);
    p->SetBorderSize(1);
    p->SetTextFont(42);
    p->SetTextSize(0.05);
    p->SetTextAlign(12);
    p->AddText(text.c_str());
    p->Draw();
}

inline double SingleNom(const Row& r, Qty q) { return q == Qty::JES ? r.rjes_nom : r.rjer_nom; }
inline double SingleNomErr(const Row& r, Qty q) { return q == Qty::JES ? r.rjes_nom_err : r.rjer_nom_err; }
inline double SingleVar(const Row& r, Qty q) { return q == Qty::JES ? r.rjes_var : r.rjer_var; }
inline double SingleVarErr(const Row& r, Qty q) { return q == Qty::JES ? r.rjes_var_err : r.rjer_var_err; }
inline double DoubleRatio(const Row& r, Qty q) { return q == Qty::JES ? r.rjes_dr : r.rjer_dr; }
inline double DoubleRatErr(const Row& r, Qty q) {
    const double unc = (q == Qty::JES ? r.rjes_dr_err : r.rjer_dr_err_unc);
    if (!r.hasFullErr) return unc; // CSV older than 44 columns: uncorrelated only
    switch (DRErrFlavor()) {
    case DRErr::Corr: {
        double v = (q == Qty::JES ? r.rjes_dr_err_corr : r.rjer_dr_err_corr);
        return v > 0 ? v : unc;
    }
    case DRErr::Full: {
        double v = (q == Qty::JES ? r.rjes_dr_err_full : r.rjer_dr_err_full);
        return v > 0 ? v : unc;
    }
    default:
        return unc;
    }
}

// The corrections compared to the nominal, and their colour on the plots.
inline std::vector<std::string> AllCorrectionSuffixes() { return {"bJER10v00", "bJER10v01", "MuonCorr"}; }

inline int CorrectionColor(const std::string& suffix) {
    if (suffix == "bJER10v00") return kRed + 1;
    if (suffix == "bJER10v01") return kBlue + 1;
    if (suffix == "MuonCorr") return kGreen + 2;
    return kMagenta + 1;
}

// Mass-cut tag used in file names; must match massCutDirNames[] in
// Response_Ratios_pTbins_Mcuts.h.
inline std::string MassTag(int idx) {
    static const char* t[3] = {"50-80GeV", "80-110GeV", "MassGE110GeV"};
    return (idx >= 0 && idx < 3) ? t[idx] : "massX";
}
constexpr int kNMass = 3;

// Rows for one selected correction, split by channel.
struct CorrData {
    std::string suffix;   // "bJER10v01"
    std::vector<Row> db;  // DB rows for this correction
    std::vector<Row> mjb; // MJB rows for this correction
};

// Graph of value +/- error vs pt_center for the rows of one mass cut, drawn
// with the given colour, marker and x offset. Returns nullptr if no row matches.
inline TGraphErrors* BuildGraph(const std::vector<Row>& rows, int massIdx,
                                const std::function<double(const Row&)>& gy,
                                const std::function<double(const Row&)>& ge, int color, int marker,
                                double xOffset = 0.0) {
    std::vector<double> x, y, ex, ey;
    for (const auto& r : rows) {
        if (MassBinIndex(r.mass_bin) != massIdx) continue;
        if (!r.hasSingles) continue; // need the single-ratio columns
        x.push_back(r.pt_center + xOffset);
        ex.push_back(0.0);
        y.push_back(gy(r));
        ey.push_back(ge(r));
    }
    if (x.empty()) return nullptr;
    auto* g = new TGraphErrors((int)x.size(), x.data(), y.data(), ex.data(), ey.data());
    g->SetMarkerColor(color);
    g->SetLineColor(color);
    g->SetMarkerStyle(marker);
    g->SetMarkerSize(1.1);
    g->SetLineWidth(2);
    return g;
}

// Expand [lo,hi] to include every point (+/- its error) of g.
inline void ExpandRange(TGraphErrors* g, double& lo, double& hi) {
    if (!g) return;
    for (int i = 0; i < g->GetN(); ++i) {
        double y = g->GetY()[i], e = g->GetEY()[i];
        lo = std::min(lo, y - e);
        hi = std::max(hi, y + e);
    }
}

// Optional fixed y range for a panel, e.g. to give Run 2 and Run 3 plots the
// same scale. EASYBJETS_YRANGE_<key>="lo,hi" replaces the automatic range;
// if unset, the automatic range is kept.
//   key = <quantity>_<panel>, quantity JES, JER, JMDP or JMDR (the last two
//   are used by Mass_Summary_Combined.h), panel SR (single ratio),
//   DR (double ratio) or COMP (response components: b, or mu / sigma in GeV).
// Example: export EASYBJETS_YRANGE_JES_DR="0.97,1.03"
inline void ApplySummaryYRange(const std::string& key, double& lo, double& hi) {
    const char* e = std::getenv(("EASYBJETS_YRANGE_" + key).c_str());
    if (!e || !*e) return;
    double a = 0.0, b = 0.0;
    char sep = 0;
    std::stringstream ss(e);
    if ((ss >> a >> sep >> b) && sep == ',' && b > a) {
        lo = a;
        hi = b;
    } else
        std::cerr << "[Summary] WARN: bad EASYBJETS_YRANGE_" << key << "='" << e << "' (expected lo,hi)\n";
}

// Same as BuildGraph, for the b_JES / b_JER components of the response.
inline TGraphErrors* BuildRespGraph(const std::vector<Row>& rows, int massIdx,
                                    const std::function<double(const Row&)>& gy,
                                    const std::function<double(const Row&)>& ge, int color, int marker,
                                    double xOffset = 0.0) {
    std::vector<double> x, y, ex, ey;
    for (const auto& r : rows) {
        if (MassBinIndex(r.mass_bin) != massIdx) continue;
        if (!r.hasResponse) continue; // need the b-component columns
        x.push_back(r.pt_center + xOffset);
        ex.push_back(0.0);
        y.push_back(gy(r));
        ey.push_back(ge(r));
    }
    if (x.empty()) return nullptr;
    auto* g = new TGraphErrors((int)x.size(), x.data(), y.data(), ex.data(), ey.data());
    g->SetMarkerColor(color);
    g->SetLineColor(color);
    g->SetMarkerStyle(marker);
    g->SetMarkerSize(1.1);
    g->SetLineWidth(2);
    return g;
}

// Marker per channel. Not used in this file at the moment
// (Response_Ratios_Combined.h defines its own copy).
constexpr int kMarkerDB = 20;  // closed circle
constexpr int kMarkerMJB = 22; // closed triangle

// Half-width of the non-closure band around 1 in the double-ratio panel.
//   R_JES: +/-1 %
//   R_JER: no band. Setting kDrawJERBand = true draws a band of half-width
//          RMS(DR - 1) of the MuonCorr double ratio, over one mass cut
//          (massIdx >= 0) or all of them (massIdx < 0).
// Returns -1 when no band is drawn.
inline double NonClosureHalfWidth(Qty q, int massIdx, const std::vector<CorrData>& corrs, bool wantDB,
                                  bool wantMJB) {
    if (q == Qty::JES) return 0.01;
    static const bool kDrawJERBand = false;
    if (!kDrawJERBand) return -1.0;
    const CorrData* mu = nullptr; // resolution -> muon-in-jet
    for (const auto& c : corrs)
        if (c.suffix.find("MuonCorr") != std::string::npos) {
            mu = &c;
            break;
        }
    if (!mu) return -1.0; // no muon series -> no band
    double s2 = 0.0;
    int n = 0;
    auto scan = [&](const std::vector<Row>& rows) {
        for (const auto& r : rows) {
            if (massIdx >= 0 && MassBinIndex(r.mass_bin) != massIdx) continue;
            if (!r.hasSingles) continue;
            double dr = DoubleRatio(r, Qty::JER);
            if (dr <= 0) continue;
            s2 += (dr - 1.0) * (dr - 1.0);
            ++n;
        }
    };
    if (wantDB) scan(mu->db);
    if (wantMJB) scan(mu->mjb);
    return n ? std::sqrt(s2 / n) : -1.0; // RMS deviation from unity
}

// One summary figure for one quantity, mass cut and topology:
//   top     data/MC single ratio, nominal (black circles) and corrections
//           (coloured triangles)
//   bottom  double ratio correction / nominal
// corrs:  the corrections to draw (at least one). perCorrection = true when
//         there is only one; the bottom axis is then named after it.
// topology: "DB", "MJB" or "combined".
inline void DrawSummaryFigure(Qty q, int massIdx, const std::string& wp_label, const std::string& topology,
                              const std::vector<CorrData>& corrs, const std::string& outPath,
                              bool perCorrection) {
    if (corrs.empty()) return;
    const bool wantDB = (topology == "DB" || topology == "combined");
    const bool wantMJB = (topology == "MJB" || topology == "combined");
    const bool combined = (topology == "combined");

    const double xlo = wantDB ? 200.0 : 450.0;
    const double xhi = wantMJB ? 1200.0 : 450.0;

    auto valNom = [&](const Row& r) {
        return SingleNom(r, q);
    };
    auto errNom = [&](const Row& r) {
        return SingleNomErr(r, q);
    };
    auto valVar = [&](const Row& r) {
        return SingleVar(r, q);
    };
    auto errVar = [&](const Row& r) {
        return SingleVarErr(r, q);
    };
    auto valDR = [&](const Row& r) {
        return DoubleRatio(r, q);
    };
    auto errDR = [&](const Row& r) {
        return DoubleRatErr(r, q);
    };

    // Build all graphs first, so the y ranges can be set from them.
    std::vector<TGraphErrors*> topGraphs, botGraphs;
    std::vector<std::pair<TGraphErrors*, std::string>> legTop;

    // Shift each series a little in pT so markers at the same pT do not
    // overlap. The shift scales with the width of the channel's pT range.
    const int nSeries = 1 + (int)corrs.size(); // nominal + corrections
    auto regW = [](bool isDB) {
        return isDB ? 250.0 : 750.0;
    };
    auto seriesOff = [&](int sidx, bool isDB) {
        return regW(isDB) * 0.03 * (sidx - (nSeries - 1) / 2.0);
    };

    // Nominal single ratio (black), taken from the first correction's rows
    // (every row also carries the nominal). One legend entry for both channels.
    const CorrData& first = corrs.front();
    bool nomLegendAdded = false;
    auto addNominal = [&](const std::vector<Row>& rows, int marker, bool isDB) {
        auto* gN = BuildGraph(rows, massIdx, valNom, errNom, kBlack, marker, seriesOff(0, isDB));
        if (gN) {
            gN->SetLineStyle(2); // dashed error bars for the nominal
            topGraphs.push_back(gN);
            if (!nomLegendAdded) {
                legTop.push_back({gN, "Nominal"});
                nomLegendAdded = true;
            }
        }
    };
    // Circle = nominal, triangle = correction. The channels separated apart by
    // their pT range (DB 200-450 GeV, MJB 450-1200 GeV).
    constexpr int kMkNom = 20, kMkCorr = 22;
    if (wantDB) addNominal(first.db, kMkNom, true);
    if (wantMJB) addNominal(first.mjb, kMkNom, false);

    // Each correction: single ratio (top) and double ratio (bottom).
    for (int ci = 0; ci < (int)corrs.size(); ++ci) {
        const auto& c = corrs[ci];
        const int col = CorrectionColor(c.suffix);
        const int sidx = ci + 1;
        TGraphErrors* legRef = nullptr;
        if (wantDB) {
            double off = seriesOff(sidx, true);
            auto* gV = BuildGraph(c.db, massIdx, valVar, errVar, col, kMkCorr, off);
            auto* gD = BuildGraph(c.db, massIdx, valDR, errDR, col, kMkCorr, off);
            if (gV) {
                topGraphs.push_back(gV);
                legRef = gV;
            }
            if (gD) botGraphs.push_back(gD);
        }
        if (wantMJB) {
            double off = seriesOff(sidx, false);
            auto* gV = BuildGraph(c.mjb, massIdx, valVar, errVar, col, kMkCorr, off);
            auto* gD = BuildGraph(c.mjb, massIdx, valDR, errDR, col, kMkCorr, off);
            if (gV) {
                topGraphs.push_back(gV);
                if (!legRef) legRef = gV;
            }
            if (gD) botGraphs.push_back(gD);
        }
        if (legRef) legTop.push_back({legRef, CorrDisp(c.suffix)});
    }

    if (topGraphs.empty()) {
        std::cerr << "[Summary] no data for " << QtyTag(q) << " mass " << MassTag(massIdx) << " (" << topology
                  << ") — skipping " << outPath << "\n";
        for (auto* g : botGraphs) delete g;
        return;
    }

    // y ranges from the points, with some margin.
    double tLo = 1e9, tHi = -1e9, bLo = 1e9, bHi = -1e9;
    for (auto* g : topGraphs) ExpandRange(g, tLo, tHi);
    for (auto* g : botGraphs) ExpandRange(g, bLo, bHi);
    if (tLo > tHi) {
        tLo = 0.9;
        tHi = 1.1;
    }
    if (bLo > bHi) {
        bLo = 0.9;
        bHi = 1.1;
    }
    auto pad = [](double& lo, double& hi) {
        double d = 0.15 * (hi - lo) + 1e-6;
        lo -= d;
        hi += d;
    };
    pad(tLo, tHi);
    pad(bLo, bHi);
    tHi += 0.35 * (tHi - tLo); // room for the legend at the top right
    {
        const std::string qt = (q == Qty::JES) ? "JES" : "JER"; // optional fixed y-range (env)
        ApplySummaryYRange(qt + "_SR", tLo, tHi);
        ApplySummaryYRange(qt + "_DR", bLo, bHi);
    }

    // Canvas: single ratio on top, double ratio below.
    TCanvas* c = new TCanvas(("c_" + outPath).c_str(), wp_label.c_str(), 800, 800);

    TPad* pad1 = new TPad("pad_top", "pad_top", 0.0, 0.32, 1.0, 1.0);
    pad1->SetBottomMargin(0.02);
    pad1->SetLeftMargin(0.14);
    pad1->SetTopMargin(0.06);
    pad1->SetRightMargin(0.05);
    pad1->Draw();
    TPad* pad2 = new TPad("pad_bot", "pad_bot", 0.0, 0.0, 1.0, 0.32);
    pad2->SetTopMargin(0.03);
    pad2->SetBottomMargin(0.34);
    pad2->SetLeftMargin(0.14);
    pad2->SetRightMargin(0.05);
    pad2->Draw();

    // Top panel
    pad1->cd();
    TH1F* frT = new TH1F(("frT_" + outPath).c_str(), "", 100, xlo, xhi);
    frT->SetStats(0);
    frT->SetMinimum(tLo);
    frT->SetMaximum(tHi);
    frT->GetYaxis()->SetTitle(QtyAxis(q).c_str());
    frT->GetYaxis()->SetTitleSize(0.055);
    frT->GetYaxis()->SetTitleOffset(1.15);
    frT->GetYaxis()->SetLabelSize(0.045);
    frT->GetXaxis()->SetLabelSize(0);
    frT->Draw("AXIS");

    TLine* refT = new TLine(xlo, 1.0, xhi, 1.0);
    refT->SetLineStyle(2);
    refT->SetLineColor(kGray + 2);
    refT->Draw("same");
    if (combined) {
        TLine* sep = new TLine(450.0, tLo, 450.0, tHi);
        sep->SetLineStyle(2);
        sep->SetLineColor(kGray + 2);
        sep->Draw("same");
    }
    for (auto* g : topGraphs) g->Draw("P same");

    // In combined mode, write the channel name under each pT range.
    if (combined) {
        TLatex reg;
        reg.SetTextAlign(22);
        reg.SetTextSize(0.036);
        reg.SetTextColor(kGray + 2);
        reg.DrawLatex(0.5 * (xlo + 450.0), tLo + 0.05 * (tHi - tLo), "Z(#rightarrow b#bar{b}) + #gamma");
        reg.DrawLatex(0.5 * (450.0 + xhi), tLo + 0.05 * (tHi - tLo), "Z(#rightarrow b#bar{b}) + jets");
    }

    // ATLAS label, lumi, mass cut, WP
    TLatex tex;
    tex.SetNDC();
    tex.SetTextFont(42);
    tex.SetTextSize(0.052);
    tex.DrawLatex(0.17, 0.86, "#font[72]{ATLAS} Preliminary");
    tex.SetTextSize(0.040);
    tex.DrawLatex(0.17, 0.79, RunConfig::LumiLabel().c_str());
    tex.DrawLatex(0.17, 0.73, MassBinDisplayLabel(massIdx).c_str());
    tex.DrawLatex(0.17, 0.67, ("WP: " + wp_label).c_str());

    // Legend at the top right, one line per series.
    double legH = 0.06 * (double)legTop.size();
    TLegend* legA = new TLegend(0.77, 0.88 - 0.06 * legTop.size(), 0.95, 0.88);
    legA->SetBorderSize(0);
    legA->SetFillStyle(0);
    legA->SetFillColor(kWhite);
    legA->SetTextSize(0.045);
    for (auto& e : legTop) legA->AddEntry(e.first, e.second.c_str(), "lep");
    legA->Draw();

    // Bottom panel: double ratio
    pad2->cd();
    TH1F* frB = new TH1F(("frB_" + outPath).c_str(), "", 100, xlo, xhi);
    frB->SetStats(0);
    frB->SetMinimum(bLo);
    frB->SetMaximum(bHi);
    std::string drTitle = DoubleRatioAxis(q, perCorrection ? corrs.front().suffix : "Corrected");
    frB->GetYaxis()->SetTitle(drTitle.c_str());
    frB->GetYaxis()->SetNdivisions(506);
    frB->GetYaxis()->SetTitleSize(0.11);
    frB->GetYaxis()->SetTitleOffset(0.55);
    frB->GetYaxis()->SetLabelSize(0.095);
    frB->GetXaxis()->SetTitle("Large-#it{R} (#it{b}#bar{#it{b}})-jet #it{p}_{T} [GeV]");
    frB->GetXaxis()->SetTitleSize(0.12);
    frB->GetXaxis()->SetTitleOffset(1.15);
    frB->GetXaxis()->SetLabelSize(0.095);
    frB->Draw("AXIS");

    // Non-closure band around 1 (see NonClosureHalfWidth: R_JES only).
    double ncHalf = NonClosureHalfWidth(q, massIdx, corrs, wantDB, wantMJB);
    if (ncHalf > 0) {
        TBox* ncB = new TBox(xlo, 1.0 - ncHalf, xhi, 1.0 + ncHalf);
        ncB->SetFillColorAlpha(kGray + 1, 0.35);
        ncB->SetLineWidth(0);
        ncB->Draw("same");
        // Labelled as an assigned uncertainty
        TLegend* legBand = new TLegend(0.5, 0.80, 0.95, 0.97);
        legBand->SetBorderSize(0);
        legBand->SetFillStyle(0);
        legBand->SetTextSize(0.72 * frB->GetXaxis()->GetTitleSize());
        legBand->AddEntry(ncB, Form("Assigned non-closure unc. (#pm%g%%)", ncHalf * 100.0), "f");
        legBand->Draw();
    }
    TLine* refB = new TLine(xlo, 1.0, xhi, 1.0);
    refB->SetLineStyle(2);
    refB->SetLineColor(kGray + 2);
    refB->Draw("same");
    if (combined) {
        TLine* sepB = new TLine(450.0, bLo, 450.0, bHi);
        sepB->SetLineStyle(2);
        sepB->SetLineColor(kGray + 2);
        sepB->Draw("same");
    }
    for (auto* g : botGraphs) g->Draw("P same");

    c->SaveAs(outPath.c_str());
    std::cout << "[Summary] Saved: " << outPath << "\n";

    for (auto* g : topGraphs) delete g;
    for (auto* g : botGraphs) delete g;
    delete c;
}

// Same summary figure for one correction, with the three mass cuts overlaid
// Top: single ratio, bottom: double ratio.
inline void DrawSummaryFigureMassOverlay(Qty q, const std::string& wp_label, const std::string& topology,
                                         const CorrData& corr, const std::string& outPath) {
    const bool wantDB = (topology == "DB" || topology == "combined");
    const bool wantMJB = (topology == "MJB" || topology == "combined");
    const bool combined = (topology == "combined");
    const double xlo = wantDB ? 200.0 : 450.0;
    const double xhi = wantMJB ? 1200.0 : 450.0;

    auto valNom = [&](const Row& r) {
        return SingleNom(r, q);
    };
    auto errNom = [&](const Row& r) {
        return SingleNomErr(r, q);
    };
    auto valVar = [&](const Row& r) {
        return SingleVar(r, q);
    };
    auto errVar = [&](const Row& r) {
        return SingleVarErr(r, q);
    };
    auto valDR = [&](const Row& r) {
        return DoubleRatio(r, q);
    };
    auto errDR = [&](const Row& r) {
        return DoubleRatErr(r, q);
    };

    auto regW = [](bool isDB) {
        return isDB ? 250.0 : 750.0;
    };
    // Shift in pT: -3 %, 0, +3 % of the channel's pT range per mass cut, plus a
    // small extra shift between nominal and correction.
    auto massOff = [&](int mi, bool isDB, bool isVar) {
        static const double f[3] = {-0.03, 0.0, 0.03};
        return regW(isDB) * f[std::min(std::max(mi, 0), 2)] + (isVar ? +0.008 : -0.008) * regW(isDB);
    };
    // Circle = nominal, triangle = correction, in both channels
    auto markerFor = [&](bool /*isDB*/, bool isVar) {
        return isVar ? 22 : 20; // triangle = correction, circle = nominal
    };

    std::vector<TGraphErrors*> topGraphs, botGraphs;
    auto addChannel = [&](const std::vector<Row>& rows, bool isDB) {
        for (int mi = 0; mi < kNMass; ++mi) {
            int col = MassCutColor(mi);
            auto* gN =
                BuildGraph(rows, mi, valNom, errNom, col, markerFor(isDB, false), massOff(mi, isDB, false));
            auto* gV =
                BuildGraph(rows, mi, valVar, errVar, col, markerFor(isDB, true), massOff(mi, isDB, true));
            auto* gD =
                BuildGraph(rows, mi, valDR, errDR, col, markerFor(isDB, true), massOff(mi, isDB, true));
            if (gN) {
                gN->SetLineStyle(2);
                topGraphs.push_back(gN);
            } // dashed: nominal error bars
            if (gV) topGraphs.push_back(gV);
            if (gD) botGraphs.push_back(gD);
        }
    };
    if (wantDB) addChannel(corr.db, true);
    if (wantMJB) addChannel(corr.mjb, false);

    if (topGraphs.empty()) {
        std::cerr << "[Summary] no data for " << QtyTag(q) << " AllMass " << corr.suffix << " (" << topology
                  << ") — skipping " << outPath << "\n";
        for (auto* g : botGraphs) delete g;
        return;
    }

    // y ranges from the plotted graphs only (the legend markers below are
    // placed off-screen and must not count).
    double tLo = 1e9, tHi = -1e9, bLo = 1e9, bHi = -1e9;
    for (auto* g : topGraphs) ExpandRange(g, tLo, tHi);
    for (auto* g : botGraphs) ExpandRange(g, bLo, bHi);
    if (tLo > tHi) {
        tLo = 0.9;
        tHi = 1.1;
    }
    if (bLo > bHi) {
        bLo = 0.9;
        bHi = 1.1;
    }
    auto pad = [](double& lo, double& hi) {
        double d = 0.15 * (hi - lo) + 1e-6;
        lo -= d;
        hi += d;
    };
    pad(tLo, tHi);
    pad(bLo, bHi);
    tHi += 0.35 * (tHi - tLo); // room for the legend at the top right
    {
        const std::string qt = (q == Qty::JES) ? "JES" : "JER"; // optional fixed y-range (env)
        ApplySummaryYRange(qt + "_SR", tLo, tHi);
        ApplySummaryYRange(qt + "_DR", bLo, bHi);
    }

    // Legend markers, one colour per mass cut. Placed at (-1e9, -1e9) so they
    // only appear in the legend, never on the plot.
    std::vector<TGraphErrors*> legendDummies;
    std::vector<std::pair<TGraphErrors*, std::string>> legMass;
    for (int mi = 0; mi < kNMass; ++mi) {
        auto* gl = new TGraphErrors(1);
        gl->SetPoint(0, -1e9, -1e9);
        gl->SetMarkerColor(MassCutColor(mi));
        gl->SetLineColor(MassCutColor(mi));
        gl->SetMarkerStyle(20);
        legendDummies.push_back(gl);
        legMass.push_back({gl, MassBinDisplayLabel(mi)});
    }

    TCanvas* c = new TCanvas(("c_" + outPath).c_str(), wp_label.c_str(), 800, 800);
    TPad* pad1 = new TPad("pad_top", "pad_top", 0.0, 0.32, 1.0, 1.0);
    pad1->SetBottomMargin(0.02);
    pad1->SetLeftMargin(0.14);
    pad1->SetTopMargin(0.06);
    pad1->SetRightMargin(0.05);
    pad1->Draw();
    TPad* pad2 = new TPad("pad_bot", "pad_bot", 0.0, 0.0, 1.0, 0.32);
    pad2->SetTopMargin(0.03);
    pad2->SetBottomMargin(0.34);
    pad2->SetLeftMargin(0.14);
    pad2->SetRightMargin(0.05);
    pad2->Draw();

    // Top panel: single ratio
    pad1->cd();
    TH1F* frT = new TH1F(("frTo_" + outPath).c_str(), "", 100, xlo, xhi);
    frT->SetStats(0);
    frT->SetMinimum(tLo);
    frT->SetMaximum(tHi);
    frT->GetYaxis()->SetTitle(QtyAxis(q).c_str());
    frT->GetYaxis()->SetTitleSize(0.055);
    frT->GetYaxis()->SetTitleOffset(1.15);
    frT->GetYaxis()->SetLabelSize(0.045);
    frT->GetXaxis()->SetLabelSize(0);
    frT->Draw("AXIS");
    TLine* refT = new TLine(xlo, 1.0, xhi, 1.0);
    refT->SetLineStyle(2);
    refT->SetLineColor(kGray + 2);
    refT->Draw("same");
    if (combined) {
        TLine* s = new TLine(450.0, tLo, 450.0, tHi);
        s->SetLineStyle(2);
        s->SetLineColor(kGray + 2);
        s->Draw("same");
    }
    for (auto* g : topGraphs) g->Draw("P same");
    if (combined) {
        TLatex reg;
        reg.SetTextAlign(22);
        reg.SetTextSize(0.036);
        reg.SetTextColor(kGray + 2);
        reg.DrawLatex(0.5 * (xlo + 450.0), tLo + 0.05 * (tHi - tLo), "Z(#rightarrow b#bar{b}) + #gamma");
        reg.DrawLatex(0.5 * (450.0 + xhi), tLo + 0.05 * (tHi - tLo), "Z(#rightarrow b#bar{b}) + jets");
    }
    TLatex tex;
    tex.SetNDC();
    tex.SetTextFont(42);
    tex.SetTextSize(0.052);
    tex.DrawLatex(0.17, 0.86, "#font[72]{ATLAS} Preliminary");
    tex.SetTextSize(0.038);
    tex.DrawLatex(0.17, 0.79, RunConfig::LumiLabel().c_str());
    tex.DrawLatex(0.17, 0.74, ("WP: " + wp_label).c_str());

    // One legend in two columns: left the mass cuts (colours), right the
    // nominal (circle, dashed) and correction (triangle, solid) markers.
    // "lp" shows the error bar style.
    auto* gNomKey = new TGraphErrors(1);
    gNomKey->SetPoint(0, -1e9, -1e9);
    gNomKey->SetMarkerStyle(20);
    gNomKey->SetMarkerColor(kBlack);
    gNomKey->SetLineColor(kBlack);
    gNomKey->SetLineStyle(2); // dashed error bars
    auto* gCorrKey = new TGraphErrors(1);
    gCorrKey->SetPoint(0, -1e9, -1e9);
    gCorrKey->SetMarkerStyle(22);
    gCorrKey->SetMarkerColor(kBlack);
    gCorrKey->SetLineColor(kBlack);
    gCorrKey->SetLineStyle(1); // solid error bars
    legendDummies.push_back(gNomKey);
    legendDummies.push_back(gCorrKey);

    TLegend* leg = new TLegend(0.46, 0.66, 0.93, 0.89);
    leg->SetBorderSize(1);
    leg->SetFillStyle(1001);
    leg->SetFillColor(kWhite);
    leg->SetTextSize(0.034);
    leg->SetNColumns(2);

    std::vector<std::pair<TObject*, std::string>> leftCol, rightCol;
    for (auto& e : legMass) leftCol.push_back({e.first, e.second});
    rightCol.push_back({gNomKey, "Nominal"});
    rightCol.push_back({gCorrKey, CorrDisp(corr.suffix)});
    size_t nrow = std::max(leftCol.size(), rightCol.size());
    for (size_t r = 0; r < nrow; ++r) {
        if (r < leftCol.size())
            leg->AddEntry(leftCol[r].first, leftCol[r].second.c_str(), "p");
        else
            leg->AddEntry((TObject*)nullptr, "", "");
        if (r < rightCol.size())
            leg->AddEntry(rightCol[r].first, rightCol[r].second.c_str(), "lp");
        else
            leg->AddEntry((TObject*)nullptr, "", "");
    }
    leg->Draw();

    // Bottom panel: double ratio
    pad2->cd();
    TH1F* frB = new TH1F(("frBo_" + outPath).c_str(), "", 100, xlo, xhi);
    frB->SetStats(0);
    frB->SetMinimum(bLo);
    frB->SetMaximum(bHi);
    frB->GetYaxis()->SetTitle(DoubleRatioAxis(q, corr.suffix).c_str());
    frB->GetYaxis()->SetNdivisions(506);
    frB->GetYaxis()->SetTitleSize(0.11);
    frB->GetYaxis()->SetTitleOffset(0.55);
    frB->GetYaxis()->SetLabelSize(0.095);
    frB->GetXaxis()->SetTitle("Large-#it{R} (#it{b}#bar{#it{b}})-jet #it{p}_{T} [GeV]");
    frB->GetXaxis()->SetTitleSize(0.12);
    frB->GetXaxis()->SetTitleOffset(1.15);
    frB->GetXaxis()->SetLabelSize(0.095);
    frB->Draw("AXIS");
    // Non-closure band around 1 (R_JES only, see NonClosureHalfWidth), over
    // all mass cuts (massIdx = -1).
    double ncHalf = NonClosureHalfWidth(q, -1, std::vector<CorrData>{corr}, wantDB, wantMJB);
    if (ncHalf > 0) {
        TBox* ncB = new TBox(xlo, 1.0 - ncHalf, xhi, 1.0 + ncHalf);
        ncB->SetFillColorAlpha(kGray + 1, 0.35);
        ncB->SetLineWidth(0);
        ncB->Draw("same");
        TLegend* legBand = new TLegend(0.5, 0.80, 0.95, 0.97);
        legBand->SetBorderSize(0);
        legBand->SetFillStyle(0);
        legBand->SetTextSize(0.72 * frB->GetXaxis()->GetTitleSize());
        legBand->AddEntry(ncB, Form("Assigned non-closure unc. (#pm%g%%)", ncHalf * 100.0), "f");
        legBand->Draw();
    }
    TLine* refB = new TLine(xlo, 1.0, xhi, 1.0);
    refB->SetLineStyle(2);
    refB->SetLineColor(kGray + 2);
    refB->Draw("same");
    if (combined) {
        TLine* s = new TLine(450.0, bLo, 450.0, bHi);
        s->SetLineStyle(2);
        s->SetLineColor(kGray + 2);
        s->Draw("same");
    }
    for (auto* g : botGraphs) g->Draw("P same");

    c->SaveAs(outPath.c_str());
    std::cout << "[Summary] Saved: " << outPath << "\n";

    for (auto* g : topGraphs) delete g;
    for (auto* g : botGraphs) delete g;
    for (auto* g : legendDummies) delete g;
    delete c;
}

// Response figure for one quantity, mass cut and correction: the fitted
// b^Data and b^MC (b_JES or b_JER) vs pT, nominal calibration on top and
// correction below. Both panels have the same height and y range so they can
// be compared directly. Data = black circles (dashed error bars), MC = red
// triangles.
inline void DrawResponseFigure(Qty q, int massIdx, const std::string& wp_label, const std::string& topology,
                               const CorrData& corr, const std::string& outPath) {
    const bool wantDB = (topology == "DB" || topology == "combined");
    const bool wantMJB = (topology == "MJB" || topology == "combined");
    const bool combined = (topology == "combined");
    const double xlo = wantDB ? 200.0 : 450.0;
    const double xhi = wantMJB ? 1200.0 : 450.0;

    auto dataNom = [&](const Row& r) {
        return q == Qty::JES ? r.bjes_data_nom : r.bjer_data_nom;
    };
    auto dataNomE = [&](const Row& r) {
        return q == Qty::JES ? r.bjes_data_nom_err : r.bjer_data_nom_err;
    };
    auto mcNom = [&](const Row& r) {
        return q == Qty::JES ? r.bjes_mc_nom : r.bjer_mc_nom;
    };
    auto mcNomE = [&](const Row& r) {
        return q == Qty::JES ? r.bjes_mc_nom_err : r.bjer_mc_nom_err;
    };
    auto dataVar = [&](const Row& r) {
        return q == Qty::JES ? r.bjes_data_var : r.bjer_data_var;
    };
    auto dataVarE = [&](const Row& r) {
        return q == Qty::JES ? r.bjes_data_var_err : r.bjer_data_var_err;
    };
    auto mcVar = [&](const Row& r) {
        return q == Qty::JES ? r.bjes_mc_var : r.bjer_mc_var;
    };
    auto mcVarE = [&](const Row& r) {
        return q == Qty::JES ? r.bjes_mc_var_err : r.bjer_mc_var_err;
    };

    auto regW = [](bool isDB) {
        return isDB ? 250.0 : 750.0;
    };
    auto off = [&](bool isDB, bool isData) {
        return regW(isDB) * (isData ? -0.012 : +0.012);
    };
    const int kData = kBlack, kMC = kRed + 1;
    const int mData = 20, mMC = 22; // circle (Data, dashed bars) / triangle (MC, solid)

    std::vector<TGraphErrors*> topGraphs, botGraphs;
    TGraphErrors *legD = nullptr, *legM = nullptr;
    auto fillChannel = [&](const std::vector<Row>& rows, bool isDB) {
        auto* dN = BuildRespGraph(rows, massIdx, dataNom, dataNomE, kData, mData, off(isDB, true));
        auto* mN = BuildRespGraph(rows, massIdx, mcNom, mcNomE, kMC, mMC, off(isDB, false));
        auto* dV = BuildRespGraph(rows, massIdx, dataVar, dataVarE, kData, mData, off(isDB, true));
        auto* mV = BuildRespGraph(rows, massIdx, mcVar, mcVarE, kMC, mMC, off(isDB, false));
        if (dN) dN->SetLineStyle(2);
        if (dV) dV->SetLineStyle(2); // Data: dashed error bars
        if (dN) {
            topGraphs.push_back(dN);
            if (!legD) legD = dN;
        }
        if (mN) {
            topGraphs.push_back(mN);
            if (!legM) legM = mN;
        }
        if (dV) botGraphs.push_back(dV);
        if (mV) botGraphs.push_back(mV);
    };
    if (wantDB) fillChannel(corr.db, true);
    if (wantMJB) fillChannel(corr.mjb, false);

    if (topGraphs.empty() && botGraphs.empty()) {
        std::cerr << "[Summary] no b_" << (q == Qty::JES ? "JES" : "JER") << " data for " << corr.suffix
                  << " mass " << MassTag(massIdx) << " (" << topology << ") — skipping " << outPath << "\n";
        return;
    }

    // Same y range in both panels, plus room for the legend in the top one.
    double yLo = 1e9, yHi = -1e9;
    for (auto* g : topGraphs) ExpandRange(g, yLo, yHi);
    for (auto* g : botGraphs) ExpandRange(g, yLo, yHi);
    if (yLo > yHi) {
        yLo = 0.0;
        yHi = 1.0;
    }
    double dpad = 0.15 * (yHi - yLo) + 1e-6;
    yLo -= dpad;
    yHi += dpad;
    double yHiTop = yHi + 0.42 * (yHi - yLo);
    {
        const std::string qt = (q == Qty::JES) ? "JES" : "JER"; // optional fixed y-range (env)
        ApplySummaryYRange(qt + "_COMP", yLo, yHi);
        yHiTop = yHi + 0.42 * (yHi - yLo);
    }

    TCanvas* c = new TCanvas(("c_" + outPath).c_str(), wp_label.c_str(), 800, 800);
    TPad* pad1 = new TPad("pad_nom", "pad_nom", 0.0, 0.50, 1.0, 1.0);
    pad1->SetBottomMargin(0.02);
    pad1->SetTopMargin(0.08);
    pad1->SetLeftMargin(0.14);
    pad1->SetRightMargin(0.05);
    pad1->Draw();
    TPad* pad2 = new TPad("pad_cor", "pad_cor", 0.0, 0.0, 1.0, 0.50);
    pad2->SetTopMargin(0.02);
    pad2->SetBottomMargin(0.20);
    pad2->SetLeftMargin(0.14);
    pad2->SetRightMargin(0.05);
    pad2->Draw();

    auto drawDivider = [&](double lo, double hi) {
        if (!combined) return;
        TLine* s = new TLine(450.0, lo, 450.0, hi);
        s->SetLineStyle(2);
        s->SetLineColor(kGray + 2);
        s->Draw("same");
    };

    // Top panel: nominal calibration
    pad1->cd();
    TH1F* frT = new TH1F(("frRT_" + outPath).c_str(), "", 100, xlo, xhi);
    frT->SetStats(0);
    frT->SetMinimum(yLo);
    frT->SetMaximum(yHiTop);
    frT->GetYaxis()->SetTitle(RespAxis(q, "Nominal").c_str());
    frT->GetYaxis()->SetTitleSize(0.07);
    frT->GetYaxis()->SetTitleOffset(0.95);
    frT->GetYaxis()->SetLabelSize(0.055);
    frT->GetXaxis()->SetLabelSize(0);
    frT->Draw("AXIS");
    drawDivider(yLo, yHiTop);
    for (auto* g : topGraphs) g->Draw("P same");

    TLatex tex;
    tex.SetNDC();
    tex.SetTextFont(42);
    tex.SetTextSize(0.07);
    tex.DrawLatex(0.17, 0.85, "#font[72]{ATLAS} #scale[0.79]{Preliminary}");
    tex.SetTextSize(0.055);
    tex.SetTextSize(0.045);
    tex.DrawLatex(0.17, 0.77, ("WP: " + wp_label).c_str());
    tex.DrawLatex(0.17, 0.70, MassBinDisplayLabel(massIdx).c_str());
    if (combined) {
        TLatex reg;
        reg.SetTextAlign(22);
        reg.SetTextSize(0.05);
        reg.SetTextColor(kGray + 2);
        reg.DrawLatex(0.5 * (xlo + 450.0), yLo + 0.06 * (yHiTop - yLo), "Z(#rightarrow b#bar{b}) + #gamma");
        reg.DrawLatex(0.5 * (450.0 + xhi), yLo + 0.06 * (yHiTop - yLo), "Z(#rightarrow b#bar{b}) + jets");
    }
    // Panel name at the top right (no box, so it never hides error bars),
    // with the Data/MC legend below it.
    TLatex tag;
    tag.SetNDC();
    tag.SetTextFont(42);
    tag.SetTextSize(0.055);
    tag.DrawLatex(0.70, 0.85, "Nominal");
    TLegend* leg = new TLegend(0.70, 0.62, 0.92, 0.80);
    leg->SetBorderSize(1);
    leg->SetFillStyle(1001);
    leg->SetFillColor(kWhite);
    leg->SetTextSize(0.05);
    if (legD) leg->AddEntry(legD, "Data", "lep");
    if (legM) leg->AddEntry(legM, "MC", "lep");
    leg->Draw();

    // Bottom panel: correction
    pad2->cd();
    TH1F* frB = new TH1F(("frRB_" + outPath).c_str(), "", 100, xlo, xhi);
    frB->SetStats(0);
    frB->SetMinimum(yLo);
    frB->SetMaximum(yHi);
    frB->GetYaxis()->SetTitle(RespAxis(q, corr.suffix).c_str());
    frB->GetYaxis()->SetTitleSize(0.07);
    frB->GetYaxis()->SetTitleOffset(0.95);
    frB->GetYaxis()->SetLabelSize(0.055);
    frB->GetXaxis()->SetTitle("Large-#it{R} (#it{b}#bar{#it{b}})-jet #it{p}_{T} [GeV]");
    frB->GetXaxis()->SetTitleSize(0.065);
    frB->GetXaxis()->SetTitleOffset(1.25);
    frB->GetXaxis()->SetLabelSize(0.055);
    frB->Draw("AXIS");
    drawDivider(yLo, yHi);
    for (auto* g : botGraphs) g->Draw("P same");
    TLatex tagB;
    tagB.SetNDC();
    tagB.SetTextFont(42);
    tagB.SetTextSize(0.055);
    tagB.DrawLatex(0.70, 0.88, CorrDisp(corr.suffix).c_str());

    c->SaveAs(outPath.c_str());
    std::cout << "[Summary] Saved: " << outPath << "\n";

    for (auto* g : topGraphs) delete g;
    for (auto* g : botGraphs) delete g;
    delete c;
}

// Same response figure for one correction, with the three mass cuts overlaid
inline void DrawResponseFigureMassOverlay(Qty q, const std::string& wp_label, const std::string& topology,
                                          const CorrData& corr, const std::string& outPath) {
    const bool wantDB = (topology == "DB" || topology == "combined");
    const bool wantMJB = (topology == "MJB" || topology == "combined");
    const bool combined = (topology == "combined");
    const double xlo = wantDB ? 200.0 : 450.0;
    const double xhi = wantMJB ? 1200.0 : 450.0;

    auto dataNom = [&](const Row& r) {
        return q == Qty::JES ? r.bjes_data_nom : r.bjer_data_nom;
    };
    auto dataNomE = [&](const Row& r) {
        return q == Qty::JES ? r.bjes_data_nom_err : r.bjer_data_nom_err;
    };
    auto mcNom = [&](const Row& r) {
        return q == Qty::JES ? r.bjes_mc_nom : r.bjer_mc_nom;
    };
    auto mcNomE = [&](const Row& r) {
        return q == Qty::JES ? r.bjes_mc_nom_err : r.bjer_mc_nom_err;
    };
    auto dataVar = [&](const Row& r) {
        return q == Qty::JES ? r.bjes_data_var : r.bjer_data_var;
    };
    auto dataVarE = [&](const Row& r) {
        return q == Qty::JES ? r.bjes_data_var_err : r.bjer_data_var_err;
    };
    auto mcVar = [&](const Row& r) {
        return q == Qty::JES ? r.bjes_mc_var : r.bjer_mc_var;
    };
    auto mcVarE = [&](const Row& r) {
        return q == Qty::JES ? r.bjes_mc_var_err : r.bjer_mc_var_err;
    };

    auto regW = [](bool isDB) {
        return isDB ? 250.0 : 750.0;
    };
    auto roff = [&](int mi, bool isDB, bool isData) {
        static const double f[3] = {-0.03, 0.0, 0.03};
        return regW(isDB) * f[std::min(std::max(mi, 0), 2)] + (isData ? -0.008 : +0.008) * regW(isDB);
    };
    const int mData = 20, mMC = 22; // circle (Data, dashed bars) / triangle (MC, solid)

    std::vector<TGraphErrors*> topGraphs, botGraphs;
    auto addChannel = [&](const std::vector<Row>& rows, bool isDB) {
        for (int mi = 0; mi < kNMass; ++mi) {
            int col = MassCutColor(mi);
            auto* dN = BuildRespGraph(rows, mi, dataNom, dataNomE, col, mData, roff(mi, isDB, true));
            auto* mN = BuildRespGraph(rows, mi, mcNom, mcNomE, col, mMC, roff(mi, isDB, false));
            auto* dV = BuildRespGraph(rows, mi, dataVar, dataVarE, col, mData, roff(mi, isDB, true));
            auto* mV = BuildRespGraph(rows, mi, mcVar, mcVarE, col, mMC, roff(mi, isDB, false));
            if (dN) {
                dN->SetLineStyle(2);
                topGraphs.push_back(dN);
            } // Data: dashed bars
            if (mN) topGraphs.push_back(mN);
            if (dV) {
                dV->SetLineStyle(2);
                botGraphs.push_back(dV);
            }
            if (mV) botGraphs.push_back(mV);
        }
    };
    if (wantDB) addChannel(corr.db, true);
    if (wantMJB) addChannel(corr.mjb, false);

    if (topGraphs.empty() && botGraphs.empty()) {
        std::cerr << "[Summary] no b_" << (q == Qty::JES ? "JES" : "JER") << " data for " << corr.suffix
                  << " AllMass (" << topology << ") — skipping " << outPath << "\n";
        return;
    }

    double yLo = 1e9, yHi = -1e9;
    for (auto* g : topGraphs) ExpandRange(g, yLo, yHi);
    for (auto* g : botGraphs) ExpandRange(g, yLo, yHi);
    if (yLo > yHi) {
        yLo = 0.0;
        yHi = 1.0;
    }
    double dpad = 0.15 * (yHi - yLo) + 1e-6;
    yLo -= dpad;
    yHi += dpad;
    double yHiTop = yHi + 0.45 * (yHi - yLo);
    {
        const std::string qt = (q == Qty::JES) ? "JES" : "JER"; // optional fixed y-range (env)
        ApplySummaryYRange(qt + "_COMP", yLo, yHi);
        yHiTop = yHi + 0.45 * (yHi - yLo);
    }

    // Legend markers (off-screen): mass-cut colours on the left, Data / MC
    // markers on the right.
    std::vector<TGraphErrors*> legendDummies;
    std::vector<std::pair<TObject*, std::string>> leftCol, rightCol;
    for (int mi = 0; mi < kNMass; ++mi) {
        auto* gl = new TGraphErrors(1);
        gl->SetPoint(0, -1e9, -1e9);
        gl->SetMarkerColor(MassCutColor(mi));
        gl->SetLineColor(MassCutColor(mi));
        gl->SetMarkerStyle(20);
        legendDummies.push_back(gl);
        leftCol.push_back({gl, MassBinDisplayLabel(mi)});
    }
    auto* gData = new TGraphErrors(1);
    gData->SetPoint(0, -1e9, -1e9);
    gData->SetMarkerStyle(mData);
    gData->SetMarkerColor(kBlack);
    gData->SetLineColor(kBlack);
    gData->SetLineStyle(2); // Data: dashed
    auto* gMC = new TGraphErrors(1);
    gMC->SetPoint(0, -1e9, -1e9);
    gMC->SetMarkerStyle(mMC);
    gMC->SetMarkerColor(kBlack);
    gMC->SetLineColor(kBlack);
    legendDummies.push_back(gData);
    legendDummies.push_back(gMC);
    rightCol.push_back({gData, "Data"});
    rightCol.push_back({gMC, "MC"});

    TCanvas* c = new TCanvas(("c_" + outPath).c_str(), wp_label.c_str(), 800, 800);
    TPad* pad1 = new TPad("pad_nom", "pad_nom", 0.0, 0.50, 1.0, 1.0);
    pad1->SetBottomMargin(0.02);
    pad1->SetTopMargin(0.08);
    pad1->SetLeftMargin(0.14);
    pad1->SetRightMargin(0.05);
    pad1->Draw();
    TPad* pad2 = new TPad("pad_cor", "pad_cor", 0.0, 0.0, 1.0, 0.50);
    pad2->SetTopMargin(0.02);
    pad2->SetBottomMargin(0.20);
    pad2->SetLeftMargin(0.14);
    pad2->SetRightMargin(0.05);
    pad2->Draw();

    auto drawDivider = [&](double lo, double hi) {
        if (!combined) return;
        TLine* s = new TLine(450.0, lo, 450.0, hi);
        s->SetLineStyle(2);
        s->SetLineColor(kGray + 2);
        s->Draw("same");
    };

    // Top panel: nominal calibration
    pad1->cd();
    TH1F* frT = new TH1F(("frROT_" + outPath).c_str(), "", 100, xlo, xhi);
    frT->SetStats(0);
    frT->SetMinimum(yLo);
    frT->SetMaximum(yHiTop);
    frT->GetYaxis()->SetTitle(RespAxis(q, "Nominal").c_str());
    frT->GetYaxis()->SetTitleSize(0.07);
    frT->GetYaxis()->SetTitleOffset(0.95);
    frT->GetYaxis()->SetLabelSize(0.055);
    frT->GetXaxis()->SetLabelSize(0);
    frT->Draw("AXIS");
    drawDivider(yLo, yHiTop);
    for (auto* g : topGraphs) g->Draw("P same");

    TLatex tex;
    tex.SetNDC();
    tex.SetTextFont(42);
    tex.SetTextSize(0.07);
    tex.DrawLatex(0.17, 0.85, "#font[72]{ATLAS} #scale[0.79]{Preliminary}");
    tex.SetTextSize(0.055);
    tex.SetTextSize(0.045);
    tex.DrawLatex(0.17, 0.77, ("WP: " + wp_label).c_str());
    if (combined) {
        TLatex reg;
        reg.SetTextAlign(22);
        reg.SetTextSize(0.05);
        reg.SetTextColor(kGray + 2);
        reg.DrawLatex(0.5 * (xlo + 450.0), yLo + 0.06 * (yHiTop - yLo), "Z(#rightarrow b#bar{b}) + #gamma");
        reg.DrawLatex(0.5 * (450.0 + xhi), yLo + 0.06 * (yHiTop - yLo), "Z(#rightarrow b#bar{b}) + jets");
    }
    // Panel name at the top right, two-column legend below it.
    TLatex tag;
    tag.SetNDC();
    tag.SetTextFont(42);
    tag.SetTextSize(0.055);
    tag.DrawLatex(0.70, 0.85, "Nominal");
    TLegend* leg = new TLegend(0.52, 0.46, 0.92, 0.76);
    leg->SetBorderSize(1);
    leg->SetFillStyle(1001);
    leg->SetFillColor(kWhite);
    leg->SetTextSize(0.04);
    leg->SetNColumns(2);
    size_t nrow = std::max(leftCol.size(), rightCol.size());
    for (size_t r = 0; r < nrow; ++r) {
        if (r < leftCol.size())
            leg->AddEntry(leftCol[r].first, leftCol[r].second.c_str(), "p");
        else
            leg->AddEntry((TObject*)nullptr, "", "");
        if (r < rightCol.size())
            leg->AddEntry(rightCol[r].first, rightCol[r].second.c_str(), "lp");
        else
            leg->AddEntry((TObject*)nullptr, "", "");
    }
    leg->Draw();

    // Bottom panel: correction
    pad2->cd();
    TH1F* frB = new TH1F(("frROB_" + outPath).c_str(), "", 100, xlo, xhi);
    frB->SetStats(0);
    frB->SetMinimum(yLo);
    frB->SetMaximum(yHi);
    frB->GetYaxis()->SetTitle(RespAxis(q, corr.suffix).c_str());
    frB->GetYaxis()->SetTitleSize(0.07);
    frB->GetYaxis()->SetTitleOffset(0.95);
    frB->GetYaxis()->SetLabelSize(0.055);
    frB->GetXaxis()->SetTitle("Large-#it{R} (#it{b}#bar{#it{b}})-jet #it{p}_{T} [GeV]");
    frB->GetXaxis()->SetTitleSize(0.065);
    frB->GetXaxis()->SetTitleOffset(1.25);
    frB->GetXaxis()->SetLabelSize(0.055);
    frB->Draw("AXIS");
    drawDivider(yLo, yHi);
    for (auto* g : botGraphs) g->Draw("P same");
    TLatex tagB;
    tagB.SetNDC();
    tagB.SetTextFont(42);
    tagB.SetTextSize(0.055);
    tagB.DrawLatex(0.70, 0.88, corr.suffix.c_str());

    c->SaveAs(outPath.c_str());
    std::cout << "[Summary] Saved: " << outPath << "\n";

    for (auto* g : topGraphs) delete g;
    for (auto* g : botGraphs) delete g;
    for (auto* g : legendDummies) delete g;
    delete c;
}

inline std::vector<std::string> ParseCorrections(const std::string& sel) {
    if (sel.empty() || sel == "ALL" || sel == "all") return AllCorrectionSuffixes();
    std::vector<std::string> out;
    std::stringstream ss(sel);
    std::string item;
    while (std::getline(ss, item, ',')) {
        item.erase(0, item.find_first_not_of(" \t"));
        item.erase(item.find_last_not_of(" \t") + 1);
        if (!item.empty()) out.push_back(item);
    }
    return out;
}

} // namespace Summary
} // namespace RatioResponseCombined


// Called by Response_Summary_Combined.cxx in each channel's macros/ folder 
// (usage in the file header). Loops over the WPs, reads the CSV
// of each correction and draws every figure.
inline void Ratio_Response_Combined_Summary_Run(const std::string& corrections = "ALL",
                                                const std::string& topology = "combined") {
    using namespace RatioResponseCombined;
    using namespace RatioResponseCombined::Summary;
    AtlasUtils::LoadAtlasStyle();

    std::string topo = topology.empty() ? "combined" : topology;
    if (topo != "DB" && topo != "MJB" && topo != "combined") {
        std::cerr << "[Summary] Unknown topology '" << topo << "' — use DB, MJB, or combined.\n";
        return;
    }
    const bool wantDB = (topo == "DB" || topo == "combined");
    const bool wantMJB = (topo == "MJB" || topo == "combined");

    std::vector<std::string> wantCorr = ParseCorrections(corrections);
    if (wantCorr.empty()) {
        std::cerr << "[Summary] no corrections selected.\n";
        return;
    }

    for (const auto& wp : WPConfig::DefaultWPs()) {
        std::string wp_file = PathBuilder::WPtoFileWP(wp);
        std::string wp_label = PathBuilder::WPtoLabel(wp);

        // Load each correction's rows for the requested channel(s).
        std::vector<CorrData> corrs;
        for (const auto& suffix : wantCorr) {
            CorrData cd;
            cd.suffix = suffix;
            if (wantDB) {
                std::string key = "DB_" + suffix;
                cd.db = RowsForVariant(ReadCSV(CSVPath("DB", key, wp_file)), key, wp_file);
            }
            if (wantMJB) {
                std::string key = "MJB_" + suffix;
                cd.mjb = RowsForVariant(ReadCSV(CSVPath("MJB", key, wp_file)), key, wp_file);
            }
            if (cd.db.empty() && cd.mjb.empty()) {
                std::cerr << "[Summary] WP " << wp_file << ": no rows for correction '" << suffix
                          << "' — skipping it.\n";
                continue;
            }
            corrs.push_back(cd);
        }
        if (corrs.empty()) {
            std::cerr << "[Summary] WP " << wp_file
                      << ": no usable corrections (did you re-run Ratio_Response "
                      << "to write the 24-column CSVs?).\n";
            continue;
        }

        std::string outDir =
            (RunConfig::PlotDir() + "/Response/Combined_Ratios/") + wp_file + "/Summary_DataMC/";
        gSystem->Exec(("mkdir -p " + outDir).c_str());
        std::string respDir =
            (RunConfig::PlotDir() + "/Response/Combined_Ratios/") + wp_file + "/Summary_Response/";
        gSystem->Exec(("mkdir -p " + respDir).c_str());

        for (Qty q : {Qty::JES, Qty::JER}) {
            // Summary figures, one full set per double-ratio error treatment in
            // Summary_DataMC/<treatment>/. Only the double-ratio errors differ
            // between the three sets.
            for (DRErr fl : {DRErr::Uncorr, DRErr::Corr, DRErr::Full}) {
                DRErrFlavor() = fl;
                std::string fdir = outDir + DRErrTag(fl) + "/";
                gSystem->Exec(("mkdir -p " + fdir).c_str());
                for (int m = 0; m < kNMass; ++m) {
                    // One figure per correction
                    for (const auto& c : corrs) {
                        std::string out = fdir + QtyTag(q) + "_Summary_" + topo + "_" + c.suffix + "_"
                                          + MassTag(m) + "_" + wp_file + ".pdf";
                        DrawSummaryFigure(q, m, wp_label, topo, {c}, out, /*perCorrection=*/true);
                    }
                    // All corrections together (if more than one is selected)
                    if (corrs.size() > 1) {
                        std::string out = fdir + QtyTag(q) + "_Summary_" + topo + "_ALL_" + MassTag(m) + "_"
                                          + wp_file + ".pdf";
                        DrawSummaryFigure(q, m, wp_label, topo, corrs, out, /*perCorrection=*/false);
                    }
                }
                // Mass cuts overlaid, one figure per correction
                for (const auto& c : corrs) {
                    std::string out = fdir + QtyTag(q) + "_Summary_" + topo + "_" + c.suffix + "_AllMass_"
                                      + wp_file + ".pdf";
                    DrawSummaryFigureMassOverlay(q, wp_label, topo, c, out);
                }
            }

            // Response figures (b^Data and b^MC), one per correction and mass
            // cut, in Summary_Response/. They do not depend on the error treatment.
            for (int m = 0; m < kNMass; ++m) {
                for (const auto& c : corrs) {
                    std::string out = respDir + std::string(q == Qty::JES ? "bJES" : "bJER") + "_Response_"
                                      + topo + "_" + c.suffix + "_" + MassTag(m) + "_" + wp_file + ".pdf";
                    DrawResponseFigure(q, m, wp_label, topo, c, out);
                }
            }
            // Response figures with the mass cuts overlaid, one per correction
            for (const auto& c : corrs) {
                std::string out = respDir + std::string(q == Qty::JES ? "bJES" : "bJER") + "_Response_" + topo
                                  + "_" + c.suffix + "_AllMass_" + wp_file + ".pdf";
                DrawResponseFigureMassOverlay(q, wp_label, topo, c, out);
            }
        }
    }
    std::cout << "[Summary] Done.\n";
}

#endif // RATIO_RESPONSE_COMBINED_SUMMARY_H
