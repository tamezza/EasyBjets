/**
 * ============================================================================
 *  Response_Summary_Combined_Inclusive.h - DB + MJB summary of R_JES and R_JER
 *                                          without mass cuts
 * ============================================================================
 * ===========================================
 * @author:    Thiziri AMEZZA
 * @email:     thiziri.amezza@cern.ch
 * @date:      10 June 2026 
 * ===========================================
 *
 * The same figures as Response_Summary_Combined.h, without the mass-cut
 * dimension (colour => correction), for two cases:
 * 
 *   inclusive  one point per channel, read from the inclusive CSV of
 *              Response_Ratios_Inclusive.h, drawn at the centre of the
 *              channel's pT range
 * 
 *   ptbins     one point per pT bin, read from the pT-binned CSVs without mass
 *              cut of Response_Ratios_pTbins.cxx (folder R2pTbin_MassInclusive)
 * 
 * Figures, for each WP and correction:
 * 
 *   R_JES / R_JER summary   top: data/MC single ratio, nominal and correction
 *                           bottom: double ratio correction / nominal
 * 
 *   b_JES / b_JER response  top: data and MC (nominal), bottom: data and MC
 *                           (correction)
 * 
 * The summary figures are drawn once per double-ratio error treatment.
 *
 * How to run (from ZbbGamma/macros or ZbbJets/macros; the plots go to that
 * channel's Run{2,3}Plots):
 * 
 *   export EASYBJETS_RUN=run3          same as for the fits and ratios
 * 
 *   inclusive, after Response_Ratios_Inclusive.cxx in both channels:
 * 
 *     root 'Response_Summary_Combined_Inclusive.cxx()'
 *     root 'Response_Summary_Combined_Inclusive.cxx("bJER10v01,MuonCorr")'
 *     root 'Response_Summary_Combined_Inclusive.cxx("ALL","DB")'
 * 
 *   pT bins, after Response_Ratios_pTbins.cxx in both channels:
 * 
 *     root 'Response_Summary_Combined_pTbins.cxx()'
 * 
 *   arguments: corrections ("ALL" or a comma list without channel prefix),
 *              topology (combined | DB | MJB)
 *
 * Input:
 * 
 *   inclusive  ../../<channel>/Run3Plots/Response/Ratios/<ch>_Inclusive_Ratio/<wp>/
 *                  Summary_Inclusive_Ratios_<wp>.csv   (WriteSummaryInclusiveCSV)
 * 
 *   ptbins     ../../<channel>/Run3Plots/Response/Ratios/<ch>_Ratio_<variant>/<wp>/
 *                  R2pTbin_MassInclusive/Summary_Ratios_<wp>.csv
 * 
 * Output: <PlotDir>/Response/Combined_Ratios/<wp>/Summary_Inclusive/
 *         (Summary_pTbins_Inclusive/ in ptbins mode)
 * 
 *   bJES|bJER_ResponseInclusive_<topo>_<corr>_<wp>.pdf
 *   <treatment>/R_JES|R_JER_SummaryInclusive_<topo>_<corr|ALL>_<wp>.pdf
 *   <treatment> = uncorrelated | correlation_aware | full_correlation
 *
 * Options: EASYBJETS_RUN and EASYBJETS_FIT_MODEL must match the previous steps
 * (they set the Run{2,3}Plots folder).
 *
 * ==========================================================================
 * @version:  1.0
 * ==========================================================================
 */

#ifndef RATIO_RESPONSE_COMBINED_SUMMARY_INCLUSIVE_H
#define RATIO_RESPONSE_COMBINED_SUMMARY_INCLUSIVE_H

#include "Response_Summary_Combined.h" // Row, ReadCSV, RowsForVariant, Qty, the
                                       // row accessors, axis titles, colours, ...

namespace RatioResponseCombined {
namespace Summary {
namespace Inclusive {

// Inclusive CSV of a channel: one file per channel and WP, all variants inside.
inline std::string CSVPathInclusive(const std::string& channel, const std::string& wp_file) {
    return std::string("../../") + ChannelDirName(channel) + "/" + PlotsDirName() + "/Response/Ratios/"
           + channel + "_Inclusive_Ratio/" + wp_file + "/Summary_Inclusive_Ratios_" + wp_file + ".csv";
}

// Graph of value +/- error vs pt_center, with no mass-cut selection. needResp
// selects rows with the b_JES / b_JER columns, otherwise rows with the single ratios.
inline TGraphErrors* BuildInclGraph(const std::vector<Row>& rows, const std::function<double(const Row&)>& gy,
                                    const std::function<double(const Row&)>& ge, int color, int marker,
                                    double xOffset, bool needResp) {
    std::vector<double> x, y, ex, ey;
    for (const auto& r : rows) {
        if (needResp ? !r.hasResponse : !r.hasSingles) continue;
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
    g->SetMarkerSize(1.3);
    g->SetLineWidth(2);
    return g;
}

// Half-width of the non-closure band around 1 in the double-ratio panel.
//   R_JES: +/-1 %
//   R_JER: no band. Setting kDrawJERBand = true draws a band of half-width
//          RMS(DR - 1) of the MuonCorr double ratio over the plotted points.
// Returns -1 when no band is drawn.
inline double NonClosureHalfWidth(Qty q, const std::vector<CorrData>& corrs, bool wantDB, bool wantMJB) {
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
            if (!r.hasResponse) continue;
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

// Summary figure: data/MC single ratio on top (nominal: black circles,
// corrections: coloured triangles), double ratio correction / nominal below.
inline void DrawInclusiveRatioFigure(Qty q, const std::string& wp_label, const std::string& topology,
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

    const int nSeries = 1 + (int)corrs.size();
    auto regW = [](bool isDB) {
        return isDB ? 250.0 : 750.0;
    };
    auto seriesOff = [&](int sidx, bool isDB) {
        return regW(isDB) * 0.03 * (sidx - (nSeries - 1) / 2.0);
    };
    const int kMkNom = 20, kMkCorr = 22;

    std::vector<TGraphErrors*> topGraphs, botGraphs;
    std::vector<std::pair<TGraphErrors*, std::string>> legTop;
    const CorrData& first = corrs.front();
    bool nomAdded = false;
    auto addNom = [&](const std::vector<Row>& rows, bool isDB) {
        auto* g = BuildInclGraph(rows, valNom, errNom, kBlack, kMkNom, seriesOff(0, isDB), false);
        if (g) {
            g->SetLineStyle(2);
            topGraphs.push_back(g);
            if (!nomAdded) {
                legTop.push_back({g, "Nominal"});
                nomAdded = true;
            }
        }
    };
    if (wantDB) addNom(first.db, true);
    if (wantMJB) addNom(first.mjb, false);

    for (int ci = 0; ci < (int)corrs.size(); ++ci) {
        const auto& c = corrs[ci];
        const int col = CorrectionColor(c.suffix);
        const int s = ci + 1;
        TGraphErrors* legRef = nullptr;
        if (wantDB) {
            auto* gV = BuildInclGraph(c.db, valVar, errVar, col, kMkCorr, seriesOff(s, true), false);
            auto* gD = BuildInclGraph(c.db, valDR, errDR, col, kMkCorr, seriesOff(s, true), false);
            if (gV) {
                topGraphs.push_back(gV);
                legRef = gV;
            }
            if (gD) botGraphs.push_back(gD);
        }
        if (wantMJB) {
            auto* gV = BuildInclGraph(c.mjb, valVar, errVar, col, kMkCorr, seriesOff(s, false), false);
            auto* gD = BuildInclGraph(c.mjb, valDR, errDR, col, kMkCorr, seriesOff(s, false), false);
            if (gV) {
                topGraphs.push_back(gV);
                if (!legRef) legRef = gV;
            }
            if (gD) botGraphs.push_back(gD);
        }
        if (legRef) legTop.push_back({legRef, CorrDisp(c.suffix)});
    }

    if (topGraphs.empty()) {
        std::cerr << "[Summary] no inclusive data for " << QtyTag(q) << " (" << topology << ") — skipping "
                  << outPath << "\n";
        for (auto* g : botGraphs) delete g;
        return;
    }

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
    tHi += 0.35 * (tHi - tLo);

    TCanvas* c = new TCanvas(("c_" + outPath).c_str(), wp_label.c_str(), 800, 800);
    // Two panels of the same height.
    TPad* pad1 = new TPad("pad_top", "", 0.0, 0.5, 1.0, 1.0);
    pad1->SetBottomMargin(0.02);
    pad1->SetLeftMargin(0.14);
    pad1->SetTopMargin(0.08);
    pad1->SetRightMargin(0.05);
    pad1->Draw();
    TPad* pad2 = new TPad("pad_bot", "", 0.0, 0.0, 1.0, 0.5);
    pad2->SetTopMargin(0.02);
    pad2->SetBottomMargin(0.16);
    pad2->SetLeftMargin(0.14);
    pad2->SetRightMargin(0.05);
    pad2->Draw();

    pad1->cd();
    TH1F* frT = new TH1F(("frIT_" + outPath).c_str(), "", 100, xlo, xhi);
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
    tex.DrawLatex(0.17, 0.86, AtlasLabel::Title().c_str());
    tex.SetTextSize(0.040);
    tex.DrawLatex(0.17, 0.79, RunConfig::LumiLabel().c_str());
    tex.DrawLatex(0.17, 0.73, "Inclusive Mass: [50 - 150] GeV");
    tex.DrawLatex(0.17, 0.67, "#it{R} = 1.0 UFO Soft Drop");
    tex.DrawLatex(0.17, 0.61, wp_label.c_str()); // e.g. GN2X (0.46% QCD WP)
    TLegend* legA = new TLegend(0.77, 0.88 - 0.06 * legTop.size(), 0.95, 0.88);
    legA->SetBorderSize(0);
    legA->SetFillStyle(0);
    legA->SetFillColor(kWhite);
    legA->SetTextSize(0.045);
    for (auto& e : legTop) legA->AddEntry(e.first, e.second.c_str(), "lep");
    legA->Draw();

    pad2->cd();
    TH1F* frB = new TH1F(("frIB_" + outPath).c_str(), "", 100, xlo, xhi);
    frB->SetStats(0);
    frB->SetMinimum(bLo);
    frB->SetMaximum(bHi);
    frB->GetYaxis()->SetTitle(DoubleRatioAxis(q, perCorrection ? corrs.front().suffix : "Corrected").c_str());
    frB->GetYaxis()->SetNdivisions(506);
    frB->GetYaxis()->SetTitleSize(0.060);
    frB->GetYaxis()->SetTitleOffset(0.90);
    frB->GetYaxis()->SetLabelSize(0.050);
    frB->GetXaxis()->SetTitle("Large-#it{R} (#it{b}#bar{#it{b}})-jet #it{p}_{T} [GeV]");
    frB->GetXaxis()->SetTitleSize(0.060);
    frB->GetXaxis()->SetTitleOffset(1.05);
    frB->GetXaxis()->SetLabelSize(0.050);
    frB->Draw("AXIS");
    // Non-closure band around 1 (R_JES only, see NonClosureHalfWidth).
    double ncHalf = NonClosureHalfWidth(q, corrs, wantDB, wantMJB);
    if (ncHalf > 0) {
        TBox* ncB = new TBox(xlo, 1.0 - ncHalf, xhi, 1.0 + ncHalf);
        ncB->SetFillColorAlpha(kGray + 1, 0.35);
        ncB->SetLineWidth(0);
        ncB->Draw("same");
        // Labelled as an assigned uncertainty.
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
    delete c;
}

// Response figure: b^Data and b^MC for the nominal (top) and the correction (bottom). 
// Data = black circles (dashed error bars), MC = red triangles.
inline void DrawInclusiveResponseFigure(Qty q, const std::string& wp_label, const std::string& topology,
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
    const int mData = 20, mMC = 22; // circle = data (dashed), triangle = MC (solid)

    std::vector<TGraphErrors*> topGraphs, botGraphs;
    TGraphErrors *legD = nullptr, *legM = nullptr;
    auto fill = [&](const std::vector<Row>& rows, bool isDB) {
        auto* dN = BuildInclGraph(rows, dataNom, dataNomE, kData, mData, off(isDB, true), true);
        auto* mN = BuildInclGraph(rows, mcNom, mcNomE, kMC, mMC, off(isDB, false), true);
        auto* dV = BuildInclGraph(rows, dataVar, dataVarE, kData, mData, off(isDB, true), true);
        auto* mV = BuildInclGraph(rows, mcVar, mcVarE, kMC, mMC, off(isDB, false), true);
        if (dN) dN->SetLineStyle(2);
        if (dV) dV->SetLineStyle(2);
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
    if (wantDB) fill(corr.db, true);
    if (wantMJB) fill(corr.mjb, false);

    if (topGraphs.empty() && botGraphs.empty()) {
        std::cerr << "[Summary] no inclusive b_" << (q == Qty::JES ? "JES" : "JER") << " data for "
                  << corr.suffix << " (" << topology << ") — skipping " << outPath << "\n";
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
    double yHiTop = yHi + 0.42 * (yHi - yLo);

    TCanvas* c = new TCanvas(("c_" + outPath).c_str(), wp_label.c_str(), 800, 800);
    TPad* pad1 = new TPad("pad_nom", "", 0.0, 0.50, 1.0, 1.0);
    pad1->SetBottomMargin(0.02);
    pad1->SetTopMargin(0.08);
    pad1->SetLeftMargin(0.14);
    pad1->SetRightMargin(0.05);
    pad1->Draw();
    TPad* pad2 = new TPad("pad_cor", "", 0.0, 0.0, 1.0, 0.50);
    pad2->SetTopMargin(0.02);
    pad2->SetBottomMargin(0.20);
    pad2->SetLeftMargin(0.14);
    pad2->SetRightMargin(0.05);
    pad2->Draw();

    auto divider = [&](double lo, double hi) {
        if (!combined) return;
        TLine* s = new TLine(450.0, lo, 450.0, hi);
        s->SetLineStyle(2);
        s->SetLineColor(kGray + 2);
        s->Draw("same");
    };

    pad1->cd();
    TH1F* frT = new TH1F(("frIRT_" + outPath).c_str(), "", 100, xlo, xhi);
    frT->SetStats(0);
    frT->SetMinimum(yLo);
    frT->SetMaximum(yHiTop);
    frT->GetYaxis()->SetTitle(RespAxis(q, "Nominal").c_str());
    frT->GetYaxis()->SetTitleSize(0.07);
    frT->GetYaxis()->SetTitleOffset(0.95);
    frT->GetYaxis()->SetLabelSize(0.055);
    frT->GetXaxis()->SetLabelSize(0);
    frT->Draw("AXIS");
    divider(yLo, yHiTop);
    for (auto* g : topGraphs) g->Draw("P same");
    TLatex tex;
    tex.SetNDC();
    tex.SetTextFont(42);
    tex.SetTextSize(0.07);
    tex.DrawLatex(0.17, 0.85, AtlasLabel::Title(0.79).c_str());
    tex.SetTextSize(0.055);
    tex.DrawLatex(0.17, 0.77, "Inclusive");
    tex.SetTextSize(0.045);
    tex.DrawLatex(0.17, 0.70, ("WP: " + wp_label).c_str());
    if (combined) {
        TLatex reg;
        reg.SetTextAlign(22);
        reg.SetTextSize(0.05);
        reg.SetTextColor(kGray + 2);
        reg.DrawLatex(0.5 * (xlo + 450.0), yLo + 0.06 * (yHiTop - yLo), "Z(#rightarrow b#bar{b}) + #gamma");
        reg.DrawLatex(0.5 * (450.0 + xhi), yLo + 0.06 * (yHiTop - yLo), "Z(#rightarrow b#bar{b}) + jets");
    }
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

    pad2->cd();
    TH1F* frB = new TH1F(("frIRB_" + outPath).c_str(), "", 100, xlo, xhi);
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
    divider(yLo, yHi);
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

} // namespace Inclusive
} // namespace Summary
} // namespace RatioResponseCombined


// Entry point, called by Response_Summary_Combined_Inclusive.cxx (ptMode "inclusive") 
// and Response_Summary_Combined_pTbins.cxx (ptMode "ptbins").
// Usage and inputs in the file header.
inline void Ratio_Response_Combined_Summary_Inclusive_Run(const std::string& corrections = "ALL",
                                                          const std::string& topology = "combined",
                                                          const std::string& ptMode = "inclusive") {
    const bool binnedIncl = (ptMode == "ptbins");
    using namespace RatioResponseCombined;
    using namespace RatioResponseCombined::Summary;
    using namespace RatioResponseCombined::Summary::Inclusive;
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

        // Inclusive: one CSV per channel with every variant, read once here.
        // pT bins: one CSV per variant, read in the loop below.
        std::vector<Row> dbAll, mjbAll;
        if (!binnedIncl) {
            if (wantDB) dbAll = ReadCSV(CSVPathInclusive("DB", wp_file));
            if (wantMJB) mjbAll = ReadCSV(CSVPathInclusive("MJB", wp_file));
        }

        std::vector<CorrData> corrs;
        for (const auto& suffix : wantCorr) {
            CorrData cd;
            cd.suffix = suffix;
            if (binnedIncl) {
                if (wantDB)
                    cd.db = RowsForVariant(ReadCSV(CSVPathMassInclusive("DB", "DB_" + suffix, wp_file)),
                                           "DB_" + suffix, wp_file);
                if (wantMJB)
                    cd.mjb = RowsForVariant(ReadCSV(CSVPathMassInclusive("MJB", "MJB_" + suffix, wp_file)),
                                            "MJB_" + suffix, wp_file);
            } else {
                if (wantDB) cd.db = RowsForVariant(dbAll, "DB_" + suffix, wp_file);
                if (wantMJB) cd.mjb = RowsForVariant(mjbAll, "MJB_" + suffix, wp_file);
            }
            if (cd.db.empty() && cd.mjb.empty()) {
                std::cerr << "[Summary] WP " << wp_file << ": no inclusive rows for '" << suffix
                          << "' — skipping it.\n";
                continue;
            }
            corrs.push_back(cd);
        }
        if (corrs.empty()) {
            std::cerr << "[Summary] WP " << wp_file
                      << ": no usable corrections (did you re-run Ratio_Response_Inclusive?).\n";
            continue;
        }

        std::string subdir = binnedIncl ? "Summary_pTbins_Inclusive/" : "Summary_Inclusive/";
        std::string outDir = (RunConfig::PlotDir() + "/Response/Combined_Ratios/") + wp_file + "/" + subdir;
        gSystem->Exec(("mkdir -p " + outDir).c_str());

        for (Qty q : {Qty::JES, Qty::JER}) {
            // Response figures, the same for all error treatments: drawn once.
            for (const auto& c : corrs) {
                std::string outR = outDir + std::string(q == Qty::JES ? "bJES" : "bJER")
                                   + "_ResponseInclusive_" + topo + "_" + c.suffix + "_" + wp_file + ".pdf";
                DrawInclusiveResponseFigure(q, wp_label, topo, c, outR);
            }
            // Summary figures, one set per error treatment, in <treatment>/.
            for (DRErr fl : {DRErr::Uncorr, DRErr::Corr, DRErr::Full}) {
                DRErrFlavor() = fl;
                std::string fdir = outDir + DRErrTag(fl) + "/";
                gSystem->Exec(("mkdir -p " + fdir).c_str());
                for (const auto& c : corrs) {
                    std::string out = fdir + QtyTag(q) + "_SummaryInclusive_" + topo + "_" + c.suffix + "_"
                                      + wp_file + ".pdf";
                    DrawInclusiveRatioFigure(q, wp_label, topo, {c}, out, /*perCorrection=*/true);
                }
                if (corrs.size() > 1) {
                    std::string out =
                        fdir + QtyTag(q) + "_SummaryInclusive_" + topo + "_ALL_" + wp_file + ".pdf";
                    DrawInclusiveRatioFigure(q, wp_label, topo, corrs, out, /*perCorrection=*/false);
                }
            } // error treatments
        } // JES, JER
    } // WPs
    std::cout << "[Summary] Done (" << (binnedIncl ? "pT bins, inclusive mass" : "inclusive") << ").\n";
}

#endif // RATIO_RESPONSE_COMBINED_SUMMARY_INCLUSIVE_H
