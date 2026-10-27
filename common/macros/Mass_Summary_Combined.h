/**
 * ============================================================================
 *  Mass_Summary_Combined.h - DB + MJB summary of the mass peak and width
 * ============================================================================
 * ===========================================
 * @author:    Thiziri AMEZZA
 * @email:     thiziri.amezza@cern.ch
 * @date:      10 June 2026
 * ===========================================
 *
 * Mass analogue of Response_Summary_Combined.h: DB and MJB results on one pT axis,
 * read from the CSVs written by Mass_Ratios_PeakWidth.h:
 *   b-JMDP = mu_data / mu_MC        (bb)-jet mass distribution peak
 *   b-JMDR = sigma_data / sigma_MC  (bb)-jet mass distribution resolution
 * Each channel adds the pT bins of its own Mass_Ratios run (inclusive: one point,
 * ptbins: one per bin), placed at pt_center, so e.g. DB inclusive with MJB in pT
 * bins works.
 *
 * How to run (from ZbbGamma/macros or ZbbJets/macros; the plots go to that
 * channel's Run{2,3}Plots), after Mass_Ratios_PeakWidth.cxx in both channels:
 *   export EASYBJETS_RUN=run3          same as for Mass_Fit and Mass_Ratios
 *   root 'Mass_Summary_Combined.cxx()'      ALL, combined, both inclusive
 *   root 'Mass_Summary_Combined.cxx("ALL","combined","inclusive","ptbins")'
 *   root 'Mass_Summary_Combined.cxx("bJER10v01,MuonCorr","DB","inclusive")'
 *   arguments: corrections, topology (DB | MJB | combined), DB pT mode, MJB pT mode
 *
 * Input:  ../../<channel>/Run3Plots/Mass/Mass_Ratios_Width/Combined/
 *             WP_<wp>_<ptmode>_<binning>/Combined_Mass_Ratio_Summary.csv
 * Output: <PlotDir>/Mass/Combined_Mass/<wp>/
 *   peak_Response_*.pdf, resolution_Response_*.pdf
 *       data and MC mu (sigma) vs pT, nominal on top, correction below
 *   {uncorrelated,correlation_aware,full_correlation,bootstrap}/b-JMD{P,R}_Summary_*.pdf
 *       single ratio (top) and double ratio (bottom), one folder per
 *       double-ratio error treatment
 *
 * Options (environment variables):
 *   EASYBJETS_RUN, EASYBJETS_BINNING, EASYBJETS_FIT_MODEL   must match Mass_Fit
 *   EASYBJETS_YRANGE_JMDP_SR, _JMDP_DR, _JMDP_COMP (and JMDR_*) = "lo,hi"
 *       fix the y range of the single ratio, double ratio and data-vs-MC plots
 *
 * ==========================================================================
 * @version:  1.0
 * ==========================================================================
 */

#ifndef MASS_SUMMARY_COMBINED_H
#define MASS_SUMMARY_COMBINED_H

#include "Response_Summary_Combined.h" // CorrectionColor, ParseCorrections, ChannelDirName,
                                       // PlotsDirName, ExpandRange, CorrDisp, ApplySummaryYRange, DRErr
#include "BinningConfig.h"
#include <fstream>
#include <sstream>

namespace MassCombined {

enum class MQ { JMDP, JMDR }; // mass distribution parameter: {peak , resolution}

inline std::string MQName(MQ q) { return q == MQ::JMDP ? "b-JMDP" : "b-JMDR"; }
inline std::string MQAxisSingle(MQ q) {
    return q == MQ::JMDP ? "b-JMDP (#mu^{Data}/#mu^{MC})" : "b-JMDR (#sigma^{Data}/#sigma^{MC})";
}
// #frac{b-JMDP^{bJR10v01}}{b-JMDP^{Nominal}}, kept short so it renders large enough in the legend
inline std::string MQAxisDR(MQ q, const std::string& corr) {
    const std::string base = (q == MQ::JMDP) ? "b-JMDP" : "b-JMDR";
    const std::string c = RatioResponseCombined::Summary::CorrDisp(corr); // bJER->bJR (display)
    return "#frac{" + base + "^{" + c + "}}{" + base + "^{Nominal}}";
}
// Full double-ratio definition for the legend in the DR pad:
//   #frac{b-JMDP^{Corrected} (#mu^{Data}/#mu^{MC})}{b-JMDP^{Nominal} (#mu^{Data}/#mu^{MC})}
inline std::string MQDRFull(MQ q, const std::string& corr) {
    const std::string base = (q == MQ::JMDP) ? "b-JMDP" : "b-JMDR";
    const std::string comp = (q == MQ::JMDP) ? "(#mu^{Data}/#mu^{MC})" : "(#sigma^{Data}/#sigma^{MC})";
    const std::string c = RatioResponseCombined::Summary::CorrDisp(corr); // bJER->bJR (display)
    return "#frac{" + base + "^{" + c + "} " + comp + "}{" + base + "^{Nominal} " + comp + "}";
}
// Component axis (the raw fitted peak / width), e.g. #mu^{Nominal} [GeV].
inline std::string MQCompAxis(MQ q, const std::string& tag) {
    return std::string(q == MQ::JMDP ? "#mu" : "#sigma") + "^{" + tag + "} [GeV]";
}

struct MRow {
    std::string channel, wp, variant;
    int pt_bin_index = -1;
    double pt_lo = 0, pt_hi = 0, pt_center = 0;
    double jmdp_nom = 0, jmdp_nom_err = 0, jmdp_var = 0, jmdp_var_err = 0;
    double jmdr_nom = 0, jmdr_nom_err = 0, jmdr_var = 0, jmdr_var_err = 0;
    double jmdp_dr = 0, jmdp_dr_err = 0, jmdr_dr = 0, jmdr_dr_err = 0;
    double jmdp_dr_err_corr = 0, jmdp_dr_err_full = 0;
    double jmdr_dr_err_corr = 0, jmdr_dr_err_full = 0;
    double mu_data_nom = 0, mu_data_nom_err = 0, mu_mc_nom = 0, mu_mc_nom_err = 0;
    double mu_data_var = 0, mu_data_var_err = 0, mu_mc_var = 0, mu_mc_var_err = 0;
    double sig_data_nom = 0, sig_data_nom_err = 0, sig_mc_nom = 0, sig_mc_nom_err = 0;
    double sig_data_var = 0, sig_data_var_err = 0, sig_mc_var = 0, sig_mc_var_err = 0;
    bool valid = false;
};

inline std::vector<MRow> ReadMassCSV(const std::string& path) {
    std::vector<MRow> out;
    std::ifstream in(path);
    if (!in.is_open()) {
        std::cerr << "[MassCombined] WARN: cannot open " << path << "\n";
        return out;
    }
    std::string line;
    bool first = true;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        if (first) {
            first = false;
            continue;
        }
        std::vector<std::string> t;
        std::stringstream ss(line);
        std::string item;
        while (std::getline(ss, item, ',')) t.push_back(item);
        if (t.size() < 35) {
            std::cerr << "[MassCombined] WARN: short line (" << t.size() << "): " << line << "\n";
            continue;
        }
        MRow r;
        try {
            r.channel = t[0];
            r.wp = t[1];
            r.variant = t[2];
            r.pt_bin_index = std::stoi(t[3]);
            r.pt_lo = std::stod(t[4]);
            r.pt_hi = std::stod(t[5]);
            r.pt_center = std::stod(t[6]);

            r.jmdp_nom = std::stod(t[7]);
            r.jmdp_nom_err = std::stod(t[8]);
            r.jmdp_var = std::stod(t[9]);
            r.jmdp_var_err = std::stod(t[10]);

            r.jmdr_nom = std::stod(t[11]);
            r.jmdr_nom_err = std::stod(t[12]);
            r.jmdr_var = std::stod(t[13]);
            r.jmdr_var_err = std::stod(t[14]);

            r.jmdp_dr = std::stod(t[15]);
            r.jmdp_dr_err = std::stod(t[16]);
            r.jmdr_dr = std::stod(t[17]);
            r.jmdr_dr_err = std::stod(t[18]);

            r.mu_data_nom = std::stod(t[19]);
            r.mu_data_nom_err = std::stod(t[20]);
            r.mu_mc_nom = std::stod(t[21]);
            r.mu_mc_nom_err = std::stod(t[22]);

            r.mu_data_var = std::stod(t[23]);
            r.mu_data_var_err = std::stod(t[24]);
            r.mu_mc_var = std::stod(t[25]);
            r.mu_mc_var_err = std::stod(t[26]);

            r.sig_data_nom = std::stod(t[27]);
            r.sig_data_nom_err = std::stod(t[28]);
            r.sig_mc_nom = std::stod(t[29]);
            r.sig_mc_nom_err = std::stod(t[30]);

            r.sig_data_var = std::stod(t[31]);
            r.sig_data_var_err = std::stod(t[32]);
            r.sig_mc_var = std::stod(t[33]);
            r.sig_mc_var_err = std::stod(t[34]);

            if (t.size() >= 39) {
                r.jmdp_dr_err_corr = std::stod(t[35]);
                r.jmdp_dr_err_full = std::stod(t[36]);
                r.jmdr_dr_err_corr = std::stod(t[37]);
                r.jmdr_dr_err_full = std::stod(t[38]);
            }
            r.valid = true;
        } catch (const std::exception&) {
            std::cerr << "[MassCombined] parse error: " << line << "\n";
            continue;
        }
        out.push_back(r);
    }
    return out;
}

// Mass CSV path. ptTag = "inclusive" | "ptbins"; binTag from BinningConfig.
inline std::string MassCSVPath(const std::string& channel, const std::string& wp_file,
                               const std::string& ptTag, const std::string& binTag) {
    return std::string("../../") + RatioResponseCombined::ChannelDirName(channel) + "/"
           + RatioResponseCombined::PlotsDirName() + "/Mass/Mass_Ratios_Width/Combined/WP_" + wp_file + "_"
           + ptTag + "_" + binTag + "/Combined_Mass_Ratio_Summary.csv";
}

inline std::vector<MRow> RowsForVariant(const std::vector<MRow>& all, const std::string& key) {
    std::vector<MRow> sel;
    for (const auto& r : all)
        if (r.variant == key) sel.push_back(r);
    return sel;
}

struct MCorr {
    std::string suffix;
    std::vector<MRow> db, mjb;
};

// How the double-ratio error is computed. The CSV holds three versions of
// it, and the summary draws a full set of figures for each (one folder per
// treatment). The one in use is set with DRErrFlavor() and read by sDRE.
//   Uncorr  errors added in quadrature, no correlation (largest error)
//           CSV columns 16 / 18 (b-JMDP / b-JMDR, counted from 0)
//   Corr    variant and nominal fully correlated, rho = 1      columns 35 / 37
//   Full    rho from correlations.txt                          columns 36 / 38
//   Boot    spread of the double ratio over the bootstrap replicas
//           (dr_errors.csv, BootstrapErrors.h); Full where not measured
// Corr and Full fall back to Uncorr when their value is missing or 0
// (CSVs written before these columns existed), so all three can always be drawn.
enum class DRErr { Uncorr, Corr, Full, Boot };
inline DRErr& DRErrFlavor() {
    static DRErr f = DRErr::Uncorr;
    return f;
}
inline const char* DRErrTag(DRErr f) {
    switch (f) {
    case DRErr::Corr: return "correlation_aware";
    case DRErr::Full: return "full_correlation";
    case DRErr::Boot: return "bootstrap";
    default: return "uncorrelated";
    }
}
inline const char* DRErrLabel(DRErr f) {
    switch (f) {
    case DRErr::Corr: return "Correlation-aware errors (#rho=1)";
    case DRErr::Full: return "Full-correlation errors (#rho from file)";
    case DRErr::Boot: return "Bootstrap errors";
    default: return "Uncorrelated errors";
    }
}

inline double sNom(const MRow& r, MQ q) { return q == MQ::JMDP ? r.jmdp_nom : r.jmdr_nom; }
inline double sNomE(const MRow& r, MQ q) { return q == MQ::JMDP ? r.jmdp_nom_err : r.jmdr_nom_err; }
inline double sVar(const MRow& r, MQ q) { return q == MQ::JMDP ? r.jmdp_var : r.jmdr_var; }
inline double sVarE(const MRow& r, MQ q) { return q == MQ::JMDP ? r.jmdp_var_err : r.jmdr_var_err; }
inline double sDR(const MRow& r, MQ q) { return q == MQ::JMDP ? r.jmdp_dr : r.jmdr_dr; }
inline double sDRE(const MRow& r, MQ q) {
    const double uncorr = (q == MQ::JMDP ? r.jmdp_dr_err : r.jmdr_dr_err);
    switch (DRErrFlavor()) {
    case DRErr::Corr: {
        double v = (q == MQ::JMDP ? r.jmdp_dr_err_corr : r.jmdr_dr_err_corr);
        return v > 0 ? v : uncorr;
    }
    case DRErr::Full: {
        double v = (q == MQ::JMDP ? r.jmdp_dr_err_full : r.jmdr_dr_err_full);
        return v > 0 ? v : uncorr;
    }
    case DRErr::Boot: {
        const std::string bin = BootstrapErrors::MassBinLabel(r.pt_lo, r.pt_hi);
        const double rel =
            BootstrapErrors::RelErr(r.channel, q == MQ::JMDP ? "bJMDP" : "bJMDR", r.variant, r.wp, bin);
        if (std::isfinite(rel)) return std::abs(sDR(r, q)) * rel;
        double v = (q == MQ::JMDP ? r.jmdp_dr_err_full : r.jmdr_dr_err_full);  // not measured: full
        return v > 0 ? v : uncorr;
    }
    default:
        return uncorr;
    }
}
inline double cDataNom(const MRow& r, MQ q) { return q == MQ::JMDP ? r.mu_data_nom : r.sig_data_nom; }
inline double cDataNomE(const MRow& r, MQ q) {
    return q == MQ::JMDP ? r.mu_data_nom_err : r.sig_data_nom_err;
}
inline double cMcNom(const MRow& r, MQ q) { return q == MQ::JMDP ? r.mu_mc_nom : r.sig_mc_nom; }
inline double cMcNomE(const MRow& r, MQ q) { return q == MQ::JMDP ? r.mu_mc_nom_err : r.sig_mc_nom_err; }
inline double cDataVar(const MRow& r, MQ q) { return q == MQ::JMDP ? r.mu_data_var : r.sig_data_var; }
inline double cDataVarE(const MRow& r, MQ q) {
    return q == MQ::JMDP ? r.mu_data_var_err : r.sig_data_var_err;
}
inline double cMcVar(const MRow& r, MQ q) { return q == MQ::JMDP ? r.mu_mc_var : r.sig_mc_var; }
inline double cMcVarE(const MRow& r, MQ q) { return q == MQ::JMDP ? r.mu_mc_var_err : r.sig_mc_var_err; }

inline TGraphErrors* MGraph(const std::vector<MRow>& rows, const std::function<double(const MRow&)>& gy,
                            const std::function<double(const MRow&)>& ge, int color, int marker,
                            double xOffset) {
    std::vector<double> x, y, ex, ey;
    for (const auto& r : rows) {
        if (!r.valid) continue;
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
    g->SetMarkerSize(1.2);
    g->SetLineWidth(2);
    return g;
}

// Half-width of the non-closure band around 1 in the double-ratio panel.
//   b-JMDP: +/-2 %
//   b-JMDR: no band. Setting kDrawJMDRBand = true draws a band of half-width
//           RMS(DR - 1) of the MuonCorr double ratio over the plotted pT bins.
// Returns -1 when no band is drawn.
inline double MassNonClosureHalfWidth(MQ q, const std::vector<MCorr>& corrs, bool wantDB, bool wantMJB) {
    if (q == MQ::JMDP) return 0.02;
    static const bool kDrawJMDRBand = false;
    if (!kDrawJMDRBand) return -1.0;
    const MCorr* mu = nullptr;
    for (const auto& c : corrs)
        if (c.suffix.find("MuonCorr") != std::string::npos) {
            mu = &c;
            break;
        }
    if (!mu) return -1.0;
    double s2 = 0.0;
    int n = 0;
    auto scan = [&](const std::vector<MRow>& rows) {
        for (const auto& r : rows) {
            if (!r.valid) continue;
            double dr = sDR(r, MQ::JMDR);
            if (dr <= 0) continue;
            s2 += (dr - 1.0) * (dr - 1.0);
            ++n;
        }
    };
    if (wantDB) scan(mu->db);
    if (wantMJB) scan(mu->mjb);
    return n ? std::sqrt(s2 / n) : -1.0;
}

// Single (b-JMDP or b-JMDR) ratio data/MC (top, nominal + corrections) and double ratio (bottom).
// ----------------------------------------------------------------------------
inline void DrawMassRatioFigure(MQ q, const std::string& wp_label, const std::string& topology,
                                const std::vector<MCorr>& corrs, const std::string& outPath,
                                bool perCorrection) {
    if (corrs.empty()) return;
    using namespace RatioResponseCombined;          // ChannelDirName, PlotsDirName
    using namespace RatioResponseCombined::Summary; // CorrectionColor, ExpandRange
    const bool wantDB = (topology == "DB" || topology == "combined");
    const bool wantMJB = (topology == "MJB" || topology == "combined");
    const bool combined = (topology == "combined");
    const double xlo = wantDB ? 200.0 : 450.0, xhi = wantMJB ? 1200.0 : 450.0;

    auto vN = [&](const MRow& r) {
        return sNom(r, q);
    };
    auto eN = [&](const MRow& r) {
        return sNomE(r, q);
    };
    auto vV = [&](const MRow& r) {
        return sVar(r, q);
    };
    auto eV = [&](const MRow& r) {
        return sVarE(r, q);
    };
    auto vD = [&](const MRow& r) {
        return sDR(r, q);
    };
    auto eD = [&](const MRow& r) {
        return sDRE(r, q);
    };
    const int nSeries = 1 + (int)corrs.size();
    auto regW = [](bool isDB) {
        return isDB ? 250.0 : 750.0;
    };
    auto so = [&](int s, bool isDB) {
        return regW(isDB) * 0.03 * (s - (nSeries - 1) / 2.0);
    };
    const int kNom = 20, kCorr = 22;

    std::vector<TGraphErrors*> top, bot;
    std::vector<std::pair<TGraphErrors*, std::string>> leg;
    const MCorr& first = corrs.front();
    bool nomAdded = false;
    auto addNom = [&](const std::vector<MRow>& rows, bool isDB) {
        auto* g = MGraph(rows, vN, eN, kBlack, kNom, so(0, isDB));
        if (g) {
            g->SetLineStyle(2);
            top.push_back(g);
            if (!nomAdded) {
                leg.push_back({g, "Nominal"});
                nomAdded = true;
            }
        }
    };
    if (wantDB) addNom(first.db, true);
    if (wantMJB) addNom(first.mjb, false);
    for (int ci = 0; ci < (int)corrs.size(); ++ci) {
        const auto& c = corrs[ci];
        int col = CorrectionColor(c.suffix);
        int s = ci + 1;
        TGraphErrors* lref = nullptr;
        if (wantDB) {
            auto* gv = MGraph(c.db, vV, eV, col, kCorr, so(s, true));
            auto* gd = MGraph(c.db, vD, eD, col, kCorr, so(s, true));
            if (gv) {
                top.push_back(gv);
                lref = gv;
            }
            if (gd) bot.push_back(gd);
        }
        if (wantMJB) {
            auto* gv = MGraph(c.mjb, vV, eV, col, kCorr, so(s, false));
            auto* gd = MGraph(c.mjb, vD, eD, col, kCorr, so(s, false));
            if (gv) {
                top.push_back(gv);
                if (!lref) lref = gv;
            }
            if (gd) bot.push_back(gd);
        }
        if (lref) leg.push_back({lref, CorrDisp(c.suffix)});
    }

    if (top.empty()) {
        std::cerr << "[MassCombined] no data for " << MQName(q) << " (" << topology << ") - skipping "
                  << outPath << "\n";
        for (auto* g : bot) delete g;
        return;
    }

    double tLo = 1e9, tHi = -1e9, bLo = 1e9, bHi = -1e9;
    for (auto* g : top) ExpandRange(g, tLo, tHi);
    for (auto* g : bot) ExpandRange(g, bLo, bHi);
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
    // Optional fixed y-range so Run 2 / Run 3 use the same scale (env override).
    {
        const std::string mq = (q == MQ::JMDP) ? "JMDP" : "JMDR";
        RatioResponseCombined::Summary::ApplySummaryYRange(mq + "_SR", tLo, tHi);
        RatioResponseCombined::Summary::ApplySummaryYRange(mq + "_DR", bLo, bHi);
    }

    TCanvas* c = new TCanvas(("c_" + outPath).c_str(), wp_label.c_str(), 800, 800);
    // Equal-size panels (50/50): the bottom double-ratio panel and the top single-ratio panel have same height.
    TPad* p1 = new TPad("pt", "", 0, 0.5, 1, 1);
    p1->SetBottomMargin(0.02);
    p1->SetLeftMargin(0.14);
    p1->SetTopMargin(0.08);
    p1->SetRightMargin(0.05);
    p1->Draw();
    TPad* p2 = new TPad("pb", "", 0, 0, 1, 0.5);
    p2->SetTopMargin(0.02);
    p2->SetBottomMargin(0.16);
    p2->SetLeftMargin(0.14);
    p2->SetRightMargin(0.05);
    p2->Draw();

    p1->cd();
    TH1F* fT = new TH1F(("fmt_" + outPath).c_str(), "", 100, xlo, xhi);
    fT->SetStats(0);
    fT->SetMinimum(tLo);
    fT->SetMaximum(tHi);
    fT->GetYaxis()->SetTitle(MQAxisSingle(q).c_str());
    fT->GetYaxis()->SetTitleSize(0.055);
    fT->GetYaxis()->SetTitleOffset(1.15);
    fT->GetYaxis()->SetLabelSize(0.045);
    fT->GetXaxis()->SetLabelSize(0);
    fT->Draw("AXIS");
    TLine* rT = new TLine(xlo, 1, xhi, 1);
    rT->SetLineStyle(2);
    rT->SetLineColor(kGray + 2);
    rT->Draw("same");
    if (combined) {
        TLine* s = new TLine(450, tLo, 450, tHi);
        s->SetLineStyle(2);
        s->SetLineColor(kGray + 2);
        s->Draw("same");
    }
    for (auto* g : top) g->Draw("P same");
    if (combined) {
        TLatex rg;
        rg.SetTextAlign(22);
        rg.SetTextSize(0.036);
        rg.SetTextColor(kGray + 2);
        rg.DrawLatex(0.5 * (xlo + 450.0), tLo + 0.05 * (tHi - tLo), "Z(#rightarrow b#bar{b}) + #gamma");
        rg.DrawLatex(0.5 * (450.0 + xhi), tLo + 0.05 * (tHi - tLo), "Z(#rightarrow b#bar{b}) + jets");
    }
    TLatex tx;
    tx.SetNDC();
    tx.SetTextFont(42);
    tx.SetTextSize(0.052);
    tx.DrawLatex(0.17, 0.86, AtlasLabel::Title().c_str());
    tx.SetTextSize(0.040);
    tx.DrawLatex(0.17, 0.79, RunConfig::LumiLabel().c_str());
    tx.DrawLatex(0.17, 0.73, "Mass fit");
    tx.DrawLatex(0.17, 0.67, "#it{R} = 1.0 UFO Soft Drop");
    tx.DrawLatex(0.17, 0.61, wp_label.c_str());
    TLegend* lg = new TLegend(0.77, 0.88 - 0.06 * leg.size(), 0.95, 0.88);
    lg->SetBorderSize(0);
    lg->SetFillStyle(0);
    lg->SetFillColor(kWhite);
    lg->SetTextSize(0.045);
    for (auto& e : leg) lg->AddEntry(e.first, e.second.c_str(), "lep");
    lg->Draw();

    p2->cd();
    TH1F* fB = new TH1F(("fmb_" + outPath).c_str(), "", 100, xlo, xhi);
    fB->SetStats(0);
    fB->SetMinimum(bLo);
    fB->SetMaximum(bHi);
    fB->GetYaxis()->SetTitle(MQAxisDR(q, perCorrection ? corrs.front().suffix : "Corrected").c_str());
    fB->GetYaxis()->SetNdivisions(510);
    fB->GetYaxis()->SetTitleSize(0.060);
    fB->GetYaxis()->SetTitleOffset(0.90);
    fB->GetYaxis()->SetLabelSize(0.050);
    fB->GetXaxis()->SetTitle("Large-#it{R} (#it{b}#bar{#it{b}})-jet #it{p}_{T} [GeV]");
    fB->GetXaxis()->SetTitleSize(0.060);
    fB->GetXaxis()->SetTitleOffset(1.05);
    fB->GetXaxis()->SetLabelSize(0.050);
    fB->Draw("AXIS");
    //
    double ncHalf = MassNonClosureHalfWidth(q, corrs, wantDB, wantMJB);
    if (ncHalf > 0) {
        TBox* ncB = new TBox(xlo, 1.0 - ncHalf, xhi, 1.0 + ncHalf);
        ncB->SetFillColorAlpha(kGray + 1, 0.35);
        ncB->SetLineWidth(0);
        ncB->Draw("same");
        TLegend* legBand = new TLegend(0.5, 0.80, 0.95, 0.97);
        legBand->SetBorderSize(0);
        legBand->SetFillStyle(0);
        legBand->SetTextSize(0.72 * fB->GetXaxis()->GetTitleSize());
        legBand->AddEntry(ncB, Form("Assigned non-closure unc. (#pm%g%%)", ncHalf * 100.0), "f");
        legBand->Draw();
    }
    TLine* rB = new TLine(xlo, 1, xhi, 1);
    rB->SetLineStyle(2);
    rB->SetLineColor(kGray + 2);
    rB->Draw("same");
    if (combined) {
        TLine* s = new TLine(450, bLo, 450, bHi);
        s->SetLineStyle(2);
        s->SetLineColor(kGray + 2);
        s->Draw("same");
    }
    for (auto* g : bot) g->Draw("P same");

    c->SaveAs(outPath.c_str());
    std::cout << "[MassCombined] Saved: " << outPath << "\n";
    for (auto* g : top) delete g;
    for (auto* g : bot) delete g;
    delete c;
}

// ----------------------------------------------------------------------------
// Peak / resolution Data-vs-MC response figure (nominal top, correction bottom).
// ----------------------------------------------------------------------------
inline void DrawMassResponseFigure(MQ q, const std::string& wp_label, const std::string& topology,
                                   const MCorr& corr, const std::string& outPath) {
    using namespace RatioResponseCombined;
    using namespace RatioResponseCombined::Summary;
    const bool wantDB = (topology == "DB" || topology == "combined");
    const bool wantMJB = (topology == "MJB" || topology == "combined");
    const bool combined = (topology == "combined");
    const double xlo = wantDB ? 200.0 : 450.0, xhi = wantMJB ? 1200.0 : 450.0;

    auto dN = [&](const MRow& r) {
        return cDataNom(r, q);
    };
    auto dNe = [&](const MRow& r) {
        return cDataNomE(r, q);
    };
    auto mN = [&](const MRow& r) {
        return cMcNom(r, q);
    };
    auto mNe = [&](const MRow& r) {
        return cMcNomE(r, q);
    };
    auto dV = [&](const MRow& r) {
        return cDataVar(r, q);
    };
    auto dVe = [&](const MRow& r) {
        return cDataVarE(r, q);
    };
    auto mV = [&](const MRow& r) {
        return cMcVar(r, q);
    };
    auto mVe = [&](const MRow& r) {
        return cMcVarE(r, q);
    };
    auto regW = [](bool isDB) {
        return isDB ? 250.0 : 750.0;
    };
    auto off = [&](bool isDB, bool isData) {
        return regW(isDB) * (isData ? -0.012 : +0.012);
    };
    const int kData = kBlack, kMC = kRed + 1, mData = 20, mMC = 22;

    std::vector<TGraphErrors*> top, bot;
    TGraphErrors *lD = nullptr, *lM = nullptr;
    auto fill = [&](const std::vector<MRow>& rows, bool isDB) {
        auto* a = MGraph(rows, dN, dNe, kData, mData, off(isDB, true));
        auto* b = MGraph(rows, mN, mNe, kMC, mMC, off(isDB, false));
        auto* aV = MGraph(rows, dV, dVe, kData, mData, off(isDB, true));
        auto* bV = MGraph(rows, mV, mVe, kMC, mMC, off(isDB, false));
        if (a) a->SetLineStyle(2);
        if (aV) aV->SetLineStyle(2);
        if (a) {
            top.push_back(a);
            if (!lD) lD = a;
        }
        if (b) {
            top.push_back(b);
            if (!lM) lM = b;
        }
        if (aV) bot.push_back(aV);
        if (bV) bot.push_back(bV);
    };
    if (wantDB) fill(corr.db, true);
    if (wantMJB) fill(corr.mjb, false);
    if (top.empty() && bot.empty()) {
        std::cerr << "[MassCombined] no " << MQName(q) << " component data for " << corr.suffix
                  << " - skipping " << outPath << "\n";
        return;
    }

    double yLo = 1e9, yHi = -1e9;
    for (auto* g : top) ExpandRange(g, yLo, yHi);
    for (auto* g : bot) ExpandRange(g, yLo, yHi);
    if (yLo > yHi) {
        yLo = 0;
        yHi = 1;
    }
    double d = 0.15 * (yHi - yLo) + 1e-6;
    yLo -= d;
    yHi += d;
    double yHiTop = yHi + 0.42 * (yHi - yLo);
    {
        const std::string mq = (q == MQ::JMDP) ? "JMDP" : "JMDR"; // optional fixed y-range (env)
        RatioResponseCombined::Summary::ApplySummaryYRange(mq + "_COMP", yLo, yHi);
        yHiTop = yHi + 0.42 * (yHi - yLo);
    }

    TCanvas* c = new TCanvas(("c_" + outPath).c_str(), wp_label.c_str(), 800, 800);
    TPad* p1 = new TPad("pn", "", 0, 0.50, 1, 1);
    p1->SetBottomMargin(0.02);
    p1->SetTopMargin(0.08);
    p1->SetLeftMargin(0.14);
    p1->SetRightMargin(0.05);
    p1->Draw();
    TPad* p2 = new TPad("pc", "", 0, 0, 1, 0.50);
    p2->SetTopMargin(0.02);
    p2->SetBottomMargin(0.20);
    p2->SetLeftMargin(0.14);
    p2->SetRightMargin(0.05);
    p2->Draw();
    auto div = [&](double lo, double hi) {
        if (!combined) return;
        TLine* s = new TLine(450, lo, 450, hi);
        s->SetLineStyle(2);
        s->SetLineColor(kGray + 2);
        s->Draw("same");
    };

    p1->cd();
    TH1F* fT = new TH1F(("fmrt_" + outPath).c_str(), "", 100, xlo, xhi);
    fT->SetStats(0);
    fT->SetMinimum(yLo);
    fT->SetMaximum(yHiTop);
    fT->GetYaxis()->SetTitle(MQCompAxis(q, "Nominal").c_str());
    fT->GetYaxis()->SetTitleSize(0.07);
    fT->GetYaxis()->SetTitleOffset(0.95);
    fT->GetYaxis()->SetLabelSize(0.055);
    fT->GetXaxis()->SetLabelSize(0);
    fT->Draw("AXIS");
    div(yLo, yHiTop);
    for (auto* g : top) g->Draw("P same");
    TLatex tx;
    tx.SetNDC();
    tx.SetTextFont(42);
    tx.SetTextSize(0.07);
    tx.DrawLatex(0.17, 0.85, AtlasLabel::Title(0.79).c_str());
    tx.SetTextSize(0.055);
    tx.DrawLatex(0.17, 0.77, "Nominal");
    tx.SetTextSize(0.045);
    tx.DrawLatex(0.17, 0.70, ("WP: " + wp_label).c_str());
    if (combined) {
        TLatex rg;
        rg.SetTextAlign(22);
        rg.SetTextSize(0.05);
        rg.SetTextColor(kGray + 2);
        rg.DrawLatex(0.5 * (xlo + 450.0), yLo + 0.06 * (yHiTop - yLo), "Z(#rightarrow b#bar{b}) + #gamma");
        rg.DrawLatex(0.5 * (450.0 + xhi), yLo + 0.06 * (yHiTop - yLo), "Z(#rightarrow b#bar{b}) + jets");
    }
    TLatex tag;
    tag.SetNDC();
    tag.SetTextFont(42);
    tag.SetTextSize(0.055);
    tag.DrawLatex(0.70, 0.85, "Nominal");
    TLegend* lg = new TLegend(0.70, 0.62, 0.92, 0.80);
    lg->SetBorderSize(1);
    lg->SetFillStyle(1001);
    lg->SetFillColor(kWhite);
    lg->SetTextSize(0.05);
    if (lD) lg->AddEntry(lD, "Data", "lep");
    if (lM) lg->AddEntry(lM, "MC", "lep");
    lg->Draw();

    p2->cd();
    TH1F* fB = new TH1F(("fmrb_" + outPath).c_str(), "", 100, xlo, xhi);
    fB->SetStats(0);
    fB->SetMinimum(yLo);
    fB->SetMaximum(yHi);
    fB->GetYaxis()->SetTitle(MQCompAxis(q, corr.suffix).c_str());
    fB->GetYaxis()->SetTitleSize(0.07);
    fB->GetYaxis()->SetTitleOffset(0.95);
    fB->GetYaxis()->SetLabelSize(0.055);
    fB->GetXaxis()->SetTitle("Large-#it{R} (#it{b}#bar{#it{b}})-jet #it{p}_{T} [GeV]");
    fB->GetXaxis()->SetTitleSize(0.065);
    fB->GetXaxis()->SetTitleOffset(1.25);
    fB->GetXaxis()->SetLabelSize(0.055);
    fB->Draw("AXIS");
    div(yLo, yHi);
    for (auto* g : bot) g->Draw("P same");
    TLatex tb;
    tb.SetNDC();
    tb.SetTextFont(42);
    tb.SetTextSize(0.055);
    tb.DrawLatex(0.70, 0.88, CorrDisp(corr.suffix).c_str());

    c->SaveAs(outPath.c_str());
    std::cout << "[MassCombined] Saved: " << outPath << "\n";
    for (auto* g : top) delete g;
    for (auto* g : bot) delete g;
    delete c;
}

} // namespace MassCombined


// ============================================================================
// called by the Mass_Summary_Combined.cxx
//   corrections, topology(DB|MJB|combined), dbPtMode/mjbPtMode (inclusive|ptbins)
// ============================================================================
inline void Mass_Ratios_Combined_Summary_Run(const std::string& corrections = "ALL",
                                             const std::string& topology = "combined",
                                             const std::string& dbPtMode = "inclusive",
                                             const std::string& mjbPtMode = "inclusive") {
    using namespace RatioResponseCombined;
    using namespace RatioResponseCombined::Summary;
    using namespace MassCombined;
    AtlasUtils::LoadAtlasStyle();

    std::string topo = topology.empty() ? "combined" : topology;
    if (topo != "DB" && topo != "MJB" && topo != "combined") {
        std::cerr << "[MassCombined] Unknown topology '" << topo << "'.\n";
        return;
    }
    const bool wantDB = (topo == "DB" || topo == "combined"), wantMJB = (topo == "MJB" || topo == "combined");
    auto ptTag = [](const std::string& m) {
        std::string s = m;
        for (auto& ch : s) ch = (char)std::tolower((unsigned char)ch);
        return (s == "ptbins") ? std::string("ptbins") : std::string("inclusive");
    };
    const std::string dbTag = ptTag(dbPtMode), mjbTag = ptTag(mjbPtMode), binTag = BinningConfig::ModeTag();

    std::vector<std::string> wantCorr = ParseCorrections(corrections);
    if (wantCorr.empty()) {
        std::cerr << "[MassCombined] no corrections.\n";
        return;
    }

    for (const auto& wp : WPConfig::DefaultWPs()) {
        std::string wp_file = PathBuilder::WPtoFileWP(wp), wp_label = PathBuilder::WPtoLabel(wp);
        std::vector<MRow> dbAll, mjbAll;
        if (wantDB) dbAll = ReadMassCSV(MassCSVPath("DB", wp_file, dbTag, binTag));
        if (wantMJB) mjbAll = ReadMassCSV(MassCSVPath("MJB", wp_file, mjbTag, binTag));

        std::vector<MCorr> corrs;
        for (const auto& s : wantCorr) {
            MCorr c;
            c.suffix = s;
            if (wantDB) c.db = RowsForVariant(dbAll, "DB_" + s);
            if (wantMJB) c.mjb = RowsForVariant(mjbAll, "MJB_" + s);
            if (c.db.empty() && c.mjb.empty()) {
                std::cerr << "[MassCombined] WP " << wp_file << ": no rows for '" << s
                          << "' - skipping it.\n";
                continue;
            }
            corrs.push_back(c);
        }
        if (corrs.empty()) {
            std::cerr << "[MassCombined] WP " << wp_file
                      << ": no usable corrections (re-run Mass_Ratios_PeakWidth?).\n";
            continue;
        }

        std::string outDir = (RunConfig::PlotDir() + "/Mass/Combined_Mass/") + wp_file + "/";
        gSystem->Exec(("mkdir -p " + outDir).c_str());

        for (MQ q : {MQ::JMDP, MQ::JMDR}) {
            // Data-vs-MC response figures - the same for all error treatments, drawn once.
            for (const auto& c : corrs) {
                DrawMassResponseFigure(q, wp_label, topo, c,
                                       outDir
                                           + (q == MQ::JMDP ? std::string("peak") : std::string("resolution"))
                                           + "_Response_" + topo + "_" + c.suffix + "_" + wp_file + ".pdf");
            }
            // Double-ratio figures
            for (MassCombined::DRErr fl :
                 {MassCombined::DRErr::Uncorr, MassCombined::DRErr::Corr, MassCombined::DRErr::Full,
                  MassCombined::DRErr::Boot}) {
                MassCombined::DRErrFlavor() = fl;
                std::string fdir = outDir + MassCombined::DRErrTag(fl) + "/";
                gSystem->Exec(("mkdir -p " + fdir).c_str());
                for (const auto& c : corrs) {
                    DrawMassRatioFigure(q, wp_label, topo, {c},
                                        fdir + MQName(q) + "_Summary_" + topo + "_" + c.suffix + "_" + wp_file
                                            + ".pdf",
                                        true);
                }
                if (corrs.size() > 1)
                    DrawMassRatioFigure(q, wp_label, topo, corrs,
                                        fdir + MQName(q) + "_Summary_" + topo + "_ALL_" + wp_file + ".pdf",
                                        false);
            }
        }
    }
    BootstrapErrors::PrintCounts("MassCombined");
    std::cout << "[MassCombined] Done.\n";
}

#endif // MASS_SUMMARY_COMBINED_H
