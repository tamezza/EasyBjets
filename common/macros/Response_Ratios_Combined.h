/**
 * ============================================================================
 *  Response_Ratios_Combined.h - DB + MJB double ratios on one pT axis
 * ============================================================================
 * ===========================================
 * @author:    Thiziri AMEZZA
 * @email:     thiziri.amezza@cern.ch
 * @date:      10 June 2026 
 * ===========================================
 *
 * The response ratios are measured separately in each channel, in different
 * pT ranges:
 *   DB  (Z(bb)+gamma, ZbbGamma)   200 < pT < 450 GeV
 *   MJB (Z(bb)+jets,  ZbbJets)    450 < pT < 1200 GeV
 * This file reads the CSV of both channels and draws them on one canvas from
 * 200 to 1200 GeV, with a dashed line at 450 GeV between the two.
 *
 * For each correction (bJER10v00, bJER10v01, MuonCorr) found in both channels
 * and each WP, it draws 7 plots of the double ratio correction / nominal vs pT:
 * 
 *   R_JES_DR_vs_pT_<corr>_<wp>_Combined.pdf                   uncorrelated errors
 *   sigma_insitu_DR_vs_pT_<corr>_<wp>_Combined.pdf            uncorrelated errors
 *   sigma_insitu_DR_vs_pT_<corr>_<wp>_Combined_CorrAware.pdf  correlation-aware errors
 *   sigma_insitu_DR_vs_pT_<corr>_<wp>_Combined_Cmp.pdf        both errors overlaid
 *   R_JER_DR_vs_pT_<corr>_<wp>_Combined.pdf                   uncorrelated errors
 *   R_JER_DR_vs_pT_<corr>_<wp>_Combined_CorrAware.pdf         correlation-aware errors
 *   R_JER_DR_vs_pT_<corr>_<wp>_Combined_Cmp.pdf               both errors overlaid
 * 
 * into <PlotDir>/Response/Combined_Ratios/<wp>/.
 * Colour : mass cut (red 50-80, blue 80-110, green >= 110 GeV), 
 * marker : channel (circle DB, triangle MJB).
 *
 * Usage (from ZbbGamma/macros or ZbbJets/macros; the plots go to that
 * channel's Run{2,3}Plots), after Response_Ratios_pTbins_Mcuts.cxx("ALL") in
 * both channels with the same EASYBJETS_RUN:
 *   root 'Response_Ratios_Combined.cxx()'              every correction
 *   root 'Response_Ratios_Combined.cxx("bJER10v01")'   one correction
 *
 * Input: Summary_Ratios_<wp>.csv written by GenerateSummaryCSV in
 * Response_Ratios_pTbins_Mcuts.h, one per channel and correction (path in CSVPath below).
 *
 * This file also provides the CSV reader (Row, ReadCSV, RowsForVariant), the
 * CSV paths (ChannelDirName, PlotsDirName, CSVPath, CSVPathMassInclusive) and
 * the mass-cut helpers (MassBinIndex, MassBinDisplayLabel, MassCutColor) used
 * by Response_Summary_Combined.h, and through it by
 * Response_Summary_Combined_Inclusive.h and Mass_Summary_Combined.h.
 * Re-run those summaries after changing any of them.
 *
 * ==========================================================================
 * @version:  1.0
 * ==========================================================================
 */

#ifndef RATIO_RESPONSE_COMBINED_H
#define RATIO_RESPONSE_COMBINED_H

#include "VariantConfig.h"
#include "CommonUtils.h"
#include "RunEra.h"

#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>

#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TH1F.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TStyle.h>
#include <TSystem.h>


namespace RatioResponseCombined {

// One row of the CSV: one correction, mass cut and pT bin. "_dr" values are
// double ratios correction / nominal.
struct Row {
    std::string channel;  // "DB" / "MJB"
    std::string wp;       // "0p46" etc.
    std::string variant;  // DB_bJER10v00, MuonCorr, ...
    std::string mass_bin; // "50-80GeV", "80-110GeV", "#geq110 GeV"
    int pt_bin_index = -1;
    double pt_lo = 0, pt_hi = 0, pt_center = 0;
    double rjes_dr = 0, rjes_dr_err = 0; // uncorrelated error
    double sig_insitu_dr = 0, sig_insitu_dr_err_unc = 0, sig_insitu_dr_err_corr = 0;
    double rjer_dr = 0, rjer_dr_err_unc = 0, rjer_dr_err_corr = 0;
    // columns 41-44; the resolution "_full" values equal the corr-aware ones
    double rjes_dr_err_corr = 0, rjes_dr_err_full = 0;
    double sig_insitu_dr_err_full = 0, rjer_dr_err_full = 0;
    bool hasFullErr = false;
    // Single data/MC ratios, nominal and correction. 
    double rjes_nom = 0, rjes_nom_err = 0, rjes_var = 0, rjes_var_err = 0;
    double rjer_nom = 0, rjer_nom_err = 0, rjer_var = 0, rjer_var_err = 0;
    bool hasSingles = false;
    // b_JES / b_JER data and MC values, nominal and correction. 
    double bjes_data_nom = 0, bjes_data_nom_err = 0, bjes_mc_nom = 0, bjes_mc_nom_err = 0;
    double bjes_data_var = 0, bjes_data_var_err = 0, bjes_mc_var = 0, bjes_mc_var_err = 0;
    double bjer_data_nom = 0, bjer_data_nom_err = 0, bjer_mc_nom = 0, bjer_mc_nom_err = 0;
    double bjer_data_var = 0, bjer_data_var_err = 0, bjer_mc_var = 0, bjer_mc_var_err = 0;
    bool hasResponse = false;
};

// Read a Summary_Ratios CSV. The first line is the header; empty lines and
// lines starting with # are skipped. 
inline std::vector<Row> ReadCSV(const std::string& path) {
    std::vector<Row> out;
    std::ifstream in(path);
    if (!in.is_open()) {
        std::cerr << "[Combined] WARN: cannot open " << path << "\n";
        return out;
    }
    std::string line;
    bool first = true;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        if (first) {
            first = false;
            continue;
        } // header

        // Split on commas (the format has no quoted fields)
        std::vector<std::string> tok;
        std::stringstream ss(line);
        std::string item;
        while (std::getline(ss, item, ',')) tok.push_back(item);
        if (tok.size() < 16) {
            std::cerr << "[Combined] WARN: malformed CSV line (only " << tok.size() << " fields): " << line
                      << "\n";
            continue;
        }
        Row r;
        try {
            r.channel = tok[0];
            r.wp = tok[1];
            r.variant = tok[2];
            r.mass_bin = tok[3];
            r.pt_bin_index = std::stoi(tok[4]);
            r.pt_lo = std::stod(tok[5]);
            r.pt_hi = std::stod(tok[6]);
            r.pt_center = std::stod(tok[7]);
            r.rjes_dr = std::stod(tok[8]);
            r.rjes_dr_err = std::stod(tok[9]);
            r.sig_insitu_dr = std::stod(tok[10]);
            r.sig_insitu_dr_err_unc = std::stod(tok[11]);
            r.sig_insitu_dr_err_corr = std::stod(tok[12]);
            r.rjer_dr = std::stod(tok[13]);
            r.rjer_dr_err_unc = std::stod(tok[14]);
            r.rjer_dr_err_corr = std::stod(tok[15]);
            // Columns 16-23: single data/MC ratios
            if (tok.size() >= 24) {
                r.rjes_nom = std::stod(tok[16]);
                r.rjes_nom_err = std::stod(tok[17]);
                r.rjes_var = std::stod(tok[18]);
                r.rjes_var_err = std::stod(tok[19]);
                r.rjer_nom = std::stod(tok[20]);
                r.rjer_nom_err = std::stod(tok[21]);
                r.rjer_var = std::stod(tok[22]);
                r.rjer_var_err = std::stod(tok[23]);
                r.hasSingles = true;
            }
            // Columns 24-39: b_JES / b_JER data and MC values
            if (tok.size() >= 40) {
                r.bjes_data_nom = std::stod(tok[24]);
                r.bjes_data_nom_err = std::stod(tok[25]);
                r.bjes_mc_nom = std::stod(tok[26]);
                r.bjes_mc_nom_err = std::stod(tok[27]);
                r.bjes_data_var = std::stod(tok[28]);
                r.bjes_data_var_err = std::stod(tok[29]);
                r.bjes_mc_var = std::stod(tok[30]);
                r.bjes_mc_var_err = std::stod(tok[31]);
                r.bjer_data_nom = std::stod(tok[32]);
                r.bjer_data_nom_err = std::stod(tok[33]);
                r.bjer_mc_nom = std::stod(tok[34]);
                r.bjer_mc_nom_err = std::stod(tok[35]);
                r.bjer_data_var = std::stod(tok[36]);
                r.bjer_data_var_err = std::stod(tok[37]);
                r.bjer_mc_var = std::stod(tok[38]);
                r.bjer_mc_var_err = std::stod(tok[39]);
                r.hasResponse = true;
            }
            // Columns 40-43: correlated double-ratio errors
            if (tok.size() >= 44) {
                r.rjes_dr_err_corr = std::stod(tok[40]);
                r.rjes_dr_err_full = std::stod(tok[41]);
                r.sig_insitu_dr_err_full = std::stod(tok[42]);
                r.rjer_dr_err_full = std::stod(tok[43]);
                r.hasFullErr = true;
            }
        } catch (const std::exception& e) {
            std::cerr << "[Combined] WARN: parse error on line: " << line << " (" << e.what() << ")\n";
            continue;
        }
        out.push_back(r);
    }
    return out;
}

// Paths to the CSVs of both channels. RunConfig::PlotDir() ("../Run3Plots")
// points only to the channel the macro is run from, so the path to each
// channel is built explicitly, starting from <channel>/macros/:
//   ../../<ZbbGamma|ZbbJets>/<Run2Plots|Run3Plots>/Response/Ratios/
//       <channel>_Ratio_<key>/<wp>/R2pTbin/Summary_Ratios_<wp>.csv
// e.g. ../../ZbbGamma/Run3Plots/Response/Ratios/DB_Ratio_DB_bJER10v01/0p46/...
// This works from either channel's macros/ folder.
inline std::string ChannelDirName(const std::string& channel) {
    return (channel == "MJB") ? "ZbbJets" : "ZbbGamma";
}

inline std::string PlotsDirName() {
    // Last part of PlotDir(): "Run2Plots" or "Run3Plots" (with the fit-model
    // suffix if one is set), the same name in both channels.
    std::string p = RunConfig::PlotDir(); // "../Run2Plots" or "../Run3Plots"
    auto pos = p.find_last_of('/');
    return (pos == std::string::npos) ? p : p.substr(pos + 1);
}

inline std::string CSVPath(const std::string& channel, const std::string& variant_label,
                           const std::string& wp_file) {
    return std::string("../../") + ChannelDirName(channel) + "/" + PlotsDirName() + "/Response/Ratios/"
           + channel + "_Ratio_" + variant_label + "/" + wp_file + "/R2pTbin/Summary_Ratios_" + wp_file
           + ".csv";
}

// Same CSV for the run without mass cuts (Response_Ratios_pTbins): pT bins and
// a single "MassInclusive" mass bin, in R2pTbin_MassInclusive/.
inline std::string CSVPathMassInclusive(const std::string& channel, const std::string& variant_label,
                                        const std::string& wp_file) {
    return std::string("../../") + ChannelDirName(channel) + "/" + PlotsDirName() + "/Response/Ratios/"
           + channel + "_Ratio_" + variant_label + "/" + wp_file + "/R2pTbin_MassInclusive/Summary_Ratios_"
           + wp_file + ".csv";
}

// The corrections that have a CSV in both channels for a WP. Variant keys are
// "<channel>_<correction>" (e.g. DB_bJER10v00 and MJB_bJER10v00 pair as
// "bJER10v00"); the nominal keys "DB" and "MJB" are skipped, since their
// double ratio is 1 by definition.
struct VariantPair {
    std::string suffix;  // e.g. "bJER10v00"
    std::string db_key;  // e.g. "DB_bJER10v00"
    std::string mjb_key; // e.g. "MJB_bJER10v00"
};

inline std::vector<VariantPair> FindCommonVariants(const std::string& wp_file) {
    std::vector<std::string> db_keys = GetPrimaryVariantKeys("DB");
    std::vector<std::string> mjb_keys = GetPrimaryVariantKeys("MJB");

    auto exists = [&](const std::string& path) {
        std::ifstream in(path);
        return in.good();
    };
    auto suffix_of = [](const std::string& s) -> std::string {
        auto p = s.find('_');
        return (p == std::string::npos) ? std::string() : s.substr(p + 1);
    };

    std::vector<VariantPair> pairs;
    for (const auto& dk : db_keys) {
        std::string ds = suffix_of(dk);
        if (ds.empty()) continue; // skip the nominal ("DB" has no suffix)
        std::string db_csv = CSVPath("DB", dk, wp_file);
        if (!exists(db_csv)) continue;
        for (const auto& mk : mjb_keys) {
            std::string ms = suffix_of(mk);
            if (ms != ds) continue;
            std::string mjb_csv = CSVPath("MJB", mk, wp_file);
            if (!exists(mjb_csv)) continue;
            pairs.push_back({ds, dk, mk});
            break;
        }
    }
    return pairs;
}

// Rows of one variant key and WP.
inline std::vector<Row> RowsForVariant(const std::vector<Row>& all, const std::string& variantLabel,
                                       const std::string& wp_file) {
    std::vector<Row> sel;
    for (const auto& r : all) {
        if (r.variant == variantLabel && r.wp == wp_file) sel.push_back(r);
    }
    return sel;
}

// Colour of each mass cut, the same as in the per-channel plots.
inline int MassCutColor(int idx) {
    static const int cols[3] = {kRed, kBlue, kGreen + 2};
    return (idx >= 0 && idx < 3) ? cols[idx] : kBlack;
}

// Mass-cut index from the mass_bin column of the CSV, as written by
// GenerateSummaryCSV (massCutDirNames[] in Response_Ratios_pTbins_Mcuts.h):
//   0 = "50-80GeV", 1 = "80-110GeV", 2 = "MassGE110GeV". Returns -1 otherwise.
inline int MassBinIndex(const std::string& s) {
    if (s == "50-80GeV") return 0;
    if (s == "80-110GeV") return 1;
    if (s == "MassGE110GeV") return 2;
    // Also accept the display form of the labels
    if (s == "#geq110 GeV") return 2;
    if (s == "50-80 GeV") return 0;
    if (s == "80-110 GeV") return 1;
    return -1;
}

// Mass-cut label shown on the plots.
inline std::string MassBinDisplayLabel(int idx) {
    switch (idx) {
    case 0:
        return "50 < M^{J} < 80 GeV";
    case 1:
        return "80 < M^{J} < 110 GeV";
    case 2:
        return "M^{J} #geq 110 GeV";
    default:
        return "??? GeV";
    }
}

// What to draw on one combined plot: the y-axis title, where to save it, and
// functions returning the value and error from a CSV row. With
// drawComparison, a second error (getYErrAlt) is drawn next to the first.
struct ChartData {
    std::string yTitle;
    std::string outputPath;
    std::string variantLabel;
    std::string wpLabel;
    std::function<double(const Row&)> getY;
    std::function<double(const Row&)> getYErr;
    bool drawComparison = false; // also draw the error from getYErrAlt
    std::function<double(const Row&)> getYErrAlt;
};

// Draw one combined plot: DB and MJB points vs pT, one graph per channel and
// mass cut.
inline void DrawCombined(const std::vector<Row>& dbRows, const std::vector<Row>& mjbRows,
                         const ChartData& cd) {
    AtlasUtils::LoadAtlasStyle();

    TCanvas* c = new TCanvas(("c_" + cd.outputPath).c_str(), cd.wpLabel.c_str(), 900, 650);
    c->SetLeftMargin(0.13);
    c->SetRightMargin(0.04);
    c->SetTopMargin(0.08);
    c->SetBottomMargin(0.13);

    // y range from the points of both channels
    double ymin = 1e6, ymax = -1e6;
    for (const auto* rows : {&dbRows, &mjbRows}) {
        for (const auto& r : *rows) {
            double y = cd.getY(r);
            double e = cd.getYErr(r);
            if (y - e < ymin) ymin = y - e;
            if (y + e > ymax) ymax = y + e;
            if (cd.drawComparison) {
                double ea = cd.getYErrAlt(r);
                if (y - ea < ymin) ymin = y - ea;
                if (y + ea > ymax) ymax = y + ea;
            }
        }
    }
    double range = ymax - ymin;
    double pad = 0.30 * range;
    double ylo = ymin - pad;
    double yhi = ymax + pad;
    if (ylo > 0.92) ylo = 0.92;
    if (yhi < 1.10) yhi = 1.10;

    // x axis over the full DB + MJB range
    double xlo = 200.0, xhi = 1200.0;

    TH1F* frame = new TH1F(("fr_" + cd.outputPath).c_str(), "", 100, xlo, xhi);
    frame->SetMinimum(ylo);
    frame->SetMaximum(yhi);
    frame->SetStats(0);
    frame->GetXaxis()->SetTitle("p_{T} [GeV]");
    frame->GetYaxis()->SetTitle(cd.yTitle.c_str());
    frame->GetXaxis()->SetTitleSize(0.045);
    frame->GetYaxis()->SetTitleSize(0.045);
    frame->GetXaxis()->SetLabelSize(0.04);
    frame->GetYaxis()->SetLabelSize(0.04);
    frame->GetYaxis()->SetTitleOffset(1.35);
    frame->GetXaxis()->SetTitleOffset(1.20);
    frame->Draw("AXIS");

    // Line at 1
    TLine* refl = new TLine(xlo, 1.0, xhi, 1.0);
    refl->SetLineStyle(2);
    refl->SetLineColor(kBlack);
    refl->SetLineWidth(1);
    refl->Draw("same");

    // Line at 450 GeV between DB and MJB
    TLine* sep = new TLine(450.0, ylo, 450.0, yhi);
    sep->SetLineStyle(7);
    sep->SetLineColor(kGray + 2);
    sep->SetLineWidth(1);
    sep->Draw("same");

    // Channel names near the top of each pT range
    TLatex region;
    region.SetTextFont(42);
    region.SetTextSize(0.035);
    region.SetTextAlign(22);
    region.SetTextColor(kGray + 2);
    region.DrawLatex(0.5 * (xlo + 450.0), yhi - 0.03 * (yhi - ylo), "Z(#rightarrow b#bar{b}) + #gamma");
    region.DrawLatex(0.5 * (450.0 + xhi), yhi - 0.03 * (yhi - ylo), "Z(#rightarrow b#bar{b}) + jets");

    // Channel marker shapes
    constexpr int kMarkerDB = 20;  // closed circle
    constexpr int kMarkerMJB = 22; // closed triangle

    // Shift each mass cut a little in pT so the markers do not overlap:
    // -3 %, 0, +3 % of the channel's pT range.
    auto offsetForMassCut = [&](int massIdx, bool isDB) -> double {
        double regionWidth = isDB ? (450.0 - xlo) : (xhi - 450.0);
        double frac[3] = {-0.03, 0.0, +0.03};
        return regionWidth * frac[std::min(std::max(massIdx, 0), 2)];
    };

    // One graph per channel and mass cut. The graphs are kept in "alive" so
    // they exist until the canvas is saved.
    std::vector<TGraphErrors*> alive;

    auto plotChannel = [&](const std::vector<Row>& rows, bool isDB) {
        std::map<int, std::vector<const Row*>> byMass;
        for (const auto& r : rows) {
            int mi = MassBinIndex(r.mass_bin);
            if (mi >= 0) byMass[mi].push_back(&r);
        }
        for (auto& [mi, ptrs] : byMass) {
            if (ptrs.empty()) continue;
            int n = (int)ptrs.size();
            std::vector<double> x(n), y(n), ex(n, 0.0), ey(n);
            std::vector<double> ya(n); // alt error (for comparison plot)
            for (int i = 0; i < n; ++i) {
                x[i] = ptrs[i]->pt_center + offsetForMassCut(mi, isDB);
                y[i] = cd.getY(*ptrs[i]);
                ey[i] = std::abs(cd.getYErr(*ptrs[i]));
                if (cd.drawComparison) ya[i] = std::abs(cd.getYErrAlt(*ptrs[i]));
            }
            int marker = isDB ? kMarkerDB : kMarkerMJB;
            int color = MassCutColor(mi);

            auto* g = new TGraphErrors(n, x.data(), y.data(), ex.data(), ey.data());
            g->SetMarkerStyle(marker);
            g->SetMarkerColor(color);
            g->SetLineColor(color);
            g->SetMarkerSize(1.4);
            g->SetLineWidth(2);
            g->Draw("P same");
            alive.push_back(g);

            if (cd.drawComparison) {
                // Second error: dashed error bars, shifted a little, no marker
                std::vector<double> x2 = x;
                for (auto& xi : x2) xi += 0.5 * offsetForMassCut(mi, isDB);
                auto* gAlt = new TGraphErrors(n, x2.data(), y.data(), ex.data(), ya.data());
                gAlt->SetMarkerStyle(marker);
                gAlt->SetMarkerColor(color);
                gAlt->SetMarkerSize(0.001); // hide marker
                gAlt->SetLineColor(color);
                gAlt->SetLineStyle(2);
                gAlt->SetLineWidth(2);
                gAlt->Draw("|| same");
                alive.push_back(gAlt);
            }
        }
    };
    plotChannel(dbRows, /*isDB=*/true);
    plotChannel(mjbRows, /*isDB=*/false);

    // Mass-cut legend (colours) at the bottom left, away from the channel
    // names at the top.
    TLegend* legM = new TLegend(0.16, 0.16, 0.40, 0.36);
    legM->SetBorderSize(0);
    legM->SetFillStyle(0);
    legM->SetTextSize(0.030);
    for (int mi = 0; mi < 3; ++mi) {
        auto* gD = new TGraphErrors(1);
        gD->SetMarkerStyle(20);
        gD->SetMarkerColor(MassCutColor(mi));
        gD->SetLineColor(MassCutColor(mi));
        gD->SetMarkerSize(1.3);
        legM->AddEntry(gD, MassBinDisplayLabel(mi).c_str(), "lp");
        alive.push_back(gD);
    }
    legM->Draw();

    // Channel legend (marker shapes)
    TLegend* legC = new TLegend(0.42, 0.16, 0.62, 0.28);
    legC->SetBorderSize(0);
    legC->SetFillStyle(0);
    legC->SetTextSize(0.030);
    auto* gDBleg = new TGraphErrors(1);
    gDBleg->SetMarkerStyle(kMarkerDB);
    gDBleg->SetMarkerColor(kBlack);
    gDBleg->SetMarkerSize(1.4);
    auto* gMJBleg = new TGraphErrors(1);
    gMJBleg->SetMarkerStyle(kMarkerMJB);
    gMJBleg->SetMarkerColor(kBlack);
    gMJBleg->SetMarkerSize(1.4);
    legC->AddEntry(gDBleg, "DB (Z(#rightarrow b#bar{b})+#gamma)", "p");
    legC->AddEntry(gMJBleg, "MJB (Z(#rightarrow b#bar{b})+jets)", "p");
    alive.push_back(gDBleg);
    alive.push_back(gMJBleg);
    legC->Draw();

    // ATLAS label
    TLatex atlasLabel;
    atlasLabel.SetNDC();
    atlasLabel.SetTextFont(72);
    atlasLabel.SetTextSize(0.04);
    atlasLabel.DrawLatex(0.16, 0.85, "ATLAS");
    atlasLabel.SetTextFont(42);
    atlasLabel.DrawLatex(0.27, 0.85, "Internal");
    atlasLabel.SetTextSize(0.032);
    atlasLabel.DrawLatex(0.16, 0.80, RunConfig::LumiLabel().c_str());
    atlasLabel.DrawLatex(0.16, 0.76, ("Variant: " + cd.variantLabel + " / Nominal").c_str());
    atlasLabel.DrawLatex(0.16, 0.72, ("WP: " + cd.wpLabel).c_str());

    c->SaveAs(cd.outputPath.c_str());
    std::cout << "[Combined] Saved: " << cd.outputPath << "\n";

    delete c;
}

// The 7 combined plots of one correction and WP (list in the file header).
inline void RunForVariantWP(const std::string& variantLabel, const std::string& wp_file,
                            const std::string& wp_label, const std::vector<Row>& dbRows,
                            const std::vector<Row>& mjbRows, const std::string& outDir) {
    std::string vl = variantLabel;
    std::string suffix = "_" + vl + "_" + wp_file + "_Combined";
    auto path = [&](const std::string& stem) {
        return outDir + stem + suffix + ".pdf";
    };

    // 1) R_JES, uncorrelated errors (the only ones drawn here for R_JES)
    {
        ChartData cd;
        cd.yTitle = "R_{JES}^{Corrected} / R_{JES}^{Nominal}";
        cd.outputPath = path("R_JES_DR_vs_pT");
        cd.variantLabel = vl;
        cd.wpLabel = wp_label;
        cd.getY = [](const Row& r) {
            return r.rjes_dr;
        };
        cd.getYErr = [](const Row& r) {
            return r.rjes_dr_err;
        };
        DrawCombined(dbRows, mjbRows, cd);
    }

    // 2) sigma_insitu, uncorrelated errors
    {
        ChartData cd;
        cd.yTitle = "#sigma_{insitu}^{Corrected} / #sigma_{insitu}^{Nominal}";
        cd.outputPath = path("sigma_insitu_DR_vs_pT");
        cd.variantLabel = vl;
        cd.wpLabel = wp_label;
        cd.getY = [](const Row& r) {
            return r.sig_insitu_dr;
        };
        cd.getYErr = [](const Row& r) {
            return r.sig_insitu_dr_err_unc;
        };
        DrawCombined(dbRows, mjbRows, cd);
    }
    // 3) sigma_insitu, correlation-aware errors
    {
        ChartData cd;
        cd.yTitle = "#sigma_{insitu}^{Corrected} / #sigma_{insitu}^{Nominal}";
        cd.outputPath = outDir + "sigma_insitu_DR_vs_pT" + suffix + "_CorrAware.pdf";
        cd.variantLabel = vl;
        cd.wpLabel = wp_label;
        cd.getY = [](const Row& r) {
            return r.sig_insitu_dr;
        };
        cd.getYErr = [](const Row& r) {
            return r.sig_insitu_dr_err_corr;
        };
        DrawCombined(dbRows, mjbRows, cd);
    }
    // 4) sigma_insitu, both errors
    {
        ChartData cd;
        cd.yTitle = "#sigma_{insitu}^{Corrected} / #sigma_{insitu}^{Nominal}";
        cd.outputPath = outDir + "sigma_insitu_DR_vs_pT" + suffix + "_Cmp.pdf";
        cd.variantLabel = vl;
        cd.wpLabel = wp_label;
        cd.getY = [](const Row& r) {
            return r.sig_insitu_dr;
        };
        cd.getYErr = [](const Row& r) {
            return r.sig_insitu_dr_err_unc;
        };
        cd.drawComparison = true;
        cd.getYErrAlt = [](const Row& r) {
            return r.sig_insitu_dr_err_corr;
        };
        DrawCombined(dbRows, mjbRows, cd);
    }

    // 5) R_JER, uncorrelated errors
    {
        ChartData cd;
        cd.yTitle = "R_{JER}^{Corrected} / R_{JER}^{Nominal}";
        cd.outputPath = path("R_JER_DR_vs_pT");
        cd.variantLabel = vl;
        cd.wpLabel = wp_label;
        cd.getY = [](const Row& r) {
            return r.rjer_dr;
        };
        cd.getYErr = [](const Row& r) {
            return r.rjer_dr_err_unc;
        };
        DrawCombined(dbRows, mjbRows, cd);
    }
    // 6) R_JER, correlation-aware errors
    {
        ChartData cd;
        cd.yTitle = "R_{JER}^{Corrected} / R_{JER}^{Nominal}";
        cd.outputPath = outDir + "R_JER_DR_vs_pT" + suffix + "_CorrAware.pdf";
        cd.variantLabel = vl;
        cd.wpLabel = wp_label;
        cd.getY = [](const Row& r) {
            return r.rjer_dr;
        };
        cd.getYErr = [](const Row& r) {
            return r.rjer_dr_err_corr;
        };
        DrawCombined(dbRows, mjbRows, cd);
    }
    // 7) R_JER, both errors
    {
        ChartData cd;
        cd.yTitle = "R_{JER}^{Corrected} / R_{JER}^{Nominal}";
        cd.outputPath = outDir + "R_JER_DR_vs_pT" + suffix + "_Cmp.pdf";
        cd.variantLabel = vl;
        cd.wpLabel = wp_label;
        cd.getY = [](const Row& r) {
            return r.rjer_dr;
        };
        cd.getYErr = [](const Row& r) {
            return r.rjer_dr_err_unc;
        };
        cd.drawComparison = true;
        cd.getYErrAlt = [](const Row& r) {
            return r.rjer_dr_err_corr;
        };
        DrawCombined(dbRows, mjbRows, cd);
    }
}

} // namespace RatioResponseCombined


// Called by Response_Ratios_Combined.cxx in each channel's
// macros/ folder (usage in the file header).
//   variantOrAll = "ALL"        every correction with a CSV in both channels
//                  "bJER10v00"  only this correction
inline void Ratio_Response_Combined_Run(const std::string& variantOrAll = "ALL") {
    using namespace RatioResponseCombined;
    AtlasUtils::LoadAtlasStyle();

    auto wps = WPConfig::DefaultWPs();
    for (const auto& wp : wps) {
        std::string wp_file = PathBuilder::WPtoFileWP(wp);
        std::string wp_label = PathBuilder::WPtoLabel(wp);

        std::vector<VariantPair> pairs;
        if (variantOrAll == "ALL" || variantOrAll == "all") {
            pairs = FindCommonVariants(wp_file);
            if (pairs.empty()) {
                std::cout << "[Combined] WP " << wp_file << ": no variants found in BOTH channels.\n"
                          << "  Searched DB CSVs at:  " << CSVPath("DB", "DB_<key>", wp_file) << "\n"
                          << "  Searched MJB CSVs at: " << CSVPath("MJB", "MJB_<key>", wp_file) << "\n"
                          << "  Run from <channel>/macros/ where ../../<other-channel>/<plots>/ "
                          << "is a sibling of the current channel.\n";
                continue;
            }
        } else {
            // One correction requested
            VariantPair vp;
            vp.suffix = variantOrAll;
            vp.db_key = "DB_" + variantOrAll;
            vp.mjb_key = "MJB_" + variantOrAll;
            pairs.push_back(vp);
        }

        // Output in the Run{2,3}Plots of the channel the macro is run from.
        std::string outDir = (RunConfig::PlotDir() + "/Response/Combined_Ratios/") + wp_file + "/";
        gSystem->Exec(("mkdir -p " + outDir).c_str());

        for (const auto& vp : pairs) {
            // Read the CSV of both channels
            std::vector<Row> dbAll = ReadCSV(CSVPath("DB", vp.db_key, wp_file));
            std::vector<Row> mjbAll = ReadCSV(CSVPath("MJB", vp.mjb_key, wp_file));
            if (dbAll.empty() && mjbAll.empty()) {
                std::cerr << "[Combined] WARN: no rows for variant '" << vp.suffix << "' WP " << wp_file
                          << " — skipping.\n";
                continue;
            }
            // Keep the rows of this variant (the variant column holds the
            // full key, e.g. DB_bJER10v00).
            std::vector<Row> dbRows = RowsForVariant(dbAll, vp.db_key, wp_file);
            std::vector<Row> mjbRows = RowsForVariant(mjbAll, vp.mjb_key, wp_file);
            // If no row matches but the file is not empty, use the whole file
            // (in case the variant column is written differently).
            if (dbRows.empty() && !dbAll.empty()) dbRows = dbAll;
            if (mjbRows.empty() && !mjbAll.empty()) mjbRows = mjbAll;

            std::cout << "[Combined] " << vp.suffix << " WP " << wp_file << " — DB rows: " << dbRows.size()
                      << ", MJB rows: " << mjbRows.size() << "\n";
            RunForVariantWP(vp.suffix, wp_file, wp_label, dbRows, mjbRows, outDir);
        }
    }
    std::cout << "[Combined] Done.\n";
}

#endif // RATIO_RESPONSE_COMBINED_H
