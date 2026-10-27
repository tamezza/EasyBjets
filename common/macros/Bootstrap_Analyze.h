/**
 * ============================================================================
 *  Bootstrap_Analyze.h - double-ratio uncertainties and diagnostics from the
 *  bootstrap replicas
 * ============================================================================
 * ===========================================
 * @author:    Thiziri AMEZZA
 * @email:     thiziri.amezza@cern.ch
 * @date:      26 October 2026
 * ===========================================
 *
 * Bootstrap_Correlations writes, for every correction, WP and bin, the fitted
 * values of each replica to ../Run3Plots/Bootstrap/diagnostics/ (../Run2Plots
 * for run2): mu and sigma of data and MC, correction (v) and nominal (n), and
 * for the response the MC truth width sigma_t. This macro computes from them,
 * in every replica, the double ratio correction / nominal of each observable:
 *
 *   observable   double ratio DR                                   inputs
 *   b-JMDP       (mu_d / mu_m)_v / (mu_d / mu_m)_n                 4
 *   b-JMDR       (sig_d / sig_m)_v / (sig_d / sig_m)_n             4
 *   b-JES        (mu_d / mu_m)_v / (mu_d / mu_m)_n  (response R)   4
 *   sigma_insitu sqrt(sig_d^2 - sig_t^2)_v / (same)_n              4
 *   b-JER        R_JER,v / R_JER,n,
 *                R_JER = sqrt(sig_d^2 - sig_t^2) / sqrt(sig_m^2 - sig_t^2)   6
 *
 * and gives its uncertainty in four ways, all from the same replicas:
 * 
 *   bootstrap      spread of DR over the replicas: half of the central 68%
 *                  interval (16th to 84th percentile), robust to the few
 *                  replicas where a fit goes far off; all correlations
 *                  included, no linear approximation: the reference.
 *                  The standard deviation is also given (CSV, summary).
 * 
 *   full           linear propagation with the bootstrap covariance of the
 *                  inputs (what the "full" errors of the ratio macros do)
 * 
 *   corr-aware     linear propagation, rho = 1 between correction and nominal
 *                  of each input, 0 otherwise
 * 
 *   uncorrelated   linear propagation, all inputs independent
 * 
 * If "full" agrees with "bootstrap", the linear "full" error of the ratio
 * macros is right for that cell.
 *
 * Outputs in ../Run3Plots/Bootstrap/analysis/:
 *   <corr>_<wp>_<mass|resp>_<bin>/
 *     dr_<mass|resp>.pdf                 DR distribution of each observable, with
 *                                        the four uncertainties
 *     variance_decomposition_<obs>.pdf   share of the DR variance from each input
 *     corr_matrix_<mass|resp>.pdf        correlations of the fitted values
 *     marginals_<mass|resp>.pdf          distribution of each fitted value
 *     summary.txt                        the numbers of the above
 *   dr_errors.csv                        one row per cell and observable
 *                                        (relative uncertainties, bootstrap 68%
 *                                        and std, full, corr-aware, uncorrelated)
 *   dr_errors_<corr>_<obs>.pdf           relative DR uncertainty per bin, four ways
 *
 * Shared by both channels: ZbbGamma/macros/Bootstrap_Analyze.cxx and
 * ZbbJets/macros/Bootstrap_Analyze.cxx call Bootstrap_Analyze_Run(variant,
 * "DB" or "MJB").
 *
 * How to run (from ZbbGamma/macros or ZbbJets/macros, after Bootstrap_Correlations.cxx):
 *   export EASYBJETS_RUN=run3                              same era as the bootstrap
 *   root 'Bootstrap_Analyze.cxx()'                all corrections of the channel
 *   root 'Bootstrap_Analyze.cxx("DB_bJER10v01")'  one correction (MJB_... in ZbbJets)
 *
 * ==========================================================================
 * @version:  1.0
 * ==========================================================================
 */

#ifndef EASYBJETS_BOOTSTRAP_ANALYZE_H
#define EASYBJETS_BOOTSTRAP_ANALYZE_H

#include "RunEra.h" 

#include <TCanvas.h>
#include <TGraph.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TPad.h>
#include <TStyle.h>
#include <TBox.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TSystemDirectory.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace BootstrapAnalyze {

const double kNaN = std::numeric_limits<double>::quiet_NaN();

// ----------------------------------------------------------------------------
// Fitted values of one observable (mass or resp), one row per replica.
// Columns, as written by Bootstrap_Correlations (v = correction, n = nominal;
// v columns are even, their nominal partner is the next odd column):
//   0 mu_d_v   1 mu_d_n   2 mu_m_v   3 mu_m_n
//   4 sig_d_v  5 sig_d_n  6 sig_m_v  7 sig_m_n
//   8 sig_t_v  9 sig_t_n      (response only: MC truth width, TRUTH block)
// ----------------------------------------------------------------------------
struct Block {
    std::string name;                       // "mass" or "resp"
    std::vector<std::vector<double>> rows;  // rows[replica][column]
    bool empty() const { return rows.empty(); }
    int ncols() const { return rows.empty() ? 0 : (int)rows[0].size(); }
};

inline std::string ColName(int c) {  // for plots
    static const char* n[10] = {"#mu_{d}^{v}",    "#mu_{d}^{n}",    "#mu_{m}^{v}",    "#mu_{m}^{n}",
                                "#sigma_{d}^{v}", "#sigma_{d}^{n}", "#sigma_{m}^{v}", "#sigma_{m}^{n}",
                                "#sigma_{t}^{v}", "#sigma_{t}^{n}"};
    return (c >= 0 && c < 10) ? n[c] : "?";
}
inline std::string ColPlain(int c) {  // for text and CSV
    static const char* n[10] = {"mu_d_v",  "mu_d_n",  "mu_m_v",  "mu_m_n",  "sig_d_v",
                                "sig_d_n", "sig_m_v", "sig_m_n", "sig_t_v", "sig_t_n"};
    return (c >= 0 && c < 10) ? n[c] : "?";
}
inline std::string ColLong(int c) {  // for the variance bars
    static const char* n[10] = {"data peak, correction",   "data peak, nominal",
                                "MC peak, correction",     "MC peak, nominal",
                                "data width, correction",  "data width, nominal",
                                "MC width, correction",    "MC width, nominal",
                                "truth width, correction", "truth width, nominal"};
    return ColName(c) + "  " + ((c >= 0 && c < 10) ? n[c] : "?");
}

struct DiagnosticFile {
    std::string filename;  // bootstrap_<variant>_<wp>_<mass|resp>_<bin>.txt
    std::string variant;   // parsed from the file name
    std::string wp;
    std::string domain;    // "mass" or "resp" 
    std::string ptbin;     // bin label, may contain "_" (e.g. 200to300GeV_50-80GeV)
    Block mass;
    Block resp;
};

// ----------------------------------------------------------------------------
// Reads a diagnostic file. Blocks start with the comment lines
//   "# MASS block: ..."   8 columns
//   "# RESP block: ..."   8 columns
//   "# TRUTH block: ..."  2 columns (sig_t_v sig_t_n), same replicas as RESP;
//                         appended to the RESP rows as columns 8 and 9
// ----------------------------------------------------------------------------
inline DiagnosticFile ReadDiagnostic(const std::string& path) {
    DiagnosticFile d;
    d.filename = path;
    d.mass.name = "mass";
    d.resp.name = "resp";
    std::vector<std::vector<double>> truth;
    std::ifstream in(path);
    if (!in.is_open()) {
        std::cerr << "[Analyze] WARN: cannot open " << path << "\n";
        return d;
    }
    int block = 0;  // 1 mass, 2 resp, 3 truth
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        size_t i = 0;
        while (i < line.size() && std::isspace((unsigned char)line[i])) ++i;
        if (i < line.size() && line[i] == '#') {
            if (line.find("MASS block") != std::string::npos) block = 1;
            else if (line.find("RESP block") != std::string::npos) block = 2;
            else if (line.find("TRUTH block") != std::string::npos) block = 3;
            continue;
        }
        std::vector<double> row;
        std::stringstream ss(line);
        double v;
        while (ss >> v) row.push_back(v);
        if (block == 1 && row.size() == 8) d.mass.rows.push_back(row);
        if (block == 2 && row.size() == 8) d.resp.rows.push_back(row);
        if (block == 3 && row.size() == 2) truth.push_back(row);
    }
    if (!truth.empty()) {
        if (truth.size() == d.resp.rows.size()) {
            for (size_t r = 0; r < truth.size(); ++r) {
                d.resp.rows[r].push_back(truth[r][0]);
                d.resp.rows[r].push_back(truth[r][1]);
            }
        } else {
            std::cerr << "[Analyze] WARN: " << path << ": " << truth.size() << " truth rows for "
                      << d.resp.rows.size() << " response rows, truth not used\n";
        }
    }
    return d;
}

// ----------------------------------------------------------------------------
// Reads the correction, WP, observable (mass or response) and bin from the name of
// a diagnostic file:
//   bootstrap_<correction>_<wp>_<mass|resp>_<bin>.txt
//   e.g. bootstrap_DB_bJER10v01_0p46_resp_200to300GeV_50-80GeV.txt
//        -> correction DB_bJER10v01, WP 0p46, observable resp,
//           bin 200to300GeV_50-80GeV
// ----------------------------------------------------------------------------
inline void ParseFilename(DiagnosticFile& d) {
    auto slash = d.filename.find_last_of("/\\");
    std::string base = (slash == std::string::npos) ? d.filename : d.filename.substr(slash + 1);
    if (base.rfind("bootstrap_", 0) == 0) base = base.substr(10);
    if (base.size() > 4 && base.substr(base.size() - 4) == ".txt")
        base = base.substr(0, base.size() - 4);

    std::string head;  // "<variant>_<wp>"
    size_t pos = std::string::npos;
    for (const std::string dom : {"mass", "resp"}) {
        size_t p = base.find("_" + dom + "_");
        if (p != std::string::npos && (pos == std::string::npos || p < pos)) {
            pos = p;
            d.domain = dom;
        }
    }
    if (pos != std::string::npos) {
        head = base.substr(0, pos);
        d.ptbin = base.substr(pos + d.domain.size() + 2);
    } else {
        auto last = base.find_last_of('_');
        if (last == std::string::npos) {
            d.variant = base;
            return;
        }
        d.ptbin = base.substr(last + 1);
        head = base.substr(0, last);
    }
    auto last = head.find_last_of('_');
    if (last == std::string::npos) {
        d.variant = head;
        return;
    }
    d.wp = head.substr(last + 1);
    d.variant = head.substr(0, last);
}

// ----------------------------------------------------------------------------
// Statistics of the columns over the replicas.
// ----------------------------------------------------------------------------
struct Marginal {
    double mean = 0, std = 0, lo = 0, hi = 0;
    int n = 0;
};

inline Marginal StatsOf(const std::vector<double>& v) {
    Marginal m;
    if (v.empty()) return m;
    m.n = (int)v.size();
    double sum = 0, sumsq = 0;
    m.lo = 1e300;
    m.hi = -1e300;
    for (double x : v) {
        sum += x;
        sumsq += x * x;
        m.lo = std::min(m.lo, x);
        m.hi = std::max(m.hi, x);
    }
    m.mean = sum / m.n;
    m.std = std::sqrt(std::max(0.0, sumsq / m.n - m.mean * m.mean));
    return m;
}

inline Marginal Stats(const Block& b, int col) {
    std::vector<double> v;
    for (const auto& row : b.rows) v.push_back(row[col]);
    return StatsOf(v);
}

// Covariance matrix of all columns.
inline std::vector<std::vector<double>> CovMatrix(const Block& b) {
    const int n = b.ncols();
    std::vector<double> mean(n, 0.0);
    for (const auto& row : b.rows)
        for (int c = 0; c < n; ++c) mean[c] += row[c] / b.rows.size();
    std::vector<std::vector<double>> C(n, std::vector<double>(n, 0.0));
    for (const auto& row : b.rows)
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j) C[i][j] += (row[i] - mean[i]) * (row[j] - mean[j]) / b.rows.size();
    return C;
}

// Correlation matrix of all columns.
inline std::vector<std::vector<double>> CorrMatrix(const Block& b) {
    auto C = CovMatrix(b);
    const int n = (int)C.size();
    std::vector<std::vector<double>> M(n, std::vector<double>(n, 0.0));
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) {
            const double d = std::sqrt(C[i][i] * C[j][j]);
            M[i][j] = (d > 0) ? C[i][j] / d : 0.0;
        }
    return M;
}

// ----------------------------------------------------------------------------
// The observables: double ratio correction / nominal as a function of one
// replica row. Returns NaN when it cannot be computed (e.g. sigma_d < sigma_t).
// ----------------------------------------------------------------------------
struct Observable {
    std::string key;          // file tag: bJMDP, bJMDR, bJES, sigma_insitu, bJER
    std::string label;        // for plots, e.g. "b-JMDP"
    std::vector<int> inputs;  // columns it depends on
    std::function<double(const std::vector<double>&)> dr;
};

inline double Insitu(double sd, double st) {  // sqrt(sigma_d^2 - sigma_t^2)
    const double u = sd * sd - st * st;
    return (u > 0) ? std::sqrt(u) : kNaN;
}

inline double Ratio(double a, double b) { return (std::abs(b) > 1e-12) ? a / b : kNaN; }

inline std::vector<Observable> ObservablesFor(const Block& b) {
    std::vector<Observable> obs;
    auto peakDR = [](const std::vector<double>& r) { return Ratio(Ratio(r[0], r[2]), Ratio(r[1], r[3])); };
    if (b.name == "mass") {
        obs.push_back({"bJMDP", "b-JMDP", {0, 1, 2, 3}, peakDR});
        obs.push_back({"bJMDR", "b-JMDR", {4, 5, 6, 7}, [](const std::vector<double>& r) {
                           return Ratio(Ratio(r[4], r[6]), Ratio(r[5], r[7]));
                       }});
    } else {
        obs.push_back({"bJES", "b-JES", {0, 1, 2, 3}, peakDR});
        if (b.ncols() >= 10) {
            obs.push_back({"sigma_insitu", "#sigma_{insitu}", {4, 5, 8, 9}, [](const std::vector<double>& r) {
                               return Ratio(Insitu(r[4], r[8]), Insitu(r[5], r[9]));
                           }});
            obs.push_back({"bJER", "b-JER", {4, 5, 6, 7, 8, 9}, [](const std::vector<double>& r) {
                               const double v = Ratio(Insitu(r[4], r[8]), Insitu(r[6], r[8]));
                               const double n = Ratio(Insitu(r[5], r[9]), Insitu(r[7], r[9]));
                               return Ratio(v, n);
                           }});
        }
    }
    return obs;
}

// ----------------------------------------------------------------------------
// The four uncertainties of one observable's DR, all from the replicas.
//   boot    half of the 16-84% interval of DR over the replicas (robust)
//   bootStd standard deviation of DR over the replicas
//   full    sqrt(c^T Cov c): gradient c at the mean inputs, bootstrap covariance
//   corr    same, rho = 1 between correction and nominal of each input, 0 else
//   unc     same, inputs independent
// ----------------------------------------------------------------------------
struct DRErrors {
    int n = 0;
    double mean = kNaN;     // median of DR over the replicas
    double boot = kNaN;     // half of the 16-84% interval
    double bootStd = kNaN;  // standard deviation
    double full = kNaN, corr = kNaN, unc = kNaN;
    std::vector<double> values;  // DR per replica
};

inline DRErrors ComputeErrors(const Block& b, const Observable& o) {
    DRErrors e;
    for (const auto& row : b.rows) {
        const double x = o.dr(row);
        if (std::isfinite(x)) e.values.push_back(x);
    }
    e.n = (int)e.values.size();
    if (e.n < 10) return e;
    const Marginal m = StatsOf(e.values);
    std::vector<double> s = e.values;
    std::sort(s.begin(), s.end());
    auto quantile = [&](double q) {
        const double pos = q * (s.size() - 1);
        const size_t i = (size_t)pos;
        return (i + 1 < s.size()) ? s[i] + (pos - i) * (s[i + 1] - s[i]) : s.back();
    };
    e.mean = quantile(0.5);
    e.boot = 0.5 * (quantile(0.84) - quantile(0.16));
    e.bootStd = m.std;

    // gradient at the mean inputs (central differences)
    const auto C = CovMatrix(b);
    std::vector<double> mean(b.ncols(), 0.0);
    for (int c = 0; c < b.ncols(); ++c) mean[c] = Stats(b, c).mean;
    if (!std::isfinite(o.dr(mean))) return e;
    std::vector<double> grad(b.ncols(), 0.0);
    for (int c : o.inputs) {
        const double h = std::max(1e-3 * std::sqrt(C[c][c]), 1e-9 * std::abs(mean[c]) + 1e-12);
        std::vector<double> up = mean, dn = mean;
        up[c] += h;
        dn[c] -= h;
        grad[c] = (o.dr(up) - o.dr(dn)) / (2 * h);
    }
    double vFull = 0, vCorr = 0, vUnc = 0;
    for (int i : o.inputs)
        for (int j : o.inputs) {
            const double gg = grad[i] * grad[j];
            vFull += gg * C[i][j];
            if (i == j) {
                vCorr += gg * C[i][i];
                vUnc += gg * C[i][i];
            } else if ((i ^ 1) == j) {  // correction and nominal of the same input
                vCorr += gg * std::sqrt(C[i][i] * C[j][j]);
            }
        }
    e.full = std::sqrt(std::max(0.0, vFull));
    e.corr = std::sqrt(std::max(0.0, vCorr));
    e.unc = std::sqrt(std::max(0.0, vUnc));
    return e;
}

// ----------------------------------------------------------------------------
// PLOT - correlation matrix of the fitted values (8x8, 10x10 with the truth
// widths). Cell (x, y) = correlation of quantities x and y over the replicas.
// ----------------------------------------------------------------------------
inline void DrawCorrMatrix(const Block& b, const std::string& title, const std::string& outPath) {
    if (b.empty()) return;
    const auto M = CorrMatrix(b);
    const int n = b.ncols();

    TCanvas* c = new TCanvas(("c_corr_" + b.name).c_str(), title.c_str(), 820, 780);
    c->SetLeftMargin(0.17);
    c->SetRightMargin(0.21);  
    c->SetBottomMargin(0.19); 
    c->SetTopMargin(0.10);

    TH2F* h = new TH2F(("h_corr_" + b.name).c_str(), title.c_str(), n, 0, n, n, 0, n);
    for (int i = 0; i < n; ++i) {
        h->GetXaxis()->SetBinLabel(i + 1, ColName(i).c_str());
        h->GetYaxis()->SetBinLabel(n - i, ColName(i).c_str());  
        for (int j = 0; j < n; ++j) h->SetBinContent(i + 1, n - j, M[j][i]);
    }
    h->SetMinimum(-1.0);
    h->SetMaximum(1.0);
    h->SetStats(0);
    h->GetXaxis()->SetLabelSize(0.045);
    h->GetYaxis()->SetLabelSize(0.045);
    h->GetXaxis()->SetTickLength(0);
    h->GetYaxis()->SetTickLength(0);
    h->GetXaxis()->SetTitle("Fitted quantity");
    h->GetYaxis()->SetTitle("Fitted quantity");
    h->GetXaxis()->SetTitleSize(0.040);
    h->GetYaxis()->SetTitleSize(0.040);
    h->GetXaxis()->SetTitleOffset(1.35);
    h->GetYaxis()->SetTitleOffset(1.75);
    h->GetZaxis()->SetTitle("Correlation over the replicas");
    h->GetZaxis()->SetTitleSize(0.035);
    h->GetZaxis()->SetTitleOffset(1.55);
    gStyle->SetPalette(kBlueRedYellow);
    gStyle->SetPaintTextFormat("5.2f");
    h->SetMarkerSize(n > 8 ? 1.0 : 1.2);
    h->Draw("COLZ TEXT");

    TLatex leg;
    leg.SetNDC();
    leg.SetTextFont(42);
    leg.SetTextSize(0.026);
    leg.DrawLatex(0.17, 0.035, n > 8 ? "#mu, #sigma: fitted peak and width;  d = data, m = MC, t = MC truth;"
                                       "  v = correction, n = nominal"
                                     : "#mu, #sigma: fitted peak and width;  d = data, m = MC;"
                                       "  v = correction, n = nominal");
    c->SaveAs(outPath.c_str());
    delete c;
    delete h;
}

// ----------------------------------------------------------------------------
// PLOT - distribution of each fitted value over the replicas.
// ----------------------------------------------------------------------------
inline void DrawMarginals(const Block& b, const std::string& title, const std::string& outPath) {
    if (b.empty()) return;
    const int n = b.ncols();
    const int nx = (n + 1) / 2;

    TCanvas* c = new TCanvas(("c_marg_" + b.name).c_str(), title.c_str(), 300 * nx, 700);
    c->Divide(nx, 2, 0.001, 0.001);
    for (int col = 0; col < n; ++col) {
        c->cd(col + 1);
        gPad->SetLeftMargin(0.18);
        gPad->SetBottomMargin(0.15);
        gPad->SetTopMargin(0.10);

        const Marginal m = Stats(b, col);
        double xlo = m.mean - 5 * m.std;
        double xhi = m.mean + 5 * m.std;
        if (xhi - xlo < 1e-6) {
            xlo = m.lo - 0.01 * std::abs(m.lo + 1e-6);
            xhi = m.hi + 0.01 * std::abs(m.hi + 1e-6);
        }
        TH1F* h = new TH1F(Form("h_marg_%s_%d", b.name.c_str(), col),
                           Form("%s;%s;Replicas", ColName(col).c_str(), ColName(col).c_str()), 40, xlo, xhi);
        for (const auto& row : b.rows) h->Fill(row[col]);
        h->SetFillColorAlpha(kAzure - 4, 0.5);
        h->SetLineColor(kAzure + 2);
        h->SetLineWidth(2);
        h->SetStats(0);
        h->GetXaxis()->SetTitleSize(0.06);
        h->GetYaxis()->SetTitleSize(0.06);
        h->GetXaxis()->SetLabelSize(0.05);
        h->GetYaxis()->SetLabelSize(0.05);
        h->Draw("HIST");

        TLine* lm = new TLine(m.mean, 0, m.mean, h->GetMaximum());
        lm->SetLineColor(kRed);
        lm->SetLineWidth(2);
        lm->SetLineStyle(2);
        lm->Draw("same");
        TBox* band = new TBox(m.mean - m.std, 0, m.mean + m.std, h->GetMaximum() * 0.05);
        band->SetFillColorAlpha(kRed, 0.4);
        band->Draw("same");

        TLatex* t = new TLatex();
        t->SetNDC();
        t->SetTextSize(0.06);
        t->SetTextFont(42);
        t->DrawLatex(0.20, 0.83, Form("mean = %.4f", m.mean));
        t->DrawLatex(0.20, 0.76, Form("std = %.4f", m.std));
        t->DrawLatex(0.20, 0.69, Form("N = %d", m.n));
    }
    c->SaveAs(outPath.c_str());
    delete c;
}

// ----------------------------------------------------------------------------
// PLOT - DR of each observable over the replicas, one panel per observable,
// with the bootstrap spread and the three linear uncertainties.
// ----------------------------------------------------------------------------
inline void DrawDRDistributions(const Block& b, const std::vector<Observable>& obs,
                                const std::vector<DRErrors>& errs, const std::string& title,
                                const std::string& outPath) {
    if (b.empty() || obs.empty()) return;
    const int n = (int)obs.size();
    TCanvas* c = new TCanvas(("c_dr_" + b.name).c_str(), title.c_str(), 520 * n, 620);
    TPad* body = new TPad("body", "", 0, 0, 1, 0.93);  
    body->Draw();
    body->cd();
    body->Divide(n, 1, 0.001, 0.001);
    for (int k = 0; k < n; ++k) {
        body->cd(k + 1);
        gPad->SetLeftMargin(0.15);
        gPad->SetBottomMargin(0.13);
        gPad->SetTopMargin(0.32);
        const DRErrors& e = errs[k];
        if (e.n < 10) {
            TLatex warn;
            warn.SetNDC();
            warn.SetTextAlign(22);
            warn.DrawLatex(0.5, 0.5, "not enough valid replicas");
            continue;
        }
        const double w = std::max(5 * e.boot, 1e-6);  // range from the robust spread
        TH1F* h = new TH1F(Form("h_dr_%s_%s", b.name.c_str(), obs[k].key.c_str()),
                           Form(";DR(%s);Replicas", obs[k].label.c_str()), 40, e.mean - w, e.mean + w);
        for (double x : e.values) h->Fill(x);
        h->SetFillColorAlpha(kViolet - 3, 0.5);
        h->SetLineColor(kViolet + 1);
        h->SetLineWidth(2);
        h->SetStats(0);
        h->GetXaxis()->SetTitleSize(0.05);
        h->GetYaxis()->SetTitleSize(0.05);
        h->GetXaxis()->SetNdivisions(505);
        h->Draw("HIST");

        TLine* lm = new TLine(e.mean, 0, e.mean, h->GetMaximum());
        lm->SetLineColor(kRed);
        lm->SetLineWidth(2);
        lm->SetLineStyle(2);
        lm->Draw("same");
        TBox* band = new TBox(e.mean - e.boot, 0, e.mean + e.boot, h->GetMaximum() * 0.04);
        band->SetFillColorAlpha(kRed, 0.5);
        band->Draw("same");

        TLatex t;
        t.SetNDC();
        t.SetTextFont(42);
        t.SetTextSize(0.050);
        t.DrawLatex(0.15, 0.94, Form("DR(%s) = %.4f  (median)", obs[k].label.c_str(), e.mean));
        t.SetTextSize(0.040);
        auto rel = [&](double s) { return 100.0 * s / std::abs(e.mean); };
        t.DrawLatex(0.15, 0.885, Form("bootstrap 68%%:  #pm %.4f  (%.2f %%)", e.boot, rel(e.boot)));
        t.DrawLatex(0.15, 0.840, Form("bootstrap std:   #pm %.4f  (%.2f %%)", e.bootStd, rel(e.bootStd)));
        t.DrawLatex(0.15, 0.795, Form("full (linear):    #pm %.4f  (%.2f %%)", e.full, rel(e.full)));
        t.DrawLatex(0.15, 0.750, Form("corr-aware:      #pm %.4f  (%.2f %%)", e.corr, rel(e.corr)));
        t.DrawLatex(0.15, 0.705, Form("uncorrelated:   #pm %.4f  (%.2f %%)", e.unc, rel(e.unc)));
    }
    c->cd();
    TLatex tt;
    tt.SetNDC();
    tt.SetTextFont(42);
    tt.SetTextSize(0.035);
    tt.DrawLatex(0.01, 0.955, title.c_str());
    c->SaveAs(outPath.c_str());
    delete c;
}

// ----------------------------------------------------------------------------
// PLOT - variance decomposition of one observable, one bar per input.
//
// The DR of the observable is computed in every replica; its variance over
// the replicas is the total variance. For each input, DR is recomputed with
// only that input taking its replica values and the others fixed at their
// mean: the bar is this variance / the total variance.
//
// The bars add up to 100% only if the inputs fluctuate independently. The
// correction and the nominal are fitted on the same events and move together,
// which partly cancels in DR, so the bars usually add up to more than 100%;
// the ratio printed under the plot shows by how much.
// ----------------------------------------------------------------------------
inline void DrawVarianceDecomposition(const Block& b, const Observable& o, const DRErrors& e,
                                      const std::string& title, const std::string& outPath) {
    if (b.empty() || e.n < 10 || !(e.bootStd > 0)) return;
    const double full_var = e.bootStd * e.bootStd;  // variances
    std::vector<double> mean(b.ncols(), 0.0);
    for (int c = 0; c < b.ncols(); ++c) mean[c] = Stats(b, c).mean;

    const int ni = (int)o.inputs.size();
    std::vector<double> partial(ni, 0.0);
    for (int k = 0; k < ni; ++k) {
        std::vector<double> v;
        for (const auto& row : b.rows) {
            std::vector<double> r = mean;
            r[o.inputs[k]] = row[o.inputs[k]];
            const double x = o.dr(r);
            if (std::isfinite(x)) v.push_back(x);
        }
        const Marginal m = StatsOf(v);
        partial[k] = m.std * m.std;
    }
    double sum_partial = 0;
    for (double x : partial) sum_partial += x;

    const int cols[10] = {kRed - 4,   kAzure - 4, kGreen - 3, kOrange - 3, kMagenta - 7,
                          kCyan - 3,  kYellow - 6, kSpring - 5, kGray + 1, kPink - 4};
    TCanvas* c = new TCanvas(("c_var_" + o.key).c_str(), title.c_str(), 900, 360 + 70 * ni);
    c->SetLeftMargin(0.36);
    c->SetRightMargin(0.05);
    c->SetTopMargin(0.16);
    c->SetBottomMargin(0.26);

    TH1F* frame = new TH1F(("h_var_" + o.key).c_str(), "", ni, 0, ni);
    double ymax = 0;
    for (int k = 0; k < ni; ++k) {
        const double pct = 100.0 * partial[k] / full_var;
        frame->SetBinContent(ni - k, pct);  // first input on top
        frame->GetXaxis()->SetBinLabel(ni - k, ColLong(o.inputs[k]).c_str());
        ymax = std::max(ymax, pct);
    }
    frame->SetStats(0);
    frame->SetMinimum(0);
    frame->SetMaximum(std::max(1.0, ymax * 1.25));
    frame->GetYaxis()->SetTitle(Form("Share of the DR(%s) variance [%%]", o.label.c_str()));
    frame->GetXaxis()->SetTitle("Input that fluctuates alone");
    frame->GetXaxis()->SetLabelSize(0.040);
    frame->GetYaxis()->SetLabelSize(0.040);
    frame->GetXaxis()->SetTitleSize(0.042);
    frame->GetYaxis()->SetTitleSize(0.042);
    frame->GetXaxis()->SetTitleOffset(3.6);
    frame->GetYaxis()->SetTitleOffset(1.1);
    frame->SetBarWidth(0.6);
    frame->SetBarOffset(0.2);
    frame->SetFillColor(kGray);
    frame->Draw("HBAR");
    TLatex val;
    val.SetTextFont(42);
    val.SetTextSize(0.036);
    val.SetTextAlign(12);
    for (int k = 0; k < ni; ++k) {
        TH1F* hb = (TH1F*)frame->Clone(Form("h_var_%s_%d", o.key.c_str(), k));
        hb->Reset();
        hb->SetBinContent(ni - k, frame->GetBinContent(ni - k));
        hb->SetFillColor(cols[o.inputs[k] % 10]);
        hb->Draw("HBAR SAME");
        val.DrawLatex(frame->GetBinContent(ni - k) + 0.01 * frame->GetMaximum(), ni - k - 0.5,
                      Form("%.1f %%", frame->GetBinContent(ni - k)));
    }

    TLatex t;
    t.SetNDC();
    t.SetTextFont(42);
    t.SetTextSize(0.042);
    t.DrawLatex(0.05, 0.93, Form("DR(%s) variance decomposition", o.label.c_str()));
    t.SetTextSize(0.034);
    t.DrawLatex(0.05, 0.875, title.c_str());
    t.SetTextSize(0.032);
    t.DrawLatex(0.05, 0.075,
                Form("#sigma^{2}_{all} = %.2e,   #Sigma #sigma^{2}_{alone} = %.2e,   ratio = %.2f", full_var,
                     sum_partial, sum_partial / full_var));
    t.DrawLatex(0.05, 0.030,
                "Each bar: variance when only that input fluctuates. Ratio #neq 1: correlated inputs.");
    c->SaveAs(outPath.c_str());
    delete c;
}

// ----------------------------------------------------------------------------
// Per-cell summary text: fitted values, correlations, DR uncertainties.
// ----------------------------------------------------------------------------
inline void WriteSummaryText(const DiagnosticFile& d, const std::string& outPath) {
    std::ofstream out(outPath);
    out << "# Bootstrap summary\n";
    out << "# Correction: " << d.variant << "\n";
    out << "# WP:         " << d.wp << "\n";
    out << "# Bin:        " << d.ptbin << "\n\n";

    auto writeBlock = [&](const Block& b) {
        if (b.empty()) return;
        const int n = b.ncols();
        out << "===== " << b.name << " =====\n";
        out << "Replicas: " << b.rows.size() << "\n\n";
        out << std::left << std::setw(10) << "quantity" << std::right << std::setw(14) << "mean"
            << std::setw(14) << "std" << std::setw(14) << "min" << std::setw(14) << "max" << "\n";
        for (int c = 0; c < n; ++c) {
            const Marginal m = Stats(b, c);
            out << std::left << std::setw(10) << ColPlain(c) << std::right << std::fixed
                << std::setprecision(6) << std::setw(14) << m.mean << std::setw(14) << m.std
                << std::setw(14) << m.lo << std::setw(14) << m.hi << "\n";
        }
        out << "\nCorrelation matrix:\n" << std::setw(10) << "";
        for (int c = 0; c < n; ++c) out << std::setw(9) << ColPlain(c);
        out << "\n";
        const auto M = CorrMatrix(b);
        for (int i = 0; i < n; ++i) {
            out << std::left << std::setw(10) << ColPlain(i) << std::right;
            for (int j = 0; j < n; ++j) out << std::setw(9) << std::fixed << std::setprecision(3) << M[i][j];
            out << "\n";
        }
        out << "\nDouble-ratio uncertainties (relative, %):\n";
        out << std::left << std::setw(14) << "observable" << std::right << std::setw(10) << "DR"
            << std::setw(12) << "boot 68%" << std::setw(12) << "boot std" << std::setw(12) << "full"
            << std::setw(12) << "corr-aware"
            << std::setw(14) << "uncorrelated" << std::setw(10) << "N" << "\n";
        for (const auto& o : ObservablesFor(b)) {
            const DRErrors e = ComputeErrors(b, o);
            auto rel = [&](double s) { return 100.0 * s / std::abs(e.mean); };
            out << std::left << std::setw(14) << o.key << std::right << std::fixed << std::setprecision(4)
                << std::setw(10) << e.mean << std::setprecision(3) << std::setw(12) << rel(e.boot)
                << std::setw(12) << rel(e.bootStd) << std::setw(12) << rel(e.full) << std::setw(12)
                << rel(e.corr) << std::setw(14) << rel(e.unc) << std::setw(10) << e.n << "\n";
        }
        out << "\n";
    };
    writeBlock(d.mass);
    writeBlock(d.resp);
}

// ----------------------------------------------------------------------------
// One row of dr_errors.csv per cell and observable.
// ----------------------------------------------------------------------------
struct ErrorRow {
    std::string channel, observable, variant, wp, bin, key;
    DRErrors e;
};

inline void WriteErrorCSV(const std::vector<ErrorRow>& rows, const std::string& path) {
    std::ofstream out(path);
    out << "channel,observable,dr,correction,wp,bin,n_replicas,dr_median,"
           "err_bootstrap68,err_bootstrap_std,err_full,err_corr_aware,err_uncorrelated,"
           "rel_bootstrap68,rel_bootstrap_std,rel_full,rel_corr_aware,rel_uncorrelated\n";
    for (const auto& r : rows) {
        const double a = std::abs(r.e.mean);
        out << r.channel << "," << r.observable << "," << r.key << "," << r.variant << "," << r.wp << ","
            << r.bin << "," << r.e.n << std::fixed << std::setprecision(6) << "," << r.e.mean << ","
            << r.e.boot << "," << r.e.bootStd << "," << r.e.full << "," << r.e.corr << "," << r.e.unc << ","
            << r.e.boot / a << "," << r.e.bootStd / a << "," << r.e.full / a << "," << r.e.corr / a << ","
            << r.e.unc / a << "\n";
    }
}

// ----------------------------------------------------------------------------
// PLOT - relative DR uncertainty per bin of one correction and observable,
// the four ways side by side.
// ----------------------------------------------------------------------------
inline void DrawErrorOverview(const std::vector<ErrorRow>& rows, const std::string& label,
                              const std::string& title, const std::string& outPath) {
    const int n = (int)rows.size();
    if (n == 0) return;
    TCanvas* c = new TCanvas(("c_err_" + outPath).c_str(), title.c_str(), std::max(800, 110 * n + 300), 700);
    c->SetBottomMargin(0.36);  
    c->SetLeftMargin(0.12);
    c->SetRightMargin(0.04);
    c->SetTopMargin(0.20);

    double ymax = 0;
    for (const auto& r : rows) {
        const double a = std::abs(r.e.mean);
        for (double s : {r.e.boot, r.e.full, r.e.corr, r.e.unc})
            if (std::isfinite(s)) ymax = std::max(ymax, 100.0 * s / a);
    }
    TH1F* frame = new TH1F(("h_err_" + outPath).c_str(), "", n, 0, n);
    for (int i = 0; i < n; ++i) frame->GetXaxis()->SetBinLabel(i + 1, rows[i].bin.c_str());
    frame->SetMinimum(0);
    frame->SetMaximum(std::max(0.1, ymax * 1.30));
    frame->SetStats(0);
    frame->GetYaxis()->SetTitle(Form("Relative uncertainty of DR(%s) [%%]", label.c_str()));
    frame->GetYaxis()->SetTitleSize(0.045);
    frame->GetYaxis()->SetTitleOffset(1.1);
    frame->GetXaxis()->SetLabelSize(0.034);
    frame->GetXaxis()->LabelsOption("v");
    frame->Draw("AXIS");

    struct Way { const char* name; int col; int mk; double off; };
    const Way ways[4] = {{"bootstrap (68% interval of DR)", kBlack, 20, -0.24},
                         {"full (bootstrap #rho, linear)", kRed + 1, 21, -0.08},
                         {"corr-aware (#rho = 1)", kAzure + 2, 22, 0.08},
                         {"uncorrelated", kGray + 1, 24, 0.24}};
    TLegend* leg = new TLegend(0.12, 0.86, 0.96, 0.98);
    leg->SetNColumns(2);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.035);
    for (int w = 0; w < 4; ++w) {
        TGraph* g = new TGraph();
        for (int i = 0; i < n; ++i) {
            const DRErrors& e = rows[i].e;
            const double s = (w == 0) ? e.boot : (w == 1) ? e.full : (w == 2) ? e.corr : e.unc;
            if (std::isfinite(s)) g->SetPoint(g->GetN(), i + 0.5 + ways[w].off, 100.0 * s / std::abs(e.mean));
        }
        g->SetMarkerStyle(ways[w].mk);
        g->SetMarkerSize(1.4);
        g->SetMarkerColor(ways[w].col);
        g->Draw("P SAME");
        leg->AddEntry(g, ways[w].name, "p");
    }
    leg->Draw();
    TLatex t;
    t.SetNDC();
    t.SetTextFont(42);
    t.SetTextSize(0.035);
    t.DrawLatex(0.12, 0.815, title.c_str());
    c->SaveAs(outPath.c_str());
    delete c;
}

// ---------------------------------------------------------------------------------------------
// Folder of the bootstrap outputs of this channel and era (as written byBootstrap_Correlations).
// ---------------------------------------------------------------------------------------------
inline std::string BootstrapDir() {
    return std::string(RunConfig::IsRun3() ? "../Run3Plots" : "../Run2Plots") + "/Bootstrap/";
}

// -------------------------------------------------------------------------------------------
// Diagnostic files in <BootstrapDir>/diagnostics/, optionally of one variant, sorted by name.
// -------------------------------------------------------------------------------------------
inline std::vector<std::string> FindDiagnostics(const std::string& filter_variant) {
    std::vector<std::string> out;
    std::string dir = BootstrapDir() + "diagnostics/";
    TSystemDirectory d(dir.c_str(), dir.c_str());
    TList* files = d.GetListOfFiles();
    if (!files) {
        std::cerr << "[Analyze] WARN: cannot list " << dir << "\n";
        return out;
    }
    TIter next(files);
    while (TObject* obj = next()) {
        std::string fn = obj->GetName();
        if (fn.rfind("bootstrap_", 0) != 0) continue;
        if (fn.size() < 4 || fn.substr(fn.size() - 4) != ".txt") continue;
        const bool all = (filter_variant == "ALL" || filter_variant == "all");
        if (!all && fn.find(filter_variant) == std::string::npos) continue;
        out.push_back(dir + fn);
    }
    std::sort(out.begin(), out.end());
    return out;
}

}  // namespace BootstrapAnalyze

// =============================================================================
// Entry point, called by Bootstrap_Analyze.cxx in each channel's macros/ folder.
// variant: "ALL" or one correction key (e.g. DB_bJER10v01); channel: DB or MJB.
// =============================================================================
inline void Bootstrap_Analyze_Run(const std::string& variant, const std::string& channel) {
    using namespace BootstrapAnalyze;
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);
    TH1::AddDirectory(false);  // the per-cell histograms are re-created for every cell

    const std::vector<std::string> paths = FindDiagnostics(variant);
    if (paths.empty()) {
        std::cerr << "[Analyze] No diagnostic files found in " << BootstrapDir() << "diagnostics/"
                  << " for variant filter '" << variant << "'.\n"
                  << "Run Bootstrap_Correlations.cxx first, with the same EASYBJETS_RUN.\n";
        return;
    }
    std::cout << "[Analyze] Found " << paths.size() << " diagnostic file(s).\n";

    const std::string outRoot = BootstrapDir() + "analysis/";
    gSystem->Exec(("mkdir -p " + outRoot).c_str());

    std::vector<ErrorRow> errorRows;
    int processed = 0;
    for (const auto& path : paths) {
        DiagnosticFile d = ReadDiagnostic(path);
        ParseFilename(d);
        if (d.mass.empty() && d.resp.empty()) {
            std::cerr << "[Analyze] WARN: " << path << " has no readable rows\n";
            continue;
        }
        const std::string domTag = d.domain.empty() ? "" : d.domain + "_";
        const std::string cellTag = d.variant + "_" + d.wp + "_" + domTag + d.ptbin;
        const std::string cellDir = outRoot + cellTag + "/";
        gSystem->Exec(("mkdir -p " + cellDir).c_str());
        const std::string title = d.variant + ", WP " + d.wp + ", " + d.ptbin;

        for (const Block* b : {&d.mass, &d.resp}) {
            if (b->empty()) continue;
            const std::string tag = b->name;
            const std::string what = (tag == "mass") ? "Mass" : "Response";
            const auto obs = ObservablesFor(*b);
            std::vector<DRErrors> errs;
            for (const auto& o : obs) {
                errs.push_back(ComputeErrors(*b, o));
                errorRows.push_back({channel, tag, d.variant, d.wp, d.ptbin, o.key, errs.back()});
            }
            DrawDRDistributions(*b, obs, errs, what + " - " + title, cellDir + "dr_" + tag + ".pdf");
            for (size_t k = 0; k < obs.size(); ++k)
                DrawVarianceDecomposition(*b, obs[k], errs[k], what + " - " + title,
                                          cellDir + "variance_decomposition_" + obs[k].key + ".pdf");
            DrawCorrMatrix(*b, what + " - " + title, cellDir + "corr_matrix_" + tag + ".pdf");
            DrawMarginals(*b, what + " - " + title, cellDir + "marginals_" + tag + ".pdf");
        }
        WriteSummaryText(d, cellDir + "summary.txt");
        ++processed;
        std::cout << "[Analyze] " << cellTag << ": " << (d.mass.empty() ? 0 : (int)d.mass.rows.size())
                  << " mass + " << (d.resp.empty() ? 0 : (int)d.resp.rows.size()) << " response replicas"
                  << (d.resp.ncols() >= 10 ? " (with truth)" : "") << "\n";
    }

    WriteErrorCSV(errorRows, outRoot + "dr_errors.csv");

    // Overview per correction and observable: relative uncertainty per bin
    std::map<std::string, std::vector<ErrorRow>> groups;
    for (const auto& r : errorRows) groups[r.variant + "|" + r.key].push_back(r);
    for (const auto& kv : groups) {
        const ErrorRow& r0 = kv.second.front();
        std::string label = r0.key;
        if (label == "bJMDP") label = "b-JMDP";
        else if (label == "bJMDR") label = "b-JMDR";
        else if (label == "bJES") label = "b-JES";
        else if (label == "bJER") label = "b-JER";
        else if (label == "sigma_insitu") label = "#sigma_{insitu}";
        DrawErrorOverview(kv.second, label, r0.variant + ", WP " + r0.wp,
                          outRoot + "dr_errors_" + r0.variant + "_" + r0.key + ".pdf");
    }

    std::cout << "\n[Analyze] Done. " << processed << " cell(s) processed.\n"
              << "  Per cell:        " << outRoot << "<correction>_<wp>_<mass|resp>_<bin>/\n"
              << "  DR uncertainties: " << outRoot << "dr_errors.csv and dr_errors_<correction>_<obs>.pdf\n";
}

#endif  // EASYBJETS_BOOTSTRAP_ANALYZE_H
