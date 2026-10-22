/**
 * ============================================================================
 *  MassFit_BackgroundStability.C - background-modeling stability of the
 *  Zbb peak fit (Mass_Fit), per reviewer request
 * ============================================================================
 * ===========================================
 * @author:    Thiziri AMEZZA
 * @email:     thiziri.amezza@cern.ch
 * @date:      18 April 2026 
 * ===========================================
 * 
 * Mass_Fit fits the three background families (Exp, Poly, PolyExp) and
 * writes a signal_fit_params.txt in each function's folder. This macro
 * reads those files for every (variant, WP, pT bin) and reports how much
 * the data peak mu and width sigma, and hence b-JMDP = mu_data/mu_MC and
 * b-JMDR = sigma_data/sigma_MC, move between background functions.
 *
 *  Note on MC: h_mc is the Z->bb signal MC mass template, fit DIRECTLY with a
 *  Gaussian (NO background subtraction - the MC sample is pure signal). So
 *  mc_mu/mc_sigma are identical across the three background functions; the MC
 *  panel is flat (a cross-check that MC is unaffected by the background choice),
 *  and b-JMDP/b-JMDR vary across functions only through the data.
 *
 *  USAGE  (run from EasyBjets/CrossChecks/  - also works from EasyBjets/)
 *  -----
 *    root 'MassFit_BackgroundStability.C("run2","DB")'
 *    root 'MassFit_BackgroundStability.C("run3","MJB","MJB_bJER10v01")'
 *    root 'MassFit_BackgroundStability.C("run2","DB","ALL",2.0)'
 *
 *    era       : run2 | run3            (Run2Plots / Run3Plots)
 *    channel   : DB | MJB               (DB = ZbbGamma, MJB = ZbbJets)
 *    variant   : "ALL" (default) or a specific plot_subdir (e.g. "DB_bJER10v00")
 *    chi2Gate  : data chi2/ndf threshold for the quality-gated systematic (2.0)
 *
 *  OUTPUT  (under <topology>/Run{2,3}Plots/<channel>_BackgroundStability/)
 *  ------
 *    background_stability_summary.csv   - one row per (variant, WP, pT bin)
 *    background_systematics_table.txt   - the spreads with four definitions
 *    <cell>_data_stability.{pdf,png}    - data  mu/sigma vs function
 *    <cell>_mc_stability.{pdf,png}      - MC    mu/sigma vs function (flat)
 *    <cell>_datamc_stability.{pdf,png}  - data + MC overlaid (x-offset)
 *    <cell>_bJMD_stability.{pdf,png}    - b-JMDP / b-JMDR vs function
 * 
 * ==========================================================================
 * @version:  1.0
 * ==========================================================================
 */
#include <TSystem.h>
#include <TString.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TH1F.h>
#include <TGraphErrors.h>
#include <TLine.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TError.h>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <functional>
#include <algorithm>
#include <cmath>

namespace {

    // one signal_fit_params.txt
    struct FuncResult {
        std::string func;
        bool   is_optimal = false;
        double bkg_chi2_ndf = 0, spurious = 0;
        double data_mu = 0, data_mu_err = 0, data_sigma = 0, data_sigma_err = 0, data_chi2_ndf = 0;
        double mc_mu = 0, mc_mu_err = 0, mc_sigma = 0, mc_sigma_err = 0;
        double jms = 0, jmr = 0;
        // error propagation for the data/MC ratios (mu and sigma uncorrelated)
        double jmsErr() const {
            if (data_mu == 0 || mc_mu == 0) return 0;
            return std::fabs(jms) * std::sqrt(std::pow(data_mu_err / data_mu, 2) +
                                              std::pow(mc_mu_err   / mc_mu,   2));
        }
        double jmrErr() const {
            if (data_sigma == 0 || mc_sigma == 0) return 0;
            return std::fabs(jmr) * std::sqrt(std::pow(data_sigma_err / data_sigma, 2) +
                                              std::pow(mc_sigma_err   / mc_sigma,   2));
        }
    };

    bool ParseParams(const TString& path, FuncResult& r) {
        std::ifstream in(path.Data());
        if (!in.is_open()) return false;
        std::string line;
        while (std::getline(in, line)) {
            if (line.empty() || line[0] == '#') continue;
            std::istringstream ss(line);
            std::string key; ss >> key;
            if      (key == "function")        ss >> r.func;
            else if (key == "is_optimal")      { int v=0; ss >> v; r.is_optimal = (v != 0); }
            else if (key == "bkg_chi2_ndf")    ss >> r.bkg_chi2_ndf;
            else if (key == "spurious_signal") ss >> r.spurious;
            else if (key == "data_mu")         ss >> r.data_mu;
            else if (key == "data_mu_err")     ss >> r.data_mu_err;
            else if (key == "data_sigma")      ss >> r.data_sigma;
            else if (key == "data_sigma_err")  ss >> r.data_sigma_err;
            else if (key == "data_chi2_ndf")   ss >> r.data_chi2_ndf;
            else if (key == "mc_mu")           ss >> r.mc_mu;
            else if (key == "mc_mu_err")       ss >> r.mc_mu_err;
            else if (key == "mc_sigma")        ss >> r.mc_sigma;
            else if (key == "mc_sigma_err")    ss >> r.mc_sigma_err;
            else if (key == "JMS")             ss >> r.jms;
            else if (key == "JMR")             ss >> r.jmr;
        }
        return !r.func.empty();
    }

    std::vector<TString> ListDir(const TString& dir, bool wantDirs) {
        std::vector<TString> out;
        void* d = gSystem->OpenDirectory(dir);
        if (!d) return out;
        const char* e;
        while ((e = gSystem->GetDirEntry(d))) {
            TString name(e);
            if (name == "." || name == "..") continue;
            TString full = dir + "/" + name;
            bool isDir = (gSystem->OpenDirectory(full) != nullptr);
            if (isDir == wantDirs) out.push_back(name);
        }
        gSystem->FreeDirectory(d);
        std::sort(out.begin(), out.end());
        return out;
    }

    TString FindPlotsBase(const TString& topology, const TString& eraSub) {
        const char* cands[] = {"", "../", "EasyBjets/", "../../"};
        for (const char* c : cands) {
            TString p = TString(c) + topology + "/" + eraSub;
            if (!gSystem->AccessPathName(p)) return p;
        }
        return TString("");
    }

    // Spread of a quantity across the functions, relative to the optimal one:
    //   env   max |x - opt| / |opt|
    //   half  (max - min) / (2 |opt|)
    //   rms   sqrt(mean (x - opt)^2) / |opt|
    //   gate  env over the functions with data chi2/ndf <= gate
    struct BkgStat {
        double opt = 0, lo = 0, hi = 0, span = 0;
        double env_pct = 0;   
        double half_pct = 0;  
        double rms_pct = 0;   
        double gate_pct = 0;  
        int    n = 0, n_gate = 0;
    };

    BkgStat ComputeStat(const std::vector<double>& v, const std::vector<double>& chi2,
                     double opt, double gate) {
        BkgStat s; s.opt = opt; s.n = (int)v.size();
        if (v.empty()) return s;
        s.lo = s.hi = v[0];
        double maxdev = 0, sumsq = 0;
        std::vector<double> kept;
        for (size_t i = 0; i < v.size(); ++i) {
            double x = v[i];
            s.lo = std::min(s.lo, x); s.hi = std::max(s.hi, x);
            maxdev = std::max(maxdev, std::fabs(x - opt));
            sumsq += (x - opt) * (x - opt);
            if (chi2[i] <= gate) kept.push_back(x);
        }
        s.span = s.hi - s.lo;
        double a = std::fabs(opt);
        if (a > 0) {
            s.env_pct  = 100.0 * maxdev / a;
            s.half_pct = 100.0 * s.span / (2.0 * a);
            s.rms_pct  = 100.0 * std::sqrt(sumsq / v.size()) / a;
        }
        s.n_gate = (int)kept.size();
        double gdev = 0;
        for (double x : kept) gdev = std::max(gdev, std::fabs(x - opt));
        if (a > 0 && !kept.empty()) s.gate_pct = 100.0 * gdev / a;
        return s;
    }

    // ---- Plotting --------------------------------------------------
    using Getter = std::function<double(const FuncResult&)>;

    struct Series {
        std::string label;
        int    color, marker;
        double dx;          // x offset so overlaid series do not sit on top of each other
        bool   hiOpt;       // hightlight the optimal function in red
        Getter val, err;
    };
    struct PanelSpec {
        std::string header;   
        std::string ytitle;   
        std::vector<Series> series;
        bool legend = false;
    };

    void DrawPanel(const std::vector<FuncResult>& fr, const PanelSpec& P) {
        const int n = (int)fr.size();
        gPad->SetGridy();
        gPad->SetLeftMargin(0.18); gPad->SetBottomMargin(0.14); gPad->SetTopMargin(0.08);

        const FuncResult* opt = &fr[0];
        for (const auto& r : fr) if (r.is_optimal) opt = &r;

        double lo = 1e30, hi = -1e30;
        for (const auto& s : P.series) for (int i = 0; i < n; ++i) {
            double y = s.val(fr[i]), e = s.err(fr[i]);
            lo = std::min(lo, y - e); hi = std::max(hi, y + e);
        }
        if (hi <= lo) { hi = lo + 1; }
        double plo = lo - 0.22 * (hi - lo), phi = hi + 0.22 * (hi - lo);

        TH1F* frame = new TH1F(Form("frame_%p_%s", (void*)&fr, P.ytitle.c_str()),
                               "", n, 0.5, n + 0.5);
        frame->SetStats(0);
        frame->SetTitle(Form("%s;;%s", P.header.c_str(), P.ytitle.c_str()));
        for (int i = 0; i < n; ++i)
            frame->GetXaxis()->SetBinLabel(i + 1, fr[i].func.c_str());
        frame->GetXaxis()->SetLabelSize(0.06);
        frame->GetYaxis()->SetTitleSize(0.05);
        frame->GetYaxis()->SetTitleOffset(1.7);
        frame->GetYaxis()->SetLabelSize(0.042);
        frame->SetMinimum(plo); frame->SetMaximum(phi);
        frame->Draw("AXIS");

        for (const auto& s : P.series) if (s.hiOpt) {
            double ov = s.val(*opt);
            TLine* l = new TLine(0.5, ov, n + 0.5, ov);
            l->SetLineColor(kRed); l->SetLineStyle(2); l->SetLineWidth(2); l->Draw("same");
        }

        TLegend* leg = P.legend ? new TLegend(0.56, 0.78, 0.93, 0.92) : nullptr;
        for (const auto& s : P.series) {
            TGraphErrors* g = new TGraphErrors(n);
            for (int i = 0; i < n; ++i) {
                g->SetPoint(i, i + 1 + s.dx, s.val(fr[i]));
                g->SetPointError(i, 0, s.err(fr[i]));
            }
            g->SetMarkerStyle(s.marker); g->SetMarkerColor(s.color);
            g->SetLineColor(s.color); g->SetMarkerSize(1.3); g->SetLineWidth(2);
            g->Draw("P same");
            if (leg) leg->AddEntry(g, s.label.c_str(), "lp");
            if (s.hiOpt) for (int i = 0; i < n; ++i) if (fr[i].is_optimal) {
                TGraphErrors* go = new TGraphErrors(1);
                go->SetPoint(0, i + 1 + s.dx, s.val(fr[i]));
                go->SetPointError(0, 0, s.err(fr[i]));
                go->SetMarkerStyle(s.marker); go->SetMarkerSize(1.6); go->SetMarkerColor(kRed);
                go->SetLineColor(kRed); go->SetLineWidth(2); go->Draw("P same");
            }
        }
        if (leg) { leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.05); leg->Draw(); }
    }

    void MakeCellPlot(const std::string& title, const std::string& subtitle,
                      const std::string& outPdf, const std::vector<FuncResult>& fr,
                      const PanelSpec& left, const PanelSpec& right) {
        TCanvas c("c_stab", "stab", 1000, 480);
        TPad* titlePad = new TPad("titlePad", "", 0, 0.90, 1, 1.0);
        TPad* bodyPad  = new TPad("bodyPad",  "", 0, 0.0,  1, 0.90);
        titlePad->Draw(); bodyPad->Draw();
        bodyPad->Divide(2, 1);

        titlePad->cd();
        TLatex tt; tt.SetNDC(); tt.SetTextAlign(22);
        tt.SetTextSize(0.40); tt.DrawLatex(0.5, 0.58, title.c_str());
        tt.SetTextSize(0.30); tt.DrawLatex(0.5, 0.18, subtitle.c_str());

        bodyPad->cd(1); DrawPanel(fr, left);
        bodyPad->cd(2); DrawPanel(fr, right);

        c.SaveAs(outPdf.c_str());
        std::string png = outPdf;
        if (png.size() > 4) png.replace(png.size() - 4, 4, ".png");
        c.SaveAs(png.c_str());
    }

}  // namespace

void MassFit_BackgroundStability(const char* era = "run2",
                                 const char* channel = "DB",
                                 const char* variantSel = "ALL",
                                 double chi2Gate = 2.0) {
    gErrorIgnoreLevel = kWarning;

    TString ch(channel); ch.ToUpper();
    TString topology = (ch == "MJB") ? "ZbbJets" : "ZbbGamma";
    TString eraSub   = TString(era).Contains("3") ? "Run3Plots" : "Run2Plots";
    TString fitTag   = "Mass_Fit_" + ch;

    TString base = FindPlotsBase(topology, eraSub);
    if (base.IsNull()) {
        printf("ERROR: cannot locate %s/%s from here (%s).\n",
               topology.Data(), eraSub.Data(), gSystem->pwd());
        printf("Run from EasyBjets/CrossChecks/ (or EasyBjets/).\n");
        return;
    }

    TString outDir = base + "/" + ch + "_BackgroundStability";
    gSystem->Exec("mkdir -p " + outDir);

    std::ofstream csv((outDir + "/background_stability_summary.csv").Data());
    std::ofstream tab((outDir + "/background_systematics_table.txt").Data());
    csv << "variant,wp,ptbin,n_func,optimal,"
        << "mu_opt,mu_min,mu_max,mu_span,mu_env_pct,mu_half_pct,mu_rms_pct,mu_gate_pct,mu_n_gate,"
        << "sigma_opt,sigma_min,sigma_max,sigma_span,sigma_env_pct,sigma_half_pct,sigma_rms_pct,sigma_gate_pct,sigma_n_gate,"
        << "bJMDP_opt,bJMDR_opt\n";

    tab << "Background-modeling stability - " << topology << " / " << eraSub << " / " << ch << "\n";
    tab << "Systematic = spread of the DATA peak (mu) / width (sigma) across background functions.\n";
    tab << "  env  = max|func - optimal| / optimal           (envelope, RECOMMENDED)\n";
    tab << "  half = (max - min) / (2*optimal)               (half full-span)\n";
    tab << "  rms  = sqrt(mean((func-optimal)^2)) / optimal  (rms about optimal)\n";
    tab << "  gate = env over functions with data chi2/ndf <= " << chi2Gate << "\n";
    tab << "b-JMDP = mu_data/mu_MC (scale),  b-JMDR = sigma_data/sigma_MC (resolution).\n";
    tab << "Since mc is function-independent, b-JMDP%% = peak-mu%%, b-JMDR%% = width-sigma%%.\n";
    tab << std::string(118, '=') << "\n";

    double worst_mu_env = 0, worst_sig_env = 0, worst_mu_gate = 0, worst_sig_gate = 0;
    int nCells = 0;

    for (const TString& variant : ListDir(base, true)) {
        if (TString(variantSel) != "ALL" && variant != variantSel) continue;
        TString fitRoot = base + "/" + variant + "/" + fitTag;
        if (gSystem->AccessPathName(fitRoot)) continue;

        for (const TString& wp : ListDir(fitRoot, true)) {
            TString wpDir = fitRoot + "/" + wp;
            for (const TString& ptbin : ListDir(wpDir, true)) {
                TString cellDir = wpDir + "/" + ptbin;
                std::vector<FuncResult> fr;
                for (const TString& func : ListDir(cellDir, true)) {
                    FuncResult r;
                    if (ParseParams(cellDir + "/" + func + "/signal_fit_params.txt", r))
                        fr.push_back(r);
                }
                if (fr.size() < 2) continue;

                const FuncResult* opt = &fr[0];
                for (const auto& r : fr) if (r.is_optimal) opt = &r;

                std::vector<double> mus, sigs, chi2;
                for (const auto& r : fr) {
                    mus.push_back(r.data_mu); sigs.push_back(r.data_sigma);
                    chi2.push_back(r.data_chi2_ndf);
                }
                BkgStat sMu  = ComputeStat(mus,  chi2, opt->data_mu,    chi2Gate);
                BkgStat sSig = ComputeStat(sigs, chi2, opt->data_sigma, chi2Gate);

                worst_mu_env   = std::max(worst_mu_env,   sMu.env_pct);
                worst_sig_env  = std::max(worst_sig_env,  sSig.env_pct);
                worst_mu_gate  = std::max(worst_mu_gate,  sMu.gate_pct);
                worst_sig_gate = std::max(worst_sig_gate, sSig.gate_pct);
                ++nCells;

                csv << variant << "," << wp << "," << ptbin << "," << fr.size() << ","
                    << opt->func << ","
                    << opt->data_mu << "," << sMu.lo << "," << sMu.hi << "," << sMu.span << ","
                    << sMu.env_pct << "," << sMu.half_pct << "," << sMu.rms_pct << ","
                    << sMu.gate_pct << "," << sMu.n_gate << ","
                    << opt->data_sigma << "," << sSig.lo << "," << sSig.hi << "," << sSig.span << ","
                    << sSig.env_pct << "," << sSig.half_pct << "," << sSig.rms_pct << ","
                    << sSig.gate_pct << "," << sSig.n_gate << ","
                    << opt->jms << "," << opt->jmr << "\n";

                // ---- text table ----
                tab << "\n" << variant << "  WP " << wp << "  " << ptbin
                    << "   (optimal: " << opt->func << ")\n";
                tab << "  " << std::left << std::setw(12) << "function" << std::right
                    << std::setw(11) << "mu_data" << std::setw(13) << "sigma_data"
                    << std::setw(11) << "b-JMDP" << std::setw(11) << "b-JMDR"
                    << std::setw(10) << "dataChi2" << "  opt\n";
                for (const auto& r : fr)
                    tab << "  " << std::left << std::setw(12) << r.func << std::right
                        << std::setw(11) << Form("%.3f", r.data_mu)
                        << std::setw(13) << Form("%.3f", r.data_sigma)
                        << std::setw(11) << Form("%.4f", r.jms)
                        << std::setw(11) << Form("%.4f", r.jmr)
                        << std::setw(10) << Form("%.2f", r.data_chi2_ndf)
                        << (r.is_optimal ? "   <==" : "") << "\n";
                tab << "  MC (function-independent):  mu = " << Form("%.3f", opt->mc_mu)
                    << " GeV,  sigma = " << Form("%.3f", opt->mc_sigma) << " GeV\n";
                tab << "  ---- systematics ----            env       half      rms     gate(chi2<="
                    << chi2Gate << ", n=" << sMu.n_gate << ")\n";
                tab << "  peak  mu / b-JMDP [%]" << std::right
                    << std::setw(9)  << Form("%.3f", sMu.env_pct)
                    << std::setw(10) << Form("%.3f", sMu.half_pct)
                    << std::setw(9)  << Form("%.3f", sMu.rms_pct)
                    << std::setw(11) << Form("%.3f", sMu.gate_pct) << "\n";
                tab << "  width sigma / b-JMDR [%]" << std::right
                    << std::setw(6)  << Form("%.3f", sSig.env_pct)
                    << std::setw(10) << Form("%.3f", sSig.half_pct)
                    << std::setw(9)  << Form("%.3f", sSig.rms_pct)
                    << std::setw(11) << Form("%.3f", sSig.gate_pct) << "\n";

                // ---- getters ----
                Getter gDmu  = [](const FuncResult& r){ return r.data_mu; };
                Getter gDmuE = [](const FuncResult& r){ return r.data_mu_err; };
                Getter gDsg  = [](const FuncResult& r){ return r.data_sigma; };
                Getter gDsgE = [](const FuncResult& r){ return r.data_sigma_err; };
                Getter gMmu  = [](const FuncResult& r){ return r.mc_mu; };
                Getter gMmuE = [](const FuncResult& r){ return r.mc_mu_err; };
                Getter gMsg  = [](const FuncResult& r){ return r.mc_sigma; };
                Getter gMsgE = [](const FuncResult& r){ return r.mc_sigma_err; };
                Getter gJms  = [](const FuncResult& r){ return r.jms; };
                Getter gJmsE = [](const FuncResult& r){ return r.jmsErr(); };
                Getter gJmr  = [](const FuncResult& r){ return r.jmr; };
                Getter gJmrE = [](const FuncResult& r){ return r.jmrErr(); };

                std::string cellTitle = (variant + "  " + wp + "  " + ptbin).Data();
                std::string stem = (outDir + "/" + variant + "_" + wp + "_" + ptbin).Data();

                // (1) data
                MakeCellPlot(cellTitle, "Data (extracted signal) - red = optimal function",
                    stem + "_data_stability.pdf", fr,
                    {"Peak  #mu_{data}", "#mu_{data} [GeV]",  {{"Data",kBlack,20,0.0,true,gDmu,gDmuE}}, false},
                    {"Width  #sigma_{data}","#sigma_{data} [GeV]",{{"Data",kBlack,20,0.0,true,gDsg,gDsgE}}, false});

                // (2) MC (flat: background-function independent)
                MakeCellPlot(cellTitle, "MC (signal template, fit directly) - flat: bkg-function independent",
                    stem + "_mc_stability.pdf", fr,
                    {"Peak  #mu_{MC}",  "#mu_{MC} [GeV]",  {{"MC",kBlue+1,21,0.0,false,gMmu,gMmuE}}, false},
                    {"Width  #sigma_{MC}","#sigma_{MC} [GeV]",{{"MC",kBlue+1,21,0.0,false,gMsg,gMsgE}}, false});

                // (3) data + MC overlaid (x-offset)
                MakeCellPlot(cellTitle, "Data vs MC",
                    stem + "_datamc_stability.pdf", fr,
                    {"Peak  #mu","#mu [GeV]",
                        {{"Data",kBlack,20,-0.12,true,gDmu,gDmuE},{"MC",kBlue+1,21,+0.12,false,gMmu,gMmuE}}, true},
                    {"Width  #sigma","#sigma [GeV]",
                        {{"Data",kBlack,20,-0.12,true,gDsg,gDsgE},{"MC",kBlue+1,21,+0.12,false,gMsg,gMsgE}}, true});

                // (4) observables: b-JMDP / b-JMDR vs function
                MakeCellPlot(cellTitle, "Calibration observables - red = optimal function",
                    stem + "_bJMD_stability.pdf", fr,
                    {"b-JMDP = #mu_{data}/#mu_{MC}", "b-JMDP", {{"b-JMDP",kGreen+2,20,0.0,true,gJms,gJmsE}}, false},
                    {"b-JMDR = #sigma_{data}/#sigma_{MC}","b-JMDR",{{"b-JMDR",kGreen+2,20,0.0,true,gJmr,gJmrE}}, false});
            }
        }
    }

    tab << "\n" << std::string(118, '=') << "\n";
    tab << Form("OVERALL (%d cells):\n", nCells);
    tab << Form("  worst peak-mu / b-JMDP  systematic: env = %.3f%%,  gated = %.3f%%\n", worst_mu_env, worst_mu_gate);
    tab << Form("  worst width-s / b-JMDR  systematic: env = %.3f%%,  gated = %.3f%%\n", worst_sig_env, worst_sig_gate);
    csv.close(); tab.close();

    printf("\nBackground-stability evaluation done (%d cells).\n", nCells);
    printf("  worst peak  / b-JMDP systematic: env=%.3f%%  gated=%.3f%%\n", worst_mu_env, worst_mu_gate);
    printf("  worst width / b-JMDR systematic: env=%.3f%%  gated=%.3f%%\n", worst_sig_env, worst_sig_gate);
    printf("  table: %s/background_systematics_table.txt\n", outDir.Data());
    printf("  csv:   %s/background_stability_summary.csv\n", outDir.Data());
    printf("  plots: %s/<cell>_{data,mc,datamc,bJMD}_stability.{pdf,png}\n", outDir.Data());
}
