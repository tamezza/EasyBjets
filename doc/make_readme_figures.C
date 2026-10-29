/**
 * make_readme_figures.C - toy sketches for the README
 *
 * Draws illustration plots from random numbers (fixed random seed), NOT
 * from data or MC: they only show what the analysis measures and how to read
 * its plots. Every plot is labelled "Illustration (toy, not data)".
 *
 *   figures/sketch_mass_stack.png        data vs signal + background MC, sideband SF
 *   figures/sketch_bkg_composition.png   flavour composition of the background jets
 *   figures/sketch_mass_fit.png          mass peak on a background, sidebands, mu and sigma
 *   figures/sketch_mass_extracted.png    extracted Z peak in data and in signal MC
 *   figures/sketch_mass_summary.png      b-JMDP and b-JMDR: single and double ratios vs pT
 *   figures/sketch_response_stack.png    balance R: data vs signal + background MC
 *   figures/sketch_response.png          balance R with the core fit and the truth width
 *   figures/sketch_response_summary.png  b-JES and b-JER: single and double ratios vs pT
 *   figures/sketch_bootstrap.png         correction vs nominal fit results over replicas
 *   figures/sketch_double_ratio.png      double ratio vs pT, the four error treatments
 *
 * Usage (from EasyBjets/doc):
 *   root make_readme_figures.C
 */

#include "../common/Utils/AtlasStyle.C"

#include <TCanvas.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TArrow.h>
#include <TBox.h>
#include <THStack.h>
#include <TPad.h>
#include <TRandom3.h>
#include <TSystem.h>
#include <TStyle.h>

#include <cmath>

namespace Sketch {

void Label(double x, double y, const char* extra = "", double size = 0.045) {
    TLatex t;
    t.SetNDC();
    t.SetTextFont(42);
    t.SetTextSize(size);
    t.DrawLatex(x, y, "#font[72]{Illustration}  (toy, not data)");
    if (extra[0]) {
        t.SetTextSize(0.85 * size);
        t.DrawLatex(x, y - 1.25 * size, extra);
    }
}

// Upper pad (distribution) and lower pad (data / MC) of a stack plot.
void SplitPads(TPad*& top, TPad*& bot) {
    top = new TPad("top", "", 0, 0.30, 1, 1);
    top->SetBottomMargin(0.02);
    top->Draw();
    bot = new TPad("bot", "", 0, 0, 1, 0.30);
    bot->SetTopMargin(0.03);
    bot->SetBottomMargin(0.38);
    bot->Draw();
}

// Data / MC ratio in the lower pad, with a line at 1.
void RatioPanel(TPad* bot, TH1F* hData, TH1F* hMC, const char* xtitle, double lo, double hi) {
    bot->cd();
    TH1F* r = (TH1F*)hData->Clone(Form("%s_ratio", hData->GetName()));
    r->Divide(hMC);
    r->SetTitle(Form(";%s;Data / MC", xtitle));
    r->SetMinimum(lo);
    r->SetMaximum(hi);
    r->GetYaxis()->SetNdivisions(504);
    r->GetXaxis()->SetTitleSize(0.14);
    r->GetXaxis()->SetLabelSize(0.12);
    r->GetYaxis()->SetTitleSize(0.12);
    r->GetYaxis()->SetLabelSize(0.11);
    r->GetYaxis()->SetTitleOffset(0.55);
    r->GetXaxis()->SetTitleOffset(1.2);
    r->Draw("E");
    TLine* one = new TLine(r->GetXaxis()->GetXmin(), 1, r->GetXaxis()->GetXmax(), 1);
    one->SetLineStyle(2);
    one->Draw();
}

// One summary panel pair (single ratio on top, double ratio below) as in the
// combined summaries: nominal and one correction vs pT, DB then MJB.
void SummaryPanels(TPad* pad, const char* qty, const double* srNom, const double* srCor, const double* err,
                   double srLo, double srHi, double drLo, double drHi) {
    const int n = 5;
    const double pt[n] = {250, 375, 525, 750, 1050};
    pad->cd();
    TPad* top = new TPad(Form("t_%s", qty), "", 0, 0.42, 1, 1);
    top->SetBottomMargin(0.02);
    top->SetLeftMargin(0.17);
    top->Draw();
    TPad* bot = new TPad(Form("b_%s", qty), "", 0, 0, 1, 0.42);
    bot->SetTopMargin(0.03);
    bot->SetBottomMargin(0.30);
    bot->SetLeftMargin(0.17);
    bot->Draw();

    top->cd();
    TH1F* f1 = new TH1F(Form("f1_%s", qty), Form(";;%s  data / MC", qty), 10, 200, 1200);
    f1->SetMinimum(srLo);
    f1->SetMaximum(srHi);
    f1->GetXaxis()->SetLabelSize(0);
    f1->GetYaxis()->SetTitleSize(0.075);
    f1->GetYaxis()->SetLabelSize(0.065);
    f1->GetYaxis()->SetTitleOffset(1.0);
    f1->Draw();
    TGraphErrors* gN = new TGraphErrors(n);
    TGraphErrors* gC = new TGraphErrors(n);
    for (int i = 0; i < n; ++i) {
        gN->SetPoint(i, pt[i] - 12, srNom[i]);
        gN->SetPointError(i, 0, err[i]);
        gC->SetPoint(i, pt[i] + 12, srCor[i]);
        gC->SetPointError(i, 0, err[i]);
    }
    gN->SetMarkerStyle(20);
    gC->SetMarkerStyle(21);
    gC->SetMarkerColor(kRed + 1);
    gC->SetLineColor(kRed + 1);
    gN->Draw("P SAME");
    gC->Draw("P SAME");
    TLine* l1 = new TLine(200, 1, 1200, 1);
    l1->SetLineStyle(2);
    l1->Draw();
    TLine* s1 = new TLine(450, srLo, 450, srHi);
    s1->SetLineStyle(3);
    s1->SetLineColor(kGray + 2);
    s1->Draw();
    TLegend* leg = new TLegend(0.55, 0.62, 0.92, 0.78);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.06);
    leg->AddEntry(gN, "nominal", "pe");
    leg->AddEntry(gC, "correction", "pe");
    leg->Draw();
    TLatex t;
    t.SetNDC();
    t.SetTextSize(0.06);
    t.SetTextColor(kGray + 2);
    t.DrawLatex(0.22, 0.06, "DB");
    t.DrawLatex(0.50, 0.06, "MJB");
    Label(0.21, 0.86, "", 0.065);

    bot->cd();
    TH1F* f2 = new TH1F(Form("f2_%s", qty), ";Large-R jet p_{T} [GeV];correction / nominal", 10, 200, 1200);
    f2->SetMinimum(drLo);
    f2->SetMaximum(drHi);
    f2->GetXaxis()->SetTitleSize(0.10);
    f2->GetXaxis()->SetLabelSize(0.09);
    f2->GetYaxis()->SetTitleSize(0.085);
    f2->GetYaxis()->SetLabelSize(0.085);
    f2->GetYaxis()->SetTitleOffset(0.85);
    f2->GetYaxis()->SetNdivisions(504);
    f2->Draw();
    TGraphErrors* gD = new TGraphErrors(n);
    for (int i = 0; i < n; ++i) {
        gD->SetPoint(i, pt[i], srCor[i] / srNom[i]);
        gD->SetPointError(i, 0, 0.5 * err[i]);
    }
    gD->SetMarkerStyle(21);
    gD->SetMarkerColor(kRed + 1);
    gD->SetLineColor(kRed + 1);
    gD->Draw("P SAME");
    TLine* l2 = new TLine(200, 1, 1200, 1);
    l2->SetLineStyle(2);
    l2->Draw();
    TLine* s2 = new TLine(450, drLo, 450, drHi);
    s2->SetLineStyle(3);
    s2->SetLineColor(kGray + 2);
    s2->Draw();
}

// ---------------------------------------------------------------------------
// Mass stack: data vs Z(bb) signal + main background MC, the background
// scaled to data in the sidebands (SF), with the data/MC ratio.
// ---------------------------------------------------------------------------
void MassStack(TRandom3& rng) {
    TCanvas c("c_mstack", "", 800, 700);
    TPad *top, *bot;
    SplitPads(top, bot);
    top->cd();
    TH1F* hB = new TH1F("hB", ";;Events / 5 GeV", 30, 50, 200);
    TH1F* hS = new TH1F("hS", "", 30, 50, 200);
    TH1F* hD = new TH1F("hD", "", 30, 50, 200);
    for (int b = 1; b <= 30; ++b) {
        const double x = hB->GetBinCenter(b);
        const double bkg = 900 * std::exp(-0.022 * (x - 50));
        const double sig = 250 * std::exp(-0.5 * std::pow((x - 91) / 10, 2));
        hB->SetBinContent(b, bkg);
        hS->SetBinContent(b, sig);
        hD->SetBinContent(b, rng.Poisson(bkg * (1.0 + 0.0006 * (x - 100)) + 1.05 * sig));
        hD->SetBinError(b, std::sqrt(hD->GetBinContent(b)));
    }
    hB->SetFillColor(kAzure - 9);
    hB->SetLineColor(kAzure + 2);
    hS->SetFillColor(kOrange - 3);
    hS->SetLineColor(kOrange + 7);
    THStack* st = new THStack("st", ";;Events / 5 GeV");
    st->Add(hB);
    st->Add(hS);
    const double yTop = hD->GetMaximum() * 1.55;
    st->SetMaximum(yTop);
    st->Draw("HIST");
    st->GetXaxis()->SetLabelSize(0);
    st->GetYaxis()->SetTitleSize(0.055);
    st->GetYaxis()->SetLabelSize(0.05);
    for (double lo : {50.0, 115.0}) {
        TBox* bx = new TBox(lo, 0, lo == 50 ? 70.0 : 200.0, yTop);
        bx->SetFillColorAlpha(kGray, 0.30);
        bx->Draw();
    }
    st->Draw("HIST SAME");
    hD->SetMarkerStyle(20);
    hD->Draw("E SAME");
    gPad->RedrawAxis();
    TLegend* leg = new TLegend(0.50, 0.55, 0.90, 0.76);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.040);
    leg->AddEntry(hD, "Data", "pe");
    leg->AddEntry(hS, "Z(#rightarrow bb) signal MC", "f");
    leg->AddEntry(hB, "Main background MC #times SF(m)", "f");
    leg->Draw();
    TLatex t;
    t.SetTextSize(0.04);
    t.SetTextColor(kGray + 2);
    t.DrawLatex(51.5, yTop * 0.71, "side-");
    t.DrawLatex(51.5, yTop * 0.66, "band");
    t.DrawLatex(150, yTop * 0.35, "sideband");
    Label(0.18, 0.86, "SF(m) = p_{0} + p_{1} m fitted to data / MC in the sidebands", 0.048);

    TH1F* hMC = (TH1F*)hB->Clone("hMC");
    hMC->Add(hS);
    RatioPanel(bot, hD, hMC, "Large-R jet mass m_{J} [GeV]", 0.8, 1.2);
    c.SaveAs("figures/sketch_mass_stack.png");
}

// ---------------------------------------------------------------------------
// Extracted signal: data minus the fitted background, and the signal MC,
// each with its Gaussian peak fit (mu, sigma).
// ---------------------------------------------------------------------------
void MassExtracted(TRandom3& rng) {
    TCanvas c("c_mext", "", 1200, 550);
    c.Divide(2, 1);
    struct Side { const char* name; const char* title; double mu, sg, amp, noise; int col; };
    const Side sides[2] = {{"data", "Data #minus background", 92.0, 11.0, 260, 1.0, kBlack},
                           {"mc", "Z(#rightarrow bb) signal MC", 90.5, 10.0, 240, 0.25, kOrange + 7}};
    for (int k = 0; k < 2; ++k) {
        c.cd(k + 1);
        gPad->SetLeftMargin(0.16);
        const Side& s = sides[k];
        TH1F* h = new TH1F(Form("hx_%s", s.name), ";Large-R jet mass m_{J} [GeV];Events / 5 GeV",
                           30, 50, 200);
        for (int b = 1; b <= 30; ++b) {
            const double x = h->GetBinCenter(b);
            const double v = s.amp * std::exp(-0.5 * std::pow((x - s.mu) / s.sg, 2));
            const double e = std::sqrt(v + 400 * std::exp(-0.022 * (x - 50))) * s.noise;
            h->SetBinContent(b, v + rng.Gaus(0, e));
            h->SetBinError(b, e);
        }
        TF1* g = new TF1(Form("g_%s", s.name), "gaus", s.mu - 1.5 * s.sg, s.mu + 1.5 * s.sg);
        g->SetParameters(s.amp, s.mu, s.sg);
        h->Fit(g, "RQ0N");
        h->SetMaximum(s.amp * 1.6);
        h->SetMinimum(-50);
        h->SetMarkerStyle(20);
        h->SetMarkerColor(s.col);
        h->SetLineColor(s.col);
        h->Draw("E");
        g->SetLineColor(kRed + 1);
        g->SetLineWidth(3);
        g->Draw("SAME");
        TLine* z = new TLine(50, 0, 200, 0);
        z->SetLineStyle(2);
        z->Draw();
        TLatex t;
        t.SetNDC();
        t.SetTextSize(0.045);
        t.DrawLatex(0.58, 0.70, s.title);
        const char* sub = k == 0 ? "data" : "MC";
        t.DrawLatex(0.58, 0.63, Form("#mu_{%s} = %.1f GeV", sub, g->GetParameter(1)));
        t.DrawLatex(0.58, 0.56, Form("#sigma_{%s} = %.1f GeV", sub, g->GetParameter(2)));
        Label(0.20, 0.86, "", 0.045);
    }
    c.cd();
    c.SaveAs("figures/sketch_mass_extracted.png");
}

// ---------------------------------------------------------------------------
// Summary of b-JMDP and b-JMDR (left / right).
// ---------------------------------------------------------------------------
void MassSummary() {
    TCanvas c("c_msum", "", 1300, 650);
    TPad* pl = new TPad("pl", "", 0, 0, 0.5, 1);
    TPad* pr = new TPad("pr", "", 0.5, 0, 1, 1);
    pl->Draw();
    pr->Draw();
    const double pNom[5] = {1.012, 1.008, 1.010, 1.006, 1.004};
    const double pCor[5] = {1.006, 1.004, 1.007, 1.003, 1.002};
    const double pErr[5] = {0.008, 0.010, 0.009, 0.012, 0.016};
    const double rNom[5] = {1.06, 1.05, 1.08, 1.04, 1.07};
    const double rCor[5] = {1.03, 1.02, 1.05, 1.03, 1.05};
    const double rErr[5] = {0.04, 0.05, 0.05, 0.06, 0.08};
    SummaryPanels(pl, "b-JMDP", pNom, pCor, pErr, 0.96, 1.06, 0.97, 1.03);
    SummaryPanels(pr, "b-JMDR", rNom, rCor, rErr, 0.80, 1.35, 0.85, 1.15);
    c.SaveAs("figures/sketch_mass_summary.png");
}

// ---------------------------------------------------------------------------
// Response stack: R = pT(J)/pT(ref), data vs signal + background MC.
// ---------------------------------------------------------------------------
void ResponseStack(TRandom3& rng) {
    TCanvas c("c_rstack", "", 800, 700);
    TPad *top, *bot;
    SplitPads(top, bot);
    top->cd();
    TH1F* hB = new TH1F("hRB", "", 40, 0.4, 1.6);
    TH1F* hS = new TH1F("hRS", "", 40, 0.4, 1.6);
    TH1F* hD = new TH1F("hRD", "", 40, 0.4, 1.6);
    for (int i = 0; i < 60000; ++i) {
        hB->Fill(rng.Gaus(0.95, 0.13), 0.5);
        hS->Fill(rng.Gaus(0.97, 0.11), 0.1);
    }
    for (int b = 1; b <= 40; ++b) {
        const double mc = hB->GetBinContent(b) + hS->GetBinContent(b);
        const double shift = 1.0 + 0.25 * (hB->GetBinCenter(b) - 0.96);  // data slightly to the right
        hD->SetBinContent(b, rng.Poisson(std::max(0.0, mc * shift)));
        hD->SetBinError(b, std::sqrt(hD->GetBinContent(b)));
    }
    hB->SetFillColor(kAzure - 9);
    hB->SetLineColor(kAzure + 2);
    hS->SetFillColor(kOrange - 3);
    hS->SetLineColor(kOrange + 7);
    THStack* st = new THStack("rst", ";;Events");
    st->Add(hB);
    st->Add(hS);
    st->SetMaximum(hD->GetMaximum() * 1.5);
    st->Draw("HIST");
    st->GetXaxis()->SetLabelSize(0);
    st->GetYaxis()->SetTitleSize(0.055);
    st->GetYaxis()->SetLabelSize(0.05);
    hD->SetMarkerStyle(20);
    hD->Draw("E SAME");
    TLegend* leg = new TLegend(0.60, 0.58, 0.92, 0.78);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.042);
    leg->AddEntry(hD, "Data", "pe");
    leg->AddEntry(hS, "Z(#rightarrow bb) signal MC", "f");
    leg->AddEntry(hB, "Main background MC", "f");
    leg->Draw();
    Label(0.18, 0.86, "R_{DB} = p_{T}^{J} / p_{T}^{#gamma} (DB),  R_{MJB} = p_{T}^{J} / p_{T}^{recoil} (MJB)",
          0.048);
    TH1F* hMC = (TH1F*)hB->Clone("hRMC");
    hMC->Add(hS);
    RatioPanel(bot, hD, hMC, "Balance R = p_{T}^{J} / p_{T}^{ref}", 0.6, 1.4);
    c.SaveAs("figures/sketch_response_stack.png");
}

// ---------------------------------------------------------------------------
// Summary of b-JES and b-JER (left / right).
// ---------------------------------------------------------------------------
void ResponseSummary() {
    TCanvas c("c_rsum", "", 1300, 650);
    TPad* pl = new TPad("pl", "", 0, 0, 0.5, 1);
    TPad* pr = new TPad("pr", "", 0.5, 0, 1, 1);
    pl->Draw();
    pr->Draw();
    const double jNom[5] = {1.010, 1.006, 1.004, 1.008, 1.003};
    const double jCor[5] = {1.004, 1.003, 1.001, 1.004, 1.001};
    const double jErr[5] = {0.006, 0.007, 0.006, 0.008, 0.011};
    const double rNom[5] = {1.08, 1.06, 1.10, 1.05, 1.09};
    const double rCor[5] = {1.04, 1.03, 1.06, 1.04, 1.06};
    const double rErr[5] = {0.05, 0.06, 0.06, 0.07, 0.09};
    SummaryPanels(pl, "R_{JES}", jNom, jCor, jErr, 0.97, 1.04, 0.98, 1.02);
    SummaryPanels(pr, "R_{JER}", rNom, rCor, rErr, 0.80, 1.40, 0.85, 1.15);
    c.SaveAs("figures/sketch_response_summary.png");
}

// ---------------------------------------------------------------------------
// Flavour composition of the background jets (b / c hadrons in the jet).
// ---------------------------------------------------------------------------
void BkgComposition() {
    TCanvas c("c_comp", "", 800, 600);
    const char* lab[6] = {"BB", "BC", "BL", "CC", "CL", "LL"};
    const double db[6] = {0.42, 0.14, 0.15, 0.07, 0.09, 0.13};
    const double mjb[6] = {0.48, 0.13, 0.13, 0.06, 0.08, 0.12};
    TH1F* hDB = new TH1F("hDB", ";Flavour of the background jet;Fraction of the background", 6, 0, 6);
    TH1F* hMJB = new TH1F("hMJB", "", 6, 0, 6);
    for (int i = 0; i < 6; ++i) {
        hDB->SetBinContent(i + 1, db[i]);
        hMJB->SetBinContent(i + 1, mjb[i]);
        hDB->GetXaxis()->SetBinLabel(i + 1, lab[i]);
    }
    hDB->SetMaximum(0.72);
    hDB->SetMinimum(0);
    hDB->GetXaxis()->SetLabelSize(0.06);
    hDB->SetFillColor(kAzure - 4);
    hDB->SetBarWidth(0.38);
    hDB->SetBarOffset(0.10);
    hMJB->SetFillColor(kOrange - 3);
    hMJB->SetBarWidth(0.38);
    hMJB->SetBarOffset(0.52);
    hDB->Draw("BAR");
    hMJB->Draw("BAR SAME");
    TLegend* leg = new TLegend(0.60, 0.62, 0.92, 0.76);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.04);
    leg->AddEntry(hDB, "DB:  #gamma+jets", "f");
    leg->AddEntry(hMJB, "MJB:  multijets", "f");
    leg->Draw();
    Label(0.20, 0.86, "BB #geq 2 b,  BC 1 b + c,  BL 1 b,  CC #geq 2 c,  CL 1 c,  LL light", 0.042);
    c.SaveAs("figures/sketch_bkg_composition.png");
}

// ---------------------------------------------------------------------------
// 1. Mass fit: falling background + Z peak, sidebands, background-only fit
//    under the peak, signal = data - background.
// ---------------------------------------------------------------------------
void MassFit(TRandom3& rng) {
    TCanvas c("c_mass", "", 800, 600);
    TH1F h("h_mass", ";Large-R jet mass m_{J} [GeV];Events / 5 GeV", 30, 50, 200);
    TF1 bkg("bkg", "[0]*exp(-[1]*(x-50))", 50, 200);
    bkg.SetParameters(900, 0.022);
    TF1 sig("sig", "gaus", 50, 200);
    sig.SetParameters(260, 91, 10);
    for (int b = 1; b <= h.GetNbinsX(); ++b) {
        double x = h.GetBinCenter(b);
        h.SetBinContent(b, rng.Poisson(bkg.Eval(x) + sig.Eval(x)));
        h.SetBinError(b, std::sqrt(h.GetBinContent(b)));
    }
    TF1 sb("sb", "[0]*exp(-[1]*(x-50)) + gaus(2)", 50, 200);
    sb.SetParameters(900, 0.022, 260, 91, 10);
    h.Fit(&sb, "RQ0N");
    TF1 bfit("bfit", "[0]*exp(-[1]*(x-50))", 50, 200);
    bfit.SetParameters(sb.GetParameter(0), sb.GetParameter(1));
    TF1 sfit("sfit", "gaus", 50, 200);
    sfit.SetParameters(sb.GetParameter(2), sb.GetParameter(3), sb.GetParameter(4));

    h.SetMaximum(h.GetMaximum() * 1.45);
    h.SetMinimum(0);
    h.SetMarkerStyle(20);
    h.Draw("E");

    // sidebands
    const double ymax = h.GetMaximum();
    TBox sb1(50, 0, 70, ymax), sb2(115, 0, 200, ymax);
    for (TBox* bx : {&sb1, &sb2}) {
        bx->SetFillColorAlpha(kGray, 0.35);
        bx->Draw();
    }
    h.Draw("E SAME");

    bfit.SetLineColor(kBlue + 1);
    bfit.SetLineStyle(2);
    bfit.SetLineWidth(3);
    bfit.Draw("SAME");
    sfit.SetLineColor(kOrange + 7);
    sfit.SetLineStyle(2);
    sfit.SetLineWidth(3);
    sfit.Draw("SAME");
    sb.SetLineColor(kRed + 1);
    sb.SetLineWidth(3);
    sb.Draw("SAME");

    // mu and sigma of the peak
    const double mu = sfit.GetParameter(1), sg = sfit.GetParameter(2);
    const double yPeak = bfit.Eval(mu) + sfit.Eval(mu);
    TLine lmu(mu, 0, mu, yPeak);
    lmu.SetLineStyle(3);
    lmu.SetLineWidth(2);
    lmu.Draw();
    TArrow aw(mu - sg, 0.5 * sfit.Eval(mu) + bfit.Eval(mu), mu + sg, 0.5 * sfit.Eval(mu) + bfit.Eval(mu),
              0.015, "<|>");
    aw.SetLineWidth(2);
    aw.Draw();
    TLatex t;
    t.SetTextSize(0.04);
    t.DrawLatex(mu + 2, yPeak * 1.03, "#mu");
    t.DrawLatex(mu + sg + 2, 0.5 * sfit.Eval(mu) + bfit.Eval(mu), "#pm#sigma");
    t.SetTextColor(kGray + 2);
    t.DrawLatex(52, ymax * 0.05, "sideband");
    t.DrawLatex(150, ymax * 0.30, "sideband");

    TLegend leg(0.56, 0.58, 0.92, 0.80);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    leg.SetTextSize(0.035);
    leg.AddEntry(&h, "Data", "pe");
    leg.AddEntry(&bfit, "Background (sideband fit)", "l");
    leg.AddEntry(&sfit, "Z #rightarrow bb peak", "l");
    leg.AddEntry(&sb, "Signal + background", "l");
    leg.Draw();
    Label(0.20, 0.87, "b-JMDP = #mu_{data} / #mu_{MC},   b-JMDR = #sigma_{data} / #sigma_{MC}");
    c.SaveAs("figures/sketch_mass_fit.png");
}

// ---------------------------------------------------------------------------
// 2. Response: R = pT(J)/pT(ref) with tails, Gaussian core fit, and the
//    narrower truth balance pT(truth)/pT(ref).
// ---------------------------------------------------------------------------
void Response(TRandom3& rng) {
    TCanvas c("c_resp", "", 800, 600);
    TH1F hR("hR", ";Balance R = p_{T}^{J} / p_{T}^{ref};Normalised events", 60, 0.4, 1.6);
    TH1F hT("hT", "", 60, 0.4, 1.6);
    for (int i = 0; i < 40000; ++i) {
        const double truth = rng.Gaus(1.0, 0.07);
        hT.Fill(truth);
        // reco = truth x detector response, with a low tail (FSR, neutrinos)
        double reco = truth * rng.Gaus(0.97, 0.09);
        if (rng.Uniform() < 0.12) reco -= std::abs(rng.Gaus(0, 0.25));
        hR.Fill(reco);
    }
    hR.Scale(1.0 / hR.Integral());
    hT.Scale(1.0 / hT.Integral());

    TF1 core("core", "gaus", 0.80, 1.10);
    hR.Fit(&core, "RQ0N");
    core.SetRange(core.GetParameter(1) - 1.5 * core.GetParameter(2),
                  core.GetParameter(1) + 1.5 * core.GetParameter(2));
    hR.Fit(&core, "RQ0N");

    hR.SetMaximum(std::max(hR.GetMaximum(), hT.GetMaximum()) * 1.35);
    hR.SetMarkerStyle(20);
    hR.Draw("E");
    hT.SetLineColor(kGreen + 2);
    hT.SetLineWidth(3);
    hT.SetLineStyle(1);
    hT.Draw("HIST SAME");
    core.SetLineColor(kRed + 1);
    core.SetLineWidth(3);
    core.Draw("SAME");

    TLegend leg(0.60, 0.60, 0.92, 0.80);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    leg.SetTextSize(0.035);
    leg.AddEntry(&hR, "R (reco)", "pe");
    leg.AddEntry(&core, "Gaussian core fit", "l");
    leg.AddEntry(&hT, "R_{truth} = p_{T}^{truth} / p_{T}^{ref}", "l");
    leg.Draw();
    TLatex t;
    t.SetNDC();
    t.SetTextSize(0.035);
    t.DrawLatex(0.20, 0.72, Form("#mu = %.3f,  #sigma = %.3f", core.GetParameter(1), core.GetParameter(2)));
    t.DrawLatex(0.20, 0.66, "#sigma_{insitu} = #sqrt{#sigma^{2}_{R} - #sigma^{2}_{truth}}");
    t.DrawLatex(0.20, 0.60, "R_{JES} = #mu_{data} / #mu_{MC}");
    t.DrawLatex(0.20, 0.54, "R_{JER} = #sigma_{insitu} / #sigma^{MC}_{insitu}");
    Label(0.20, 0.87);
    c.SaveAs("figures/sketch_response.png");
}

// ---------------------------------------------------------------------------
// 3. Bootstrap: fitted mu of the correction vs the nominal in each replica.
// ---------------------------------------------------------------------------
void Bootstrap(TRandom3& rng) {
    TCanvas c("c_boot", "", 700, 650);
    c.SetLeftMargin(0.20);
    c.SetRightMargin(0.06);
    const int n = 200;
    const double rho = 0.6;
    TGraph g(n);
    for (int i = 0; i < n; ++i) {
        const double a = rng.Gaus(), b = rng.Gaus();
        const double xn = 0.980 + 0.004 * a;
        const double xc = 0.995 + 0.004 * (rho * a + std::sqrt(1 - rho * rho) * b);
        g.SetPoint(i, xn, xc);
    }
    TH2F frame("frame", ";#mu_{nominal} (one replica);#mu_{correction} (same replica)", 10, 0.964, 0.996, 10,
               0.980, 1.022);
    frame.GetYaxis()->SetTitleOffset(1.9);
    frame.GetXaxis()->SetNdivisions(505);
    frame.GetYaxis()->SetNdivisions(505);
    frame.Draw();
    g.SetMarkerStyle(20);
    g.SetMarkerSize(0.8);
    g.SetMarkerColor(kAzure + 2);
    g.Draw("P SAME");
    TLatex t;
    t.SetNDC();
    t.SetTextSize(0.037);
    t.DrawLatex(0.25, 0.81, Form("%d replicas, same Poisson weight per event", n));
    t.DrawLatex(0.25, 0.76,
                Form("Pearson #rho = %.2f  #rightarrow  correlations.txt", g.GetCorrelationFactor()));
    Label(0.25, 0.87);
    c.SaveAs("figures/sketch_bootstrap.png");
}

// ---------------------------------------------------------------------------
// 4. Double ratio vs pT, DB (200-450) + MJB (450-1200), with the four error
//    treatments side by side, for a response observable (left: bootstrap and
//    full agree) and a mass observable (right: the bootstrap is larger, since
//    it includes the background subtraction that the fit errors leave out).
// ---------------------------------------------------------------------------
void DoubleRatioPanel(TRandom3& rng, const char* name, const char* ytitle, const double* unc,
                      const double* full, const double* corr, const double* boot, bool withLegend) {
    const int n = 5;
    const double pt[n] = {250, 375, 525, 750, 1050};
    const double hw[n] = {50, 75, 75, 150, 150};
    double y[n];
    for (int i = 0; i < n; ++i) y[i] = 1.0 + rng.Gaus(0, 0.006);

    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.03);
    gPad->SetTopMargin(0.20);
    TH1F* frame =
        new TH1F(Form("frame_%s", name), Form(";Large-R jet p_{T} [GeV];%s", ytitle), 10, 200, 1200);
    frame->SetMinimum(0.90);
    frame->SetMaximum(1.12);
    frame->GetYaxis()->SetTitleOffset(1.5);
    frame->Draw();

    const double* errs[4] = {unc, full, corr, boot};
    const int cols[4] = {kGray + 1, kRed + 1, kAzure + 2, kBlack};
    const int mks[4] = {24, 20, 22, 21};
    const double off[4] = {-0.33, -0.11, 0.11, 0.33};
    TGraphErrors* g[4];
    for (int k = 0; k < 4; ++k) {
        g[k] = new TGraphErrors(n);
        for (int i = 0; i < n; ++i) {
            g[k]->SetPoint(i, pt[i] + off[k] * hw[i], y[i]);
            g[k]->SetPointError(i, 0, errs[k][i]);
        }
        g[k]->SetLineColor(cols[k]);
        g[k]->SetMarkerColor(cols[k]);
        g[k]->SetMarkerStyle(mks[k]);
        g[k]->SetMarkerSize(1.1);
        g[k]->SetLineWidth(2);
        g[k]->Draw("P SAME");
    }
    TLine* one = new TLine(200, 1, 1200, 1);
    one->SetLineStyle(2);
    one->Draw();
    TLine* split = new TLine(450, 0.90, 450, 1.08);
    split->SetLineStyle(3);
    split->SetLineColor(kGray + 2);
    split->Draw();
    TLatex t;
    t.SetTextSize(0.040);
    t.SetTextColor(kGray + 2);
    t.DrawLatex(255, 0.915, "DB");
    t.DrawLatex(720, 0.915, "MJB");

    if (withLegend) {
        TLegend* leg = new TLegend(0.15, 0.81, 0.97, 0.97);
        leg->SetNColumns(2);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.040);
        leg->AddEntry(g[0], "uncorrelated", "pe");
        leg->AddEntry(g[1], "full (fit errors, bootstrap #rho)", "pe");
        leg->AddEntry(g[2], "corr-aware (#rho = 1)", "pe");
        leg->AddEntry(g[3], "bootstrap (spread of DR)", "pe");
        leg->Draw();
    }
}

void DoubleRatio(TRandom3& rng) {
    TCanvas c("c_dr", "", 1300, 620);
    TPad* top = new TPad("top", "", 0, 0.86, 1, 1);
    top->Draw();
    TPad* left = new TPad("left", "", 0, 0, 0.5, 0.86);
    TPad* right = new TPad("right", "", 0.5, 0, 1, 0.86);
    left->Draw();
    right->Draw();

    //                       DB 200-300, 300-450, MJB 450-600, 600-900, 900-1200
    const double rUnc[5]  = {0.028, 0.032, 0.026, 0.034, 0.046};
    const double rFull[5] = {0.016, 0.019, 0.015, 0.020, 0.028};
    const double rCorr[5] = {0.006, 0.007, 0.006, 0.008, 0.011};
    const double rBoot[5] = {0.017, 0.020, 0.015, 0.021, 0.030};
    const double mUnc[5]  = {0.034, 0.038, 0.030, 0.036, 0.048};
    const double mFull[5] = {0.014, 0.016, 0.012, 0.015, 0.020};
    const double mCorr[5] = {0.006, 0.007, 0.005, 0.006, 0.009};
    const double mBoot[5] = {0.032, 0.036, 0.027, 0.034, 0.045};

    left->cd();
    DoubleRatioPanel(rng, "resp", "DR(b-JES) = correction / nominal", rUnc, rFull, rCorr, rBoot, true);
    TLatex t;
    t.SetNDC();
    t.SetTextFont(42);
    t.SetTextSize(0.042);
    t.DrawLatex(0.18, 0.74, "Response: bootstrap #approx full");
    right->cd();
    DoubleRatioPanel(rng, "mass", "DR(b-JMDP) = correction / nominal", mUnc, mFull, mCorr, mBoot, false);
    t.DrawLatex(0.18, 0.74, "Mass: bootstrap > full");
    t.SetTextSize(0.036);
    t.DrawLatex(0.18, 0.69, "(fit errors miss the background subtraction)");

    top->cd();
    Label(0.02, 0.45, "", 0.32);
    c.SaveAs("figures/sketch_double_ratio.png");
}

}  // namespace Sketch

void make_readme_figures() {
    gROOT->SetBatch(kTRUE);
    SetAtlasStyle();
    gStyle->SetOptFit(0);
    gSystem->Exec("mkdir -p figures");
    TRandom3 rng(2026);
    Sketch::MassStack(rng);
    Sketch::BkgComposition();
    Sketch::MassFit(rng);
    Sketch::MassExtracted(rng);
    Sketch::MassSummary();
    Sketch::ResponseStack(rng);
    Sketch::Response(rng);
    Sketch::ResponseSummary();
    Sketch::Bootstrap(rng);
    Sketch::DoubleRatio(rng);
}
