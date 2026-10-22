/**
 * =============================================================================================
 *  CommonUtils.h — Shared utilities for EasyBjets (used by both ZbbGamma and ZbbJets)
 * =============================================================================================
 * ===========================================
 * @author:    Thiziri AMEZZA
 * @email:     thiziri.amezza@cern.ch
 * @date:      10 March 2026 
 * ===========================================
 *
 * This header collects common utility tools used across both channels of the
 * EasyBjets in-situ calibration framework (Direct Balance in ZbbGamma, Multijets
 * Balance in ZbbJets). It is meant to be included in the macros/ directory of each channel.
 *
 * Contents:
 *   Log             coloured [LEVEL] [context] messages, level set with Log::SetLevel
 *   SafeTFile       TFile that closes itself; CanvasGuard, same idea for a TCanvas
 *   SaveCanvasPdf   SaveAs(pdf) plus a .png copy when EASYBJETS_SAVE_PNG is set
 *   FileUtils       mkdir -p, FileExists
 *   HistUtils       detached histograms, MC stat band, data point style, ratio
 *   PadUtils        the upper/lower ratio-plot layout used by all stack plots
 *   Validate        tree/branch lookups that complain instead of segfaulting
 *   CutFlowUtils    cut names per channel, cutflow histograms of all samples
 *   BinConfig       pT bins and mass cuts used by the response fits
 *   FitUtils        the iterative Gaussian fit of the response peak and the
 *                   alternative models (Bukin, Crystal Ball, double CB, double
 *                   Gaussian), plus the JES / JER ratio error propagation
 *   FitOutputUtils  fit parameter text file written next to each fit plot
 *
 * Environment variables read here:
 *   EASYBJETS_SAVE_PNG              also save canvases as png
 *   EASYBJETS_FIT_MODEL             gauss (default) | all | bukin | cb | dscb | dgauss
 *   EASYBJETS_FIT_MODE              scan (default) | scan_dense | scan_log | legacy | legacy7grid
 *   EASYBJETS_FIT_NSIGMA            window half-width in legacy mode (2.0)
 *   EASYBJETS_FIT_NSIGMA_LIST       candidate windows in scan mode, "1.0,1.5:1.0,2.0"
 *   EASYBJETS_FIT_MODEL_NSIGMA_LO/HI  window of the alternative models, in sigma (3, 2)
 *   EASYBJETS_FIT_ERR_SCALE         pdg: scale fit errors by sqrt(chi2/ndf) when > 1
 *   EASYBJETS_SEED_MODE             truncated (default) | raw, see SeedMode
 *   EASYBJETS_RESPONSE_FIT_PTMIN    per-pT-bin floor of the fit window, R > ptmin / pT_lo
 * 
 * ==========================================================================
 * @version:  1.0
 * ========================================================================== 
 */

#ifndef COMMON_UTILS_H
#define COMMON_UTILS_H

#include <string>
#include <vector>
#include <memory>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cctype>
#include <cstdlib>
#include <functional>
#include <algorithm>
#include <chrono>
#include <cassert>
#include <TString.h>

#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TH2F.h>
#include <THStack.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TDirectory.h>
#include <TAxis.h>

#include "VariantConfig.h"

// ============================================================================
// 1. STANDARDIZED LOGGING
// ============================================================================
namespace Log {

    enum Level { DEBUG = 0, INFO = 1, WARNING = 2, ERROR = 3, SILENT = 4 };

    inline Level& GlobalLevel() {
        static Level level = INFO;
        return level;
    }

    inline void SetLevel(Level lvl) { GlobalLevel() = lvl; }

    inline const char* LevelName(Level lvl) {
        switch (lvl) {
            case DEBUG:   return "DEBUG";
            case INFO:    return "INFO";
            case WARNING: return "WARN";
            case ERROR:   return "ERROR";
            default:      return "???";
        }
    }

    inline const char* LevelColor(Level lvl) {
        switch (lvl) {
            case DEBUG:   return "\033[36m";
            case INFO:    return "\033[32m";
            case WARNING: return "\033[33m";
            case ERROR:   return "\033[31m";
            default:      return "\033[0m";
        }
    }

    inline void Print(Level lvl, const std::string& context, const std::string& msg) {
        if (lvl < GlobalLevel()) return;
        std::ostream& out = (lvl >= WARNING) ? std::cerr : std::cout;
        out << LevelColor(lvl) << "[" << LevelName(lvl) << "]"
            << "\033[0m";
        if (!context.empty()) out << " [" << context << "]";
        out << " " << msg << std::endl;
    }

    inline void debug(const std::string& ctx, const std::string& msg)   { Print(DEBUG, ctx, msg); }
    inline void info(const std::string& ctx, const std::string& msg)    { Print(INFO, ctx, msg); }
    inline void warning(const std::string& ctx, const std::string& msg) { Print(WARNING, ctx, msg); }
    inline void error(const std::string& ctx, const std::string& msg)   { Print(ERROR, ctx, msg); }

    // Simple timer for performance profiling
    class Timer {
        std::chrono::steady_clock::time_point start_;
        std::string label_;
    public:
        Timer(const std::string& label) : label_(label), start_(std::chrono::steady_clock::now()) {
            debug("Timer", "Started: " + label);
        }
        ~Timer() {
            auto elapsed = std::chrono::steady_clock::now() - start_;
            double sec = std::chrono::duration<double>(elapsed).count();
            info("Timer", label_ + " completed in " + std::to_string(sec) + " s");
        }
    };  

}  // namespace Log


// ============================================================================
// 2. RAII WRAPPERS FOR ROOT OBJECTS
// ============================================================================

// SafeTFile f(path); if (!f) ...; TH1F* h = f.Get<TH1F>("h");
// A zombie file counts as not open. Move-only.
class SafeTFile {
    TFile* file_ = nullptr;
public:
    SafeTFile() = default;

    SafeTFile(const std::string& path, const std::string& mode = "READ") {
        file_ = TFile::Open(path.c_str(), mode.c_str());
        if (file_ && file_->IsZombie()) {
            delete file_;
            file_ = nullptr;
        }
    }

    ~SafeTFile() { Close(); }

    SafeTFile(SafeTFile&& o) noexcept : file_(o.file_) { o.file_ = nullptr; }
    SafeTFile& operator=(SafeTFile&& o) noexcept {
        if (this != &o) { Close(); file_ = o.file_; o.file_ = nullptr; }
        return *this;
    }
    SafeTFile(const SafeTFile&) = delete;
    SafeTFile& operator=(const SafeTFile&) = delete;

    void Close() {
        if (file_) { file_->Close(); delete file_; file_ = nullptr; }
    }

    explicit operator bool() const { return file_ != nullptr; }
    TFile* Get() const { return file_; }
    TFile* operator->() const { return file_; }

    template<typename T>
    T* Get(const std::string& name) const {
        if (!file_) return nullptr;
        return static_cast<T*>(file_->Get(name.c_str()));
    }
};

// Helper to detach a histogram from its file so it can be written to another file.
template<typename T>
inline T* DetachFromDirectory(T* hist) {
    if (hist) hist->SetDirectory(nullptr);
    return hist;
}

// Check if we should save .png versions of plots.
inline bool EasyBjetsSavePng() {
    static bool cached = [](){
        const char* env = std::getenv("EASYBJETS_SAVE_PNG");
        if (!env || !*env) return false;
        std::string s(env);
        for (auto& c : s) c = (char)std::tolower((unsigned char)c);
        return !(s == "0" || s == "false" || s == "no" || s == "off");
    }();
    return cached;
}

// Use this rather than canvas->SaveAs() so every plot gets the png copy
// when EASYBJETS_SAVE_PNG is on.
inline void SaveCanvasPdf(TCanvas* c, const std::string& pdfPath) {
    if (!c || pdfPath.empty()) return;
    c->Update();
    c->SaveAs(pdfPath.c_str());
    if (EasyBjetsSavePng()) {
        std::string png = pdfPath;
        const std::string ext = ".pdf";
        if (png.size() >= ext.size() &&
            png.compare(png.size() - ext.size(), ext.size(), ext) == 0)
            png.replace(png.size() - ext.size(), ext.size(), ".png");
        else
            png += ".png";
        c->SaveAs(png.c_str());
    }
}

// RAII wrapper for a TCanvas.
// Saves the canvas to savePath (if set) and deletes it on scope exit.
class CanvasGuard {
    TCanvas* canvas_ = nullptr;
    std::string savePath_;

    void DoSave() { SaveCanvasPdf(canvas_, savePath_); }
public:
    CanvasGuard(const std::string& name, const std::string& title,
                int width, int height, const std::string& savePath = "")
        : savePath_(savePath)
    {
        canvas_ = new TCanvas(name.c_str(), title.c_str(), width, height);
    }

    ~CanvasGuard() {
        if (canvas_) {
            DoSave();
            delete canvas_;
        }
    }

    CanvasGuard(CanvasGuard&& o) noexcept : canvas_(o.canvas_), savePath_(std::move(o.savePath_)) {
        o.canvas_ = nullptr;
    }
    CanvasGuard(const CanvasGuard&) = delete;
    CanvasGuard& operator=(const CanvasGuard&) = delete;

    TCanvas* Canvas() { return canvas_; }
    TCanvas* operator->() { return canvas_; }

    void Save() { DoSave(); }
    void SetSavePath(const std::string& path) { savePath_ = path; }
};


// ============================================================================
// 3. FILE SYSTEM UTILITIES
// ============================================================================
namespace FileUtils {

    inline void MakeDir(const std::string& path) {
        gSystem->Exec(("mkdir -p " + path).c_str());
    }

    inline bool FileExists(const std::string& path) {
        if (path.find("root://") == 0) return true;
        return !gSystem->AccessPathName(path.c_str());
    }

    inline SafeTFile OpenRootFile(const std::string& path, const std::string& context = "",
                                   bool required = true) {
        if (!FileExists(path)) {
            if (required)
                Log::error(context, "File not found: " + path);
            else
                Log::debug(context, "Optional file not found: " + path);
            return SafeTFile();
        }
        SafeTFile f(path, "READ");
        if (!f) {
            if (required)
                Log::error(context, "Failed to open: " + path);
            else
                Log::info(context, "Optional file not available: " + path);
        }
        return f;
    }

}  // namespace FileUtils


// ============================================================================
// 4. HISTOGRAM UTILITIES
// ============================================================================
namespace HistUtils {

    // Fixed-size array of owned TH1F, one per sample. operator[] gives the
    // raw pointer because that is what ROOT wants everywhere.
    template<int N>
    struct HistArray {
        std::unique_ptr<TH1F> arr[N];

        TH1F* operator[](int i) const { return arr[i].get(); }
        void Set(int i, TH1F* h) { arr[i].reset(h); }
        TH1F* Raw(int i) const { return arr[i].get(); }

        struct Iter {
            const std::unique_ptr<TH1F>* ptr;
            TH1F* operator*() const { return ptr->get(); }
            Iter& operator++() { ++ptr; return *this; }
            bool operator!=(const Iter& o) const { return ptr != o.ptr; }
        };
        Iter begin() const { return {arr}; }
        Iter end()   const { return {arr + N}; }
    };

    inline void Normalize(std::vector<TH1F*>& histograms) {
        for (auto& h : histograms) {
            if (h && h->Integral() > 0) {
                h->Scale(1.0 / h->Integral());
            }
        }
    }

    inline void Normalize(TH1F* h) {
        if (h && h->Integral() > 0) h->Scale(1.0 / h->Integral());
    }

    // Create histogram detached from any TDirectory 
    inline TH1F* CreateDetached(const std::string& name, const std::string& title,
                                 int nbins, double xmin, double xmax) {
        TH1F* h = new TH1F(name.c_str(), title.c_str(), nbins, xmin, xmax);
        h->SetDirectory(nullptr);
        h->Sumw2();
        return h;
    }

    inline TH1D* CreateDetachedD(const std::string& name, const std::string& title,
                                  int nbins, double xmin, double xmax) {
        TH1D* h = new TH1D(name.c_str(), title.c_str(), nbins, xmin, xmax);
        h->SetDirectory(nullptr);
        h->Sumw2();
        return h;
    }

    // Clone and detach from directory
    template<typename T>
    inline T* CloneDetached(const T* src, const std::string& newName) {
        if (!src) return nullptr;
        T* clone = static_cast<T*>(src->Clone(newName.c_str()));
        clone->SetDirectory(nullptr);
        return clone;
    }

    // Initialize sample-indexed histogram vector
    inline std::vector<TH1F*> InitSampleHistograms(
            const std::vector<std::string>& histNames,
            const std::string& variable, const std::string& title,
            int nbins, double xmin, double xmax) {
        std::vector<TH1F*> hists;
        for (const auto& name : histNames) {
            hists.push_back(CreateDetached("h_" + name + "_" + variable, title, nbins, xmin, xmax));
        }
        return hists;
    }

    inline void ApplySampleColors(std::vector<TH1F*>& hists,
                                   const std::vector<int>& colors,
                                   const std::vector<int>& fillStyles) {
        for (size_t i = 0; i < hists.size(); ++i) {
            if (!hists[i]) continue;
            if (i < colors.size()) {
                hists[i]->SetLineColor(colors[i]);
                hists[i]->SetFillColor(colors[i]);
            }
            if (i < fillStyles.size()) {
                hists[i]->SetFillStyle(fillStyles[i]);
            }
        }
    }

    // ===========================================================
    // MC stat-uncertainty band — canonical Run-2 reference style
    // ===========================================================
    // Use these helpers for ANY stack plot with a "MC Stat. Unc." band.
    // Fixed at: kBlack hatch pattern 3254at 30% alpha.
    // Match this style across the framework to keep cutflow / mass-stack /
    // response-stack plots visually consistent. Legend label: "MC Stat. Unc."
    // ===========================================================
    inline TGraphAsymmErrors* MakeUncertaintyBand(TH1F* hMC) {
        if (!hMC) return nullptr;
        TGraphAsymmErrors* band = new TGraphAsymmErrors(hMC);
        for (int i = 0; i < band->GetN(); ++i) {
            double error = hMC->GetBinError(i + 1);
            band->SetPointEYhigh(i, error);
            band->SetPointEYlow(i, error);
        }
        band->SetFillColorAlpha(kBlack, 0.3);  
        band->SetFillStyle(3254);              
        band->SetMarkerSize(0);                
        band->SetLineWidth(0);                 
        return band;
    }

    // Build ratio uncertainty band (relative errors around y=1, for ratio panels)
    inline TGraphAsymmErrors* MakeRatioUncertaintyBand(TH1F* hMC) {
        if (!hMC) return nullptr;
        TGraphAsymmErrors* band = new TGraphAsymmErrors(hMC);
        for (int i = 0; i < band->GetN(); ++i) {
            double content = hMC->GetBinContent(i + 1);
            double error   = hMC->GetBinError(i + 1);
            double relErr  = (content > 0) ? error / content : 0;
            band->SetPointY(i, 1.0);
            band->SetPointEYhigh(i, relErr);
            band->SetPointEYlow(i, relErr);
        }
        band->SetFillColorAlpha(kBlack, 0.3);  
        band->SetFillStyle(3254);              
        band->SetMarkerSize(0);
        band->SetLineWidth(0);
        return band;
    }

    inline void StyleDataPoints(TH1F* hData) {
        if (!hData) return;
        hData->SetMarkerStyle(20);   
        hData->SetMarkerSize(0.8);   
        hData->SetMarkerColor(kBlack);
        hData->SetLineColor(kBlack);
        hData->SetLineWidth(1);
    }


    inline TH1F* MakeRatio(TH1F* hData, TH1F* hMC, const std::string& name = "h_ratio") {
        if (!hData || !hMC) return nullptr;
        TH1F* ratio = CloneDetached(hData, name);
        ratio->SetLineColor(kBlack);
        ratio->SetMarkerStyle(21);
        ratio->Divide(hMC);
        return ratio;
    }

    inline void DeleteAll(std::vector<TH1F*>& hists) {
        for (auto h : hists) delete h;
        hists.clear();
    }

}  // namespace HistUtils


// ============================================================================
// 5. CANVAS / PAD UTILITIES
// ============================================================================
namespace PadUtils {

    struct TwoPadLayout {
        TPad* upper = nullptr;
        TPad* lower = nullptr;
    };

    inline TwoPadLayout CreateRatioPads(TCanvas* canvas,
                                         double splitFrac = 0.30,
                                         const std::string& padPrefix = "pad") {
        TwoPadLayout layout;
        canvas->cd();

        layout.upper = new TPad((padPrefix + "1").c_str(), "upper", 0, splitFrac, 1, 1);
        layout.upper->SetBottomMargin(0.02);
        layout.upper->SetTopMargin(0.06);
        layout.upper->SetLeftMargin(0.14);
        layout.upper->SetRightMargin(0.04);
        layout.upper->Draw();

        canvas->cd();
        layout.lower = new TPad((padPrefix + "2").c_str(), "lower", 0, 0, 1, splitFrac);
        layout.lower->SetTopMargin(0.02);
        layout.lower->SetBottomMargin(0.35);
        layout.lower->SetLeftMargin(0.14);
        layout.lower->SetRightMargin(0.04);
        layout.lower->SetGridy();
        layout.lower->Draw();

        return layout;
    }

    inline void StyleRatioHist(TH1F* h, const std::string& xTitle = "",
                                double yMin = 0.5, double yMax = 1.5) {
        if (!h) return;
        h->SetMinimum(yMin);
        h->SetMaximum(yMax);
        h->GetYaxis()->SetTitle("Data/MC");
        h->GetYaxis()->SetTitleSize(0.15);
        h->GetYaxis()->SetLabelSize(0.1);
        h->GetYaxis()->SetTitleOffset(0.3);
        h->GetYaxis()->SetNdivisions(505);

        if (!xTitle.empty())
            h->GetXaxis()->SetTitle(xTitle.c_str());
        h->GetXaxis()->SetTitleSize(0.12);
        h->GetXaxis()->SetLabelSize(0.12);
        h->GetXaxis()->SetTitleOffset(1.25);
    }

    inline TLine* DrawRefLine(TPad* pad, double y = 1.0) {
        TLine* line = new TLine(pad->GetUxmin(), y, pad->GetUxmax(), y);
        line->SetLineColor(kBlack);
        line->SetLineWidth(2);
        line->SetLineStyle(2);
        line->Draw("same");
        return line;
    }

}  // namespace PadUtils


// ============================================================================
// 6. INPUT VALIDATION
// ============================================================================
namespace Validate {

    inline TTree* GetTree(TFile* file, const std::string& treeName,
                           const std::string& context = "") {
        if (!file) {
            Log::error(context, "Null file when looking for tree '" + treeName + "'");
            return nullptr;
        }
        TTree* tree = static_cast<TTree*>(file->Get(treeName.c_str()));
        if (!tree) {
            Log::error(context, "Tree '" + treeName + "' not found in " + file->GetName());
        }
        return tree;
    }

    template<typename T = TH1F>
    inline T* GetHist(TFile* file, const std::string& histName,
                       const std::string& context = "") {
        if (!file) return nullptr;
        T* h = static_cast<T*>(file->Get(histName.c_str()));
        if (!h) {
            Log::error(context, "Histogram '" + histName + "' not found in " + file->GetName());
        }
        return h;
    }

    inline bool HasBranch(TTree* tree, const std::string& branchName) {
        return tree && tree->GetBranch(branchName.c_str()) != nullptr;
    }

    template<typename T>
    inline bool SafeSetBranch(TTree* tree, const std::string& branchName, T* addr,
                               const std::string& context = "") {
        if (!HasBranch(tree, branchName)) {
            Log::warning(context, "Branch '" + branchName + "' not found");
            return false;
        }
        tree->SetBranchAddress(branchName.c_str(), addr);
        return true;
    }

}  // namespace Validate


// ============================================================================
// 7. CUTFLOW UTILITIES
// ============================================================================
namespace CutFlowUtils {

    // Returns cut names that match the analysis display convention.
    // Pass cfg.channel ("MJB" or "DB") and, optionally, the working-point label.
    // label string so the b-tagger entry reads "D_{xbb} 0p46% WP" etc.
    // When wp is omitted (e.g. just to get nCuts) the entry reads "D_{xbb}".
    inline std::vector<std::string> DefaultCutNames(const std::string& channel = "MJB",
                                                     const std::string& wp = "") {
        const std::string btagDB  = wp.empty() ? "D_{xbb}"
                                               : "D_{xbb} " + wp;
        const std::string btagMJB = wp.empty() ? "D_{xbb}"
                                               : "D_{xbb} " + wp;
        // DB — 12 cuts
        if (channel == "DB") {
            return {
                "Initial",                                    
                "p_{T}^{#gamma} > 150 GeV",                   
                "|#eta^{#gamma}| < 1.37",                     
                "p_{T}^{J} > 200 GeV",                        
                "#DeltaR(#gamma, J) > 0.2",                   
                "|#Delta#phi(#gamma, J)| > 2.8",              
                "JVT (p_{T}^{j} < 60, |#eta^{j}| < 2.4)",     
                "#DeltaR(#gamma, j) > 0.4",                   
                "#DeltaR(J, j) > 1.4",                        
                "Veto",                                       
                btagDB,                                       
                "M^{J}"                                       
            };
        }
        // MJB - 8 cuts
        return {
            "Initial",                                    // Cut 0
            btagMJB,                                      // Cut 1
            "p_{T}^{J} > 450",                            // Cut 2
            "|#eta^{J}| < 1.2",                           // Cut 3
            "Small-R p_{T} & #eta",                       // Cut 4
            "JVT > 0.59",                                 // Cut 5
            "|#alpha - #pi| < 0.3",                       // Cut 6
            "#Delta R (J1, DB) > 1.5"                     // Cut 7
        };
    }

    // Normalization bin for the cutflow.
    // This is the bin at which the main-background MC cutflow is normalized to
    // data before computing downstream ratios. 
    inline int NormBin(const std::string& channel = "DB") {
        if (channel == "MJB") return 1;  
        return 10;                        
    }

    inline int NCuts(const std::string& channel = "DB") {
        return static_cast<int>(DefaultCutNames(channel).size());
    }

    // Flat (mass-independent) sideband scale factor, data/MC yield ratio.
    inline double ComputeSidebandSF(TH1F* hData, TH1F* hBkg,
                                     double sb1_lo, double sb1_hi,
                                     double sb2_lo, double sb2_hi) {
        if (!hData || !hBkg) return 1.0;
        int b1lo = hBkg->FindBin(sb1_lo);
        int b1hi = hBkg->FindBin(sb1_hi) - 1;
        int b2lo = hBkg->FindBin(sb2_lo);
        int b2hi = hBkg->FindBin(sb2_hi) - 1;
        double bkgSum  = hBkg ->Integral(b1lo, b1hi) + hBkg ->Integral(b2lo, b2hi);
        double dataSum = hData->Integral(b1lo, b1hi) + hData->Integral(b2lo, b2hi);
        if (bkgSum <= 0) return 1.0;
        return dataSum / bkgSum;
    }

    inline void PrintEfficiencyTable(
            const std::string& label,
            const std::vector<std::string>& wpLabels,
            const std::vector<double>& dataEff,
            const std::vector<double>& signalEff,
            const std::vector<double>& bkgEff) {
        std::cout << "\n[" << label << "] Efficiency Summary:" << std::endl;
        std::cout << "+---------------------+---------------------+---------------------+---------------------+" << std::endl;
        std::cout << "| WP                  | Data Final Eff.     | Signal Final Eff.   | Bkg Final Eff.      |" << std::endl;
        std::cout << "+---------------------+---------------------+---------------------+---------------------+" << std::endl;
        for (size_t i = 0; i < wpLabels.size(); ++i) {
            std::cout << "| " << std::left << std::setw(8) << wpLabels[i] << " | "
                      << std::right << std::setw(16) << std::scientific << std::setprecision(4) << dataEff[i] << " | "
                      << std::setw(16) << signalEff[i] << " | "
                      << std::setw(16) << bkgEff[i] << " |" << std::endl;
        }
        std::cout << "+---------------------+---------------------+---------------------+---------------------+" << std::endl;
    }

    // CutFlowSet — RAII container for per-sample cutflow histograms.
    class CutFlowSet {
        TH1F* histos_[SampleConfig::kNSamples] = {nullptr};
        bool  loaded_[SampleConfig::kNSamples] = {false};

    public:
        CutFlowSet() = default;
        ~CutFlowSet() { Clear(); }

        CutFlowSet(CutFlowSet&& o) noexcept {
            for (int i = 0; i < SampleConfig::kNSamples; ++i) {
                histos_[i] = o.histos_[i];
                loaded_[i] = o.loaded_[i];
                o.histos_[i] = nullptr;
                o.loaded_[i] = false;
            }
        }
        CutFlowSet& operator=(CutFlowSet&& o) noexcept {
            if (this != &o) {
                Clear();
                for (int i = 0; i < SampleConfig::kNSamples; ++i) {
                    histos_[i] = o.histos_[i];
                    loaded_[i] = o.loaded_[i];
                    o.histos_[i] = nullptr;
                    o.loaded_[i] = false;
                }
            }
            return *this;
        }
        CutFlowSet(const CutFlowSet&)            = delete;
        CutFlowSet& operator=(const CutFlowSet&) = delete;

        void Clear() {
            for (int i = 0; i < SampleConfig::kNSamples; ++i) {
                if (histos_[i]) { delete histos_[i]; histos_[i] = nullptr; }
                loaded_[i] = false;
            }
        }

        bool IsLoaded(int s) const {
            return s >= 0 && s < SampleConfig::kNSamples && loaded_[s];
        }

        TH1F* operator[](int s) const {
            return (s >= 0 && s < SampleConfig::kNSamples) ? histos_[s] : nullptr;
        }

        void GetRawArrays(TH1F* h_out[], bool loaded_out[]) const {
            for (int i = 0; i < SampleConfig::kNSamples; ++i) {
                h_out[i]      = histos_[i];
                loaded_out[i] = loaded_[i];
            }
        }

        // Take ownership  
        void Set(int s, TH1F* h) {
            if (s < 0 || s >= SampleConfig::kNSamples) return;
            if (histos_[s]) delete histos_[s];
            histos_[s] = h;
            loaded_[s] = (h != nullptr);
        }
    };

    // ------------------------------------------------------------------------
    // LoadAllCutflows - open the slimmed ntuple for each sample 
    // and read its "h_cutflow" histogram into a CutFlowSet.
    // ------------------------------------------------------------------------
    inline CutFlowSet LoadAllCutflows(const VariantConfig& cfg,
                                       const std::string& wp,
                                       const std::string& label = "CutFlow") {
        CutFlowSet set;
        auto fileNames = SampleConfig::FileNames(cfg);
        std::string fileWP = PathBuilder::WPtoFileWP(wp);

        for (int s = 0; s < SampleConfig::kNSamples; ++s) {
            std::string path = PathBuilder::SamplePath(cfg, fileNames[s], wp);
            SafeTFile f(path, "READ");
            if (!f) {
                Log::warning(label, "Cannot open " + path);
                continue;
            }
            
            TH1* hRaw = f.Get<TH1>("h_cutflow");
            if (!hRaw) {
                Log::warning(label, "h_cutflow not found in " + path);
                continue;
            }

            std::string uniqueName = "h_cutflow_" + cfg.key + "_" + fileWP
                                   + "_" + fileNames[s];

            int nbins = hRaw->GetNbinsX();
            TH1F* hNew = new TH1F(uniqueName.c_str(),
                                   hRaw->GetTitle(),
                                   nbins,
                                   hRaw->GetXaxis()->GetXmin(),
                                   hRaw->GetXaxis()->GetXmax());
            hNew->SetDirectory(nullptr);  // detach from file
            for (int b = 1; b <= nbins; ++b) {
                hNew->SetBinContent(b, hRaw->GetBinContent(b));
                hNew->SetBinError  (b, hRaw->GetBinError(b));
                const char* lbl = hRaw->GetXaxis()->GetBinLabel(b);
                if (lbl && *lbl) hNew->GetXaxis()->SetBinLabel(b, lbl);
            }
            set.Set(s, hNew);
        }
        return set;
    }

}  // namespace CutFlowUtils


// ============================================================================
// 8. SHARED BIN CONFIGURATIONS: pT bins and mass cuts of the response fits
// ============================================================================
// Two edge arrays per channel. ptBins_* is used for the event selection and
// has an open last bin (1e9) so nothing above the displayed range is lost.
// ptBins_*_display has the physical upper edge and is for axes, labels and
// filenames.
namespace BinConfig {

    constexpr double kPtBinOpenEnd = 1.0e9;

    // MJB, three pT bins
    constexpr int nPtBins_3 = 3;
    constexpr double ptBins_3[4]         = {450, 600, 900, kPtBinOpenEnd};
    constexpr double ptBins_3_display[4] = {450, 600, 900, 1200};
    inline const char* PtBinLabel_3(int i) {
        static const char* labels[3] = {
            "450 < p_{T} < 600 GeV",
            "600 < p_{T} < 900 GeV",
            "900 < p_{T} < 1200 GeV"
        };
        return (i >= 0 && i < 3) ? labels[i] : "???";
    }
    inline const char* PtBinFileLabel_3(int i) {
        static const char* labels[3] = {
            "450to600GeV",
            "600to900GeV",
            "900to1200GeV"
        };
        return (i >= 0 && i < 3) ? labels[i] : "???";
    }

    // DB, two pT bins
    constexpr int nPtBins_DB = 2;
    constexpr double ptBins_DB[3]         = {200, 300, kPtBinOpenEnd};
    constexpr double ptBins_DB_display[3] = {200, 300, 450};
    inline const char* PtBinLabel_DB(int i) {
        static const char* labels[2] = {
            "200 < p_{T} < 300 GeV",
            "300 < p_{T} < 450 GeV"
        };
        return (i >= 0 && i < 2) ? labels[i] : "???";
    }

    inline const char* PtBinFileLabel_DB(int i) {
        static const char* labels[2] = {
            "200to300GeV",
            "300to450GeV"
        };
        return (i >= 0 && i < 2) ? labels[i] : "???";
    }

    // mass cuts of Response_Fit_pTbins_Mcuts; the last one is open-ended
    constexpr int nMassCuts_3 = 3;
    inline const char* MassCutLabel(int i) {
        static const char* labels[3] = {"50-80GeV", "80-110GeV", "#geq110 GeV"};
        return (i >= 0 && i < 3) ? labels[i] : "???";
    }
    inline const char* MassCutFileLabel(int i) {
        static const char* labels[3] = {"50-80GeV", "80-110GeV", "MassGE110GeV"};
        return (i >= 0 && i < 3) ? labels[i] : "???";
    }
    constexpr double massCutLow_3[3]  = {50, 80, 110};
    constexpr double massCutHigh_3[3] = {80, 110, -1};

}  // namespace BinConfig


// ============================================================================
// 9. MULTITHREADING SUPPORT
// ============================================================================
namespace MTUtils {

    inline void EnableMT() {
        #ifdef R__USE_IMT
        static bool enabled = false;
        if (!enabled) {
            ROOT::EnableImplicitMT();
            enabled = true;
            Log::info("MT", "Implicit multithreading enabled");
        }
        #endif
    }

}  // namespace MTUtils


// ============================================================================
// FitUtils: peak fits of the response distributions
// ============================================================================
namespace FitUtils {

    // One fit of a peak. For the non-Gaussian models mean and sigma are the
    // peak position and core width parameters; pars holds all parameters so
    // the real curve can be redrawn with MakeModelTF1.
    struct GausFitResult {
        double amp = 0, ampErr = 0;
        double mean = 0, meanErr = 0;
        double sigma = 0, sigmaErr = 0;
        double chi2 = 0;
        int    ndf = 0;
        double chi2ndf = 0;
        double fitLo = 0, fitHi = 0;
        double nSigma = 0;          // half-width of the window in sigma (symmetric windows)
        double nSigmaLo = 0;        // asymmetric half-width below peak (down side)
        double nSigmaHi = 0;        // asymmetric half-width above peak (up side)
        bool   valid = false;
        std::string modelName = "gaus";  // gaus, bukin, crystalball, dscb, doublegauss
        std::vector<double> pars;

        double AbsDevFrom1() const { return std::abs(chi2ndf - 1.0); }

        // PDG-standard error scaling: / With EASYBJETS_FIT_ERR_SCALE=pdg the parameter 
        // errors are inflated by sqrt(chi2/ndf) when chi2/ndf > 1 (PDG practice). 
        // Only applied where results are displayed or saved, never in the ratio arrays. 
        // This is a conservative uncertainty estimate that acknowledges model mismatch.
        // Central values are unchanged.
        double errScaleFactor() const {
            static const bool kScaleOn = [](){
                const char* e = std::getenv("EASYBJETS_FIT_ERR_SCALE");
                if (!e || !*e) return false;
                std::string s(e);
                for (auto& c : s) c = (char)std::tolower((unsigned char)c);
                return (s == "pdg" || s == "sqrt" || s == "1" || s == "on" || s == "true");
            }();
            if (!kScaleOn) return 1.0;
            if (!std::isfinite(chi2ndf) || chi2ndf <= 1.0) return 1.0;
            return std::sqrt(chi2ndf);
        }
        double meanErrScaled()  const { return meanErr  * errScaleFactor(); }
        double sigmaErrScaled() const { return sigmaErr * errScaleFactor(); }
        double ampErrScaled()   const { return ampErr   * errScaleFactor(); }

        // Window selection criterion in scan mode. chi2/ndf above 1 (the fit
        // disagrees with the data) is penalised four times more than the same
        // distance below 1 (a slightly too wide window is harmless).
        double Penalty() const {
            if (!std::isfinite(chi2ndf)) return 1e9;
            if (chi2ndf <= 1.0) return (1.0 - chi2ndf) * 0.25;
            return (chi2ndf - 1.0);
        }
    };

    // Checks how much of a fit window is actually populated.
    // A window with many empty bins gives a meaningless chi2/ndf; reject it.
    inline bool HasEnoughContentInWindow(TH1* h, double lo, double hi,
                                        double minFilledFraction = 0.7) {
        if (!h) return false;
        int blo = h->FindBin(lo), bhi = h->FindBin(hi);
        int total = bhi - blo + 1;
        if (total <= 0) return false;
        int filled = 0;
        for (int b = blo; b <= bhi; ++b)
            if (h->GetBinContent(b) > 0) ++filled;
        return (double)filled / total >= minFilledFraction;
    }

    // ============================================================================
    // Alternative fit models: Bukin and Crystal Ball
    // ============================================================================
    // The response peak has an ISR tail on the low side and the mass peak is
    // slightly asymmetric after background subtraction; a Gaussian forced on
    // them gives a sigma that depends on the window. These models are tried
    // next to (or instead of) the Gaussian when EASYBJETS_FIT_MODEL is set.
    // The default is Gaussian only.
    //
    enum class FitModelSet { GaussianOnly, AllModels, BukinOnly, CrystalBallOnly, DoubleCBOnly, DoubleGaussOnly };

    inline FitModelSet CurrentFitModelSet() {
        static FitModelSet cached = [](){
            const char* env = std::getenv("EASYBJETS_FIT_MODEL");
            if (!env || !*env) return FitModelSet::GaussianOnly;
            std::string s(env);
            for (auto& c : s) c = (char)std::tolower((unsigned char)c);
            if (s == "gauss" || s == "gaussian") return FitModelSet::GaussianOnly;
            if (s == "all" || s == "bukin-cb" || s == "extended") return FitModelSet::AllModels;
            if (s == "bukin" || s == "bukinonly" || s == "bukin-only") return FitModelSet::BukinOnly;
            if (s == "crystalball" || s == "cb" || s == "crystalballonly" || s == "cb-only")
                return FitModelSet::CrystalBallOnly;
            if (s == "dscb" || s == "doublecb" || s == "double-cb" || s == "dscbonly" ||
                s == "doublesidedcb" || s == "double-sided-cb")
                return FitModelSet::DoubleCBOnly;
            if (s == "doublegauss" || s == "doublegaussian" || s == "dgauss" ||
                s == "2gauss" || s == "twogauss" || s == "double-gauss")
                return FitModelSet::DoubleGaussOnly;
            std::cerr << "[FitModel] WARN: unknown EASYBJETS_FIT_MODEL='" << env
                      << "', using gaussian-only\n";
            return FitModelSet::GaussianOnly;
        }();
        return cached;
    }

    // Short filesystem tag for the active fit model set. 
    // Output folder suffix for the alternative models (see RunEra::FitModelDirTag).
    inline std::string CurrentFitModelTag() {
        switch (CurrentFitModelSet()) {
            case FitModelSet::BukinOnly:       return "_bukin";
            case FitModelSet::CrystalBallOnly: return "_cb";
            case FitModelSet::DoubleCBOnly:    return "_dscb";
            case FitModelSet::DoubleGaussOnly: return "_dgauss";
            case FitModelSet::AllModels:       return "_allmodels";
            case FitModelSet::GaussianOnly:    default: return "";
        }
    }

    // Bukin: skewed Gaussian core with exponential tails, same parametrisation
    // as RooBukinPdf with rho1 = rho2 = 0. [0] peak height, [1] peak position,
    // [2] width, [3] asymmetry xi. At xi = 0 it is exactly a Gaussian.
    inline double BukinFunction(double* x, double* par) {
        const double Ap    = par[0];   // peak amplitude 
        const double Xp    = par[1];   // peak position
        const double sigma = par[2];   // width 
        const double xi    = par[3];   // asymmetry parameter
        if (sigma < 1e-9) return 0.0;

        const double consts = 2.0 * std::sqrt(2.0 * std::log(2.0));   // FWHM/sigma
        const double r3     = std::log(2.0);
        const double xx     = x[0];

        // Symmetric limit: pure Gaussian
        if (std::fabs(xi) < 1e-7) {
            const double r = (xx - Xp) / sigma;
            return Ap * std::exp(-0.5 * r * r);
        }

        const double r4 = std::sqrt(xi * xi + 1.0);
        const double r1 = xi / r4;
        const double r5 = xi / std::log(r4 + xi);                    // tail-slope coeff
        const double x1 = Xp + (sigma * consts) * (r1 - 1.0) / 2.0;  // core/tail boundaries
        const double x2 = Xp + (sigma * consts) * (r1 + 1.0) / 2.0;

        double f;
        if (xx < x1) {                       // left exponential tail
            f = std::exp(-r3 + 4.0 * r3 * (xx - x1) * r5 / (sigma * consts));
        } else if (xx > x2) {                // right exponential tail
            f = std::exp(-r3 - 4.0 * r3 * (xx - x2) * r5 / (sigma * consts));
        } else {                             // skewed-Gaussian core
            const double denom = std::log(1.0 + 2.0 * xi * (xi - r4));
            if (!std::isfinite(denom) || std::fabs(denom) < 1e-9) {
                const double r = (xx - Xp) / sigma;
                return Ap * std::exp(-0.5 * r * r);
            }
            const double inner = 1.0 + 4.0 * xi * r4 * (xx - Xp) / (sigma * consts);
            if (inner <= 0) return 0.0;
            const double r2 = std::log(inner) / denom;
            f = std::exp(-r3 * r2 * r2);
        }
        return std::isfinite(f) ? Ap * f : 0.0;
    }

    // Crystal Ball function (5 params): Gaussian core + power-law tail.
    // Power-law tail on one side (left for alpha > 0).
    // Parameters: [0]=amplitude, [1]=mean, [2]=sigma, [3]=alpha, [4]=n
    inline double CrystalBallFunction(double* x, double* par) {
        double amp    = par[0];
        double mean   = par[1];
        double sigma  = par[2];
        double alpha  = par[3];   
        double n      = par[4];   

        if (sigma < 1e-9 || n <= 1.0 || std::abs(alpha) < 1e-9) return 0.0;

        double t = (x[0] - mean) / sigma;
        if (alpha < 0) t = -t;   // Ensure tail is on the left side
        double absAlpha = std::abs(alpha);

        if (t >= -absAlpha) {
            // Gaussian core
            return amp * std::exp(-0.5 * t * t);
        } else {
            // Power-law tail
            double A = std::pow(n / absAlpha, n) * std::exp(-0.5 * absAlpha * absAlpha);
            double B = n / absAlpha - absAlpha;
            double denom = B - t;
            if (denom <= 0) return 0.0;
            return amp * A * std::pow(denom, -n);
        }
    }

    // Double-sided Crystal Ball (7 params): Gaussian core + INDEPENDENT power-law
    // tails on BOTH sides. (RooCrystalBall with one sigma and two tails). 
    //   [0]=amp [1]=mean [2]=sigma [3]=alphaL [4]=nL [5]=alphaR [6]=nR
    inline double DoubleCrystalBallFunction(double* x, double* par) {
        double amp   = par[0];
        double mean  = par[1];
        double sigma = par[2];
        double aL = std::abs(par[3]), nL = par[4];
        double aR = std::abs(par[5]), nR = par[6];
        if (sigma < 1e-9 || aL < 1e-9 || aR < 1e-9 || nL <= 1.0 || nR <= 1.0) return 0.0;

        double t = (x[0] - mean) / sigma;
        if (t < -aL) {                                   // left power-law tail
            double A = std::pow(nL / aL, nL) * std::exp(-0.5 * aL * aL);
            double B = nL / aL - aL;
            double d = B - t;
            return (d > 0) ? amp * A * std::pow(d, -nL) : 0.0;
        } else if (t > aR) {                             // right power-law tail
            double A = std::pow(nR / aR, nR) * std::exp(-0.5 * aR * aR);
            double B = nR / aR - aR;
            double d = B + t;
            return (d > 0) ? amp * A * std::pow(d, -nR) : 0.0;
        }
        return amp * std::exp(-0.5 * t * t);             // Gaussian core
    }

    // Double Gaussian (5 params): two concentric Gaussians sharing a mean — a
    // narrow core plus a wide one for the tails. 
    // Its effective RMS is a real standard deviation, so it can be used
    // in the JER quadrature subtraction, unlike the CB core sigma.
    inline double DoubleGaussianFunction(double* x, double* par) {
        double amp   = par[0];
        double mean  = par[1];
        double s1    = par[2];
        double frac  = par[3];
        double rwide = par[4];
        if (s1 < 1e-9 || rwide < 1.0) return 0.0;
        double s2 = s1 * rwide;
        double t1 = (x[0] - mean) / s1;
        double t2 = (x[0] - mean) / s2;
        return amp * (frac * std::exp(-0.5 * t1 * t1) +
                      (1.0 - frac) * std::exp(-0.5 * t2 * t2));
    }

    // RMS of the mixture. frac is a peak-height fraction; the area fraction
    // of each component is proportional to frac_i * sigma_i.
    inline double DoubleGaussianSigmaEff(double s1, double frac, double rwide) {
        if (s1 <= 0 || rwide < 1.0) return s1;
        double s2 = s1 * rwide;
        double a1 = frac        * s1;        
        double a2 = (1.0 - frac) * s2;       
        double W  = a1 + a2;
        if (W <= 0) return s1;
        double p1 = a1 / W, p2 = a2 / W;
        double var = p1 * s1 * s1 + p2 * s2 * s2;
        return (var > 0) ? std::sqrt(var) : s1;
    }

    // The four Fit*OverWindow functions fit one model over [lo, hi], seeded
    // from a Gaussian result. valid is false if the fit did not converge.
    // Try to fit a Bukin function over [lo, hi] with seed values from the
    // best Gaussian candidate so far. 
    inline GausFitResult FitBukinOverWindow(TH1* hist,
                                             const std::string& nameTag,
                                             double lo, double hi,
                                             double mu_seed, double sig_seed,
                                             double amp_seed) {
        GausFitResult r;
        r.modelName = "bukin";
        r.fitLo = lo; r.fitHi = hi;

        TF1 func(Form("bukin_%s", nameTag.c_str()), BukinFunction, lo, hi, 4);
        func.SetParameters(amp_seed, mu_seed, sig_seed, 0.0);  
        func.SetParLimits(2, sig_seed * 0.3, sig_seed * 3.0);  
        func.SetParLimits(3, -0.5, 0.5);                        

        int status = hist->Fit(&func, "RQ0");
        if (status != 0) return r;

        r.amp       = func.GetParameter(0);  r.ampErr   = func.GetParError(0);
        r.mean      = func.GetParameter(1);  r.meanErr  = func.GetParError(1);
        r.sigma     = func.GetParameter(2);  r.sigmaErr = func.GetParError(2);
        r.chi2      = func.GetChisquare();   r.ndf      = func.GetNDF();
        r.chi2ndf   = (r.ndf > 0) ? r.chi2 / r.ndf : 999;
        r.pars.clear();
        for (int i = 0; i < func.GetNpar(); ++i) r.pars.push_back(func.GetParameter(i));
        r.valid     = (r.ndf >= 1 && r.sigma > 0 && std::isfinite(r.chi2ndf));
        return r;
    }

    // Crystal Ball - tail on whichever side of the peak is wider.
    inline GausFitResult FitCrystalBallOverWindow(TH1* hist,
                                                    const std::string& nameTag,
                                                    double lo, double hi,
                                                    double mu_seed, double sig_seed,
                                                    double amp_seed,
                                                    double alpha_seed = 1.5) {
        GausFitResult r;
        r.modelName = "crystalball";
        r.fitLo = lo; r.fitHi = hi;

        TF1 func(Form("cb_%s", nameTag.c_str()), CrystalBallFunction, lo, hi, 5);
        func.SetParameters(amp_seed, mu_seed, sig_seed, alpha_seed, 3.0);
        func.SetParLimits(2, sig_seed * 0.3, sig_seed * 3.0);  
        func.SetParLimits(3, 0.5, 5.0);                         
        func.SetParLimits(4, 1.1, 20.0);                        

        int status = hist->Fit(&func, "RQ0");
        if (status != 0) return r;

        r.amp       = func.GetParameter(0);  r.ampErr   = func.GetParError(0);
        r.mean      = func.GetParameter(1);  r.meanErr  = func.GetParError(1);
        r.sigma     = func.GetParameter(2);  r.sigmaErr = func.GetParError(2);
        r.chi2      = func.GetChisquare();   r.ndf      = func.GetNDF();
        r.chi2ndf   = (r.ndf > 0) ? r.chi2 / r.ndf : 999;
        r.pars.clear();
        for (int i = 0; i < func.GetNpar(); ++i) r.pars.push_back(func.GetParameter(i));
        r.valid     = (r.ndf >= 1 && r.sigma > 0 && std::isfinite(r.chi2ndf));
        return r;
    }

    // Fit a double-sided Crystal Ball over [lo, hi], seeded from the best
    // Gaussian. Maps peak->mean, width->sigma (par[2]).
    inline GausFitResult FitDoubleCBOverWindow(TH1* hist,
                                               const std::string& nameTag,
                                               double lo, double hi,
                                               double mu_seed, double sig_seed,
                                               double amp_seed) {
        GausFitResult r;
        r.modelName = "dscb";
        r.fitLo = lo; r.fitHi = hi;

        TF1 func(Form("dscb_%s", nameTag.c_str()), DoubleCrystalBallFunction, lo, hi, 7);
        func.SetParameters(amp_seed, mu_seed, sig_seed, 1.5, 3.0, 1.5, 3.0);
        func.SetParLimits(2, sig_seed * 0.3, sig_seed * 3.0);  
        func.SetParLimits(3, 0.3, 5.0);                         
        func.SetParLimits(4, 1.1, 50.0);                        
        func.SetParLimits(5, 0.3, 5.0);                         
        func.SetParLimits(6, 1.1, 50.0);                        

        int status = hist->Fit(&func, "RQ0");
        if (status != 0) return r;

        r.amp       = func.GetParameter(0);  r.ampErr   = func.GetParError(0);
        r.mean      = func.GetParameter(1);  r.meanErr  = func.GetParError(1);
        r.sigma     = func.GetParameter(2);  r.sigmaErr = func.GetParError(2);
        r.chi2      = func.GetChisquare();   r.ndf      = func.GetNDF();
        r.chi2ndf   = (r.ndf > 0) ? r.chi2 / r.ndf : 999;
        r.pars.clear();
        for (int i = 0; i < func.GetNpar(); ++i) r.pars.push_back(func.GetParameter(i));
        r.valid     = (r.ndf >= 1 && r.sigma > 0 && std::isfinite(r.chi2ndf));
        return r;
    }

    // Fit a shared-mean double Gaussian over [lo, hi], seeded from the best
    // Gaussian. 
    inline GausFitResult FitDoubleGaussianOverWindow(TH1* hist,
                                                     const std::string& nameTag,
                                                     double lo, double hi,
                                                     double mu_seed, double sig_seed,
                                                     double amp_seed) {
        GausFitResult r;
        r.modelName = "doublegauss";
        r.fitLo = lo; r.fitHi = hi;

        TF1 func(Form("dgauss_%s", nameTag.c_str()), DoubleGaussianFunction, lo, hi, 5);
        func.SetParameters(amp_seed, mu_seed, sig_seed * 0.85, 0.7, 2.0);
        func.SetParLimits(2, sig_seed * 0.2, sig_seed * 1.5);  
        func.SetParLimits(3, 0.2, 1.0);                         
        func.SetParLimits(4, 1.0, 8.0);                         

        int status = hist->Fit(&func, "RQ0");
        if (status != 0) return r;

        double s1    = func.GetParameter(2);
        double frac  = func.GetParameter(3);
        double rwide = func.GetParameter(4);

        r.amp       = func.GetParameter(0);  r.ampErr  = func.GetParError(0);
        r.mean      = func.GetParameter(1);  r.meanErr = func.GetParError(1);
        r.sigma     = DoubleGaussianSigmaEff(s1, frac, rwide);
        {
            double ds1 = func.GetParError(2), df = func.GetParError(3), dr = func.GetParError(4);
            double base = r.sigma;
            double g1 = DoubleGaussianSigmaEff(s1 + ds1, frac, rwide) - base;
            double g2 = DoubleGaussianSigmaEff(s1, std::min(1.0, frac + df), rwide) - base;
            double g3 = DoubleGaussianSigmaEff(s1, frac, rwide + dr) - base;
            r.sigmaErr = std::sqrt(g1 * g1 + g2 * g2 + g3 * g3);
        }
        r.chi2      = func.GetChisquare();   r.ndf     = func.GetNDF();
        r.chi2ndf   = (r.ndf > 0) ? r.chi2 / r.ndf : 999;
        r.pars.clear();
        for (int i = 0; i < func.GetNpar(); ++i) r.pars.push_back(func.GetParameter(i));
        r.valid     = (r.ndf >= 1 && r.sigma > 0 && std::isfinite(r.chi2ndf));
        return r;
    }

    // Rebuild the ACTUAL fitted model as a TF1 over [r.fitLo, r.fitHi], so plots
    // draw the real curve (DSCB/CB/Bukin) rather than the Gaussian approximation.
    inline TF1* MakeModelTF1(const GausFitResult& r, const std::string& name) {
        const double lo = r.fitLo, hi = r.fitHi;
        if (r.modelName == "bukin" && r.pars.size() >= 4) {
            TF1* f = new TF1(name.c_str(), BukinFunction, lo, hi, 4);
            for (int i = 0; i < 4; ++i) f->SetParameter(i, r.pars[i]);
            return f;
        }
        if (r.modelName == "crystalball" && r.pars.size() >= 5) {
            TF1* f = new TF1(name.c_str(), CrystalBallFunction, lo, hi, 5);
            for (int i = 0; i < 5; ++i) f->SetParameter(i, r.pars[i]);
            return f;
        }
        if (r.modelName == "dscb" && r.pars.size() >= 7) {
            TF1* f = new TF1(name.c_str(), DoubleCrystalBallFunction, lo, hi, 7);
            for (int i = 0; i < 7; ++i) f->SetParameter(i, r.pars[i]);
            return f;
        }
        if (r.modelName == "doublegauss" && r.pars.size() >= 5) {
            TF1* f = new TF1(name.c_str(), DoubleGaussianFunction, lo, hi, 5);
            for (int i = 0; i < 5; ++i) f->SetParameter(i, r.pars[i]);
            return f;
        }
        TF1* f = new TF1(name.c_str(), "gaus", lo, hi);   // gaus fallback
        f->SetParameters(r.amp, r.mean, r.sigma);
        return f;
    }

    // Model name for plot legends.
    inline std::string ModelDisplayName(const std::string& m) {
        if (m == "bukin")       return "Bukin fit";
        if (m == "crystalball") return "Crystal Ball fit";
        if (m == "dscb")        return "Double-sided CB fit";
        if (m == "doublegauss") return "Double-Gaussian fit";
        return "Gaussian fit";
    }


    // =================================================================================
    // FitMode — selects which window-scan + selection algorithm
    // EASYBJETS_FIT_MODE controls the iterative fit's behavior.
    // =================================================================================
    //
    // --- fit window strategy (EASYBJETS_FIT_MODE) -------------------------------------
    //
    //   scan         (default) windows {+-1, +-1.5} sigma, or EASYBJETS_FIT_NSIGMA_LIST;
    //                Selects best candidate via the asymmetric penalty function
    //
    //   scan_dense   Dense symmetric and asymmetric windows up to 3.5 sigma; selects best candidate by
    //                by minimum |log(chi2/ndf)| (log-symmetric distance from 1) among those passing the 
    //                the content cut and the ndf/chi2 floor
    //                
    //   scan_log     Same selection rule as scan_dense, but windows from 1.5 sigma up; selects best candidate by
    //                by minimum |log(chi2/ndf)| (log-symmetric distance from 1) among those passing the 
    //                the floor, no content cut
    //
    //   legacy       Single fixed window of +-EASYBJETS_FIT_NSIGMA (default 2.0) sigma 
    //
    //   legacy7grid  Scan over a fixed N-candidate grid:, best by |log(chi2/ndf)| with
    //                the ndf/chi2 floor
    //
    // In every mode sigma is limited to [0, 1] in the fits and a fit whose
    // mean lands outside its own window is discarded.
    //
    // Examples:
    //   export EASYBJETS_FIT_MODE=scan         # current default (if unset)
    //   export EASYBJETS_FIT_MODE=scan_dense   # dense grid (sym+asym) + log-distance + ndf floor
    //   export EASYBJETS_FIT_MODE=scan_log     # 9-point grid + log-distance + ndf floor
    //   export EASYBJETS_FIT_MODE=legacy       # +/-2sigma single window
    //   export EASYBJETS_FIT_MODE=legacy7grid  # 7-grid + log-distance + ndf floor
    //   export EASYBJETS_FIT_NSIGMA=2.5        # legacy mode: use +/-2.5sigma instead
    //
    enum class FitMode { Scan, ScanDense, ScanLog, Legacy, Legacy7Grid };

    inline FitMode CurrentFitMode() {
        static const FitMode m = []() {
            const char* e = std::getenv("EASYBJETS_FIT_MODE");
            if (!e || !*e) return FitMode::Scan;  
            std::string s(e);
            for (char& c : s) c = std::tolower((unsigned char)c);
            while (!s.empty() && std::isspace((unsigned char)s.front())) s.erase(s.begin());
            while (!s.empty() && std::isspace((unsigned char)s.back()))  s.pop_back();
            if (s == "scan")        return FitMode::Scan;
            if (s == "scan_dense" || s == "scandense" || s == "dense")
                                    return FitMode::ScanDense;
            if (s == "scan_log" || s == "scanlog" || s == "log")
                                    return FitMode::ScanLog;
            if (s == "legacy")      return FitMode::Legacy;
            if (s == "legacy7grid" || s == "legacy7" || s == "7grid")
                                    return FitMode::Legacy7Grid;
            std::cerr << "[FitMode] WARN: EASYBJETS_FIT_MODE='" << e
                      << "' not recognized. Valid: scan|scan_dense|scan_log|legacy|legacy7grid. "
                      << "Falling back to 'scan'.\n";
            return FitMode::Scan;
        }();
        return m;
    }

    inline const char* FitModeName(FitMode m) {
        switch (m) {
            case FitMode::Scan:        return "scan";
            case FitMode::ScanDense:   return "scan_dense";
            case FitMode::ScanLog:     return "scan_log";
            case FitMode::Legacy:      return "legacy";
            case FitMode::Legacy7Grid: return "legacy7grid";
        }
        return "unknown";
    }

    // For legacy mode: which fixed N-sigma window to use.
    // EASYBJETS_FIT_NSIGMA scalar (default 2.0).
    inline double CurrentLegacyNSigma() {
        static const double n = []() {
            const char* e = std::getenv("EASYBJETS_FIT_NSIGMA");
            if (!e || !*e) return 2.0;
            try { double v = std::stod(e); if (v > 0 && v < 10) return v; }
            catch (...) {}
            return 2.0;
        }();
        return n;
    }

    // ===========================================================
    // SeedMode - Seed window of the first fit (EASYBJETS_SEED_MODE):
    //
    //   truncated  (default) - caller provides seedLo/seedHi (typically
    //              [0.8, 1.2] for response). The seed Gaussian is fit over
    //              this fixed window. Stable and robust for symmetric peaks.
    //
    //   raw        mean +- maxNSigma * RMS of the histogram. Lets the
    //              raw histogram statistics drive the seed; better for
    //              asymmetric peaks because no truncation bias is introduced.
    //              Caller-provided seedLo/seedHi are IGNORED.
    // ===========================================================
    enum class SeedMode { Truncated, Raw };

    inline SeedMode CurrentSeedMode() {
        static const SeedMode m = []() {
            const char* e = std::getenv("EASYBJETS_SEED_MODE");
            if (!e || !*e) return SeedMode::Truncated;
            std::string s(e);
            for (char& c : s) c = std::tolower((unsigned char)c);
            if (s == "raw" || s == "atlas" || s == "reference") return SeedMode::Raw;
            if (s == "truncated" || s == "default") return SeedMode::Truncated;
            std::cerr << "[SeedMode] WARN: unknown EASYBJETS_SEED_MODE='" << e
                      << "', using 'truncated'.\n";
            return SeedMode::Truncated;
        }();
        return m;
    }

    inline const char* SeedModeName(SeedMode m) {
        return (m == SeedMode::Raw) ? "raw" : "truncated";
    }

    // Per-pT-bin physical floor for response fits:
    //   Rmin = fit_minPt / pT_bin_lo
    // where fit_minPt is set by EASYBJETS_RESPONSE_FIT_PTMIN (default unset).
    // When the env var is set and the caller passes a positive pT_bin_lo,
    // the per-pT floor is computed and added to the existing physicalLo via
    // std::max - so a pT-bin floor can only raise the floor, never lower it.
    //
    inline double CurrentFitMinPt() {
        static const double v = []() {
            const char* e = std::getenv("EASYBJETS_RESPONSE_FIT_PTMIN");
            if (!e || !*e) return -1.0;
            try { double f = std::stod(e); if (f > 0 && f < 1000) return f; }
            catch (...) {}
            return -1.0;
        }();
        return v;
    }

    inline GausFitResult IterativeGausFit(
            TH1* hist,
            const std::string& nameTag,
            const std::string& label,
            const std::string& scanPdfPath = "",
            double seedLo    = 0.8, 
            double seedHi    = 1.2, 
            double maxNSigma = 2., 
            double physicalLo = -1.0,
            double pTBinLo    = -1.0)
    {
        struct NSigPair { double lo; double hi; };
        std::vector<NSigPair> nSigmaValues;
        {
            std::vector<NSigPair> defaultList = {
                {1.0, 1.0}, {1.5, 1.5}
            };
            const char* envList = std::getenv("EASYBJETS_FIT_NSIGMA_LIST");
            if (envList && *envList) {
                defaultList.clear();
                std::string s(envList);
                std::vector<std::string> toks;
                {
                    std::string cur;
                    for (char c : s) {
                        if (c == ',') { if (!cur.empty()) toks.push_back(cur); cur.clear(); }
                        else if (!std::isspace((unsigned char)c)) cur.push_back(c);
                    }
                    if (!cur.empty()) toks.push_back(cur);
                }
                for (const auto& t : toks) {
                    auto colonPos = t.find(':');
                    try {
                        if (colonPos == std::string::npos) {
                            double v = std::stod(t);
                            if (v > 0) defaultList.push_back({v, v});
                        } else {
                            double lo = std::stod(t.substr(0, colonPos));
                            double hi = std::stod(t.substr(colonPos + 1));
                            if (lo > 0 && hi > 0) defaultList.push_back({lo, hi});
                        }
                    } catch (...) {}
                }
                if (defaultList.empty()) defaultList = {{1.0,1.0}, {1.5,1.5}, {2.0,2.0}};
            }
            for (const auto& p : defaultList)
                if (p.lo <= maxNSigma + 1e-9 && p.hi <= maxNSigma + 1e-9)
                    nSigmaValues.push_back(p);
        }

        // Mode override: The nSigmaValues built above is the user-supplied or 
        // default grid for 'scan' mode. For 'legacy' and 'legacy7grid' modes, override here.
        // The mode is read once per process (cached static in CurrentFitMode).
        const FitMode fitMode = CurrentFitMode();
        if (fitMode == FitMode::Legacy) {
            // Single fixed window - N-sigma scalar via EASYBJETS_FIT_NSIGMA.
            double N = CurrentLegacyNSigma();
            // maxNSigma cap from the call site for safety.
            if (N > maxNSigma) N = maxNSigma;
            nSigmaValues.clear();
            nSigmaValues.push_back({N, N});
        } else if (fitMode == FitMode::Legacy7Grid) {
            // Fixed N-candidate grid (asymmetric on the low side to capture ISR).
            std::vector<NSigPair> grid = {
                //{1.0, 1.0}, {1.5, 1.0}, {1.5, 1.5}, {2.0, 1.5},
                //{2.0, 2.0}, {2.5, 1.5}, {3.0, 2.0}
                {1.0, 1.0}, {1.25, 1.0}, {1.5, 1.0}, {1.5, 1.25},
                {1.5, 1.5}, {2.0, 1.5},
                {2.0, 2.0}, {2.5, 1.5}, {3.0, 2.0}
            };
            nSigmaValues.clear();
            for (const auto& p : grid) {
                if (p.lo <= maxNSigma + 1e-9 && p.hi <= maxNSigma + 1e-9)
                    nSigmaValues.push_back(p);
            }
        } else if (fitMode == FitMode::ScanDense) {
            // Dense grid: the 9 symmetric windows
            // PLUS asymmetric windows (wider on the low-R side for the ISR tail).
            // Selected by |log(Chi^2/ndf)| (log-symmetric distance from 1), not the
            // asymmetric penalty, with the under-fit floor guarding narrow
            // candidates the dense grid would otherwise expose more often.
            std::vector<NSigPair> grid = {
                // Symmetric
                {1.00, 1.00}, {1.25, 1.25}, {1.50, 1.50},
                {1.75, 1.75}, {2.00, 2.00}, {2.25, 2.25},
                {2.50, 2.50}, {2.75, 2.75}, {3.00, 3.00},
                // Asymmetric - wider on low-R (down-sigma) for the ISR tail
                {1.50, 1.00}, {1.75, 1.50}, {2.00, 1.50}, {2.00, 1.75},
                {2.50, 1.75}, {2.50, 2.00}, {3.00, 2.00}, {3.00, 2.50},
                {3.50, 2.50}
            };
            nSigmaValues.clear();
            for (const auto& p : grid) {
                if (p.lo <= maxNSigma + 1e-9 && p.hi <= maxNSigma + 1e-9)
                    nSigmaValues.push_back(p);
            }
        } else if (fitMode == FitMode::ScanLog) {
            // Same N-point grid as ScanDense; differs only in selection rule
            std::vector<NSigPair> grid = {
                // Symmetric 
                {1.50, 1.50},
                {1.75, 1.75}, {2.00, 2.00}, {2.25, 2.25},
                {2.50, 2.50}, {2.75, 2.75}, {3.00, 3.00},
                // Asymmetric - wider on low-R for ISR tail
                {1.75, 1.50}, {2.00, 1.50},
                {2.00, 1.75}, {2.50, 1.75},
                {2.50, 2.00}, {3.00, 2.00}, {3.00, 2.50},
                {3.50, 2.50}
            };
            nSigmaValues.clear();
            for (const auto& p : grid) {
                if (p.lo <= maxNSigma + 1e-9 && p.hi <= maxNSigma + 1e-9)
                    nSigmaValues.push_back(p);
            }
        }
        // FitMode::Scan: keep nSigmaValues as built above.
        const int    maxIter  = 15;
        const double convTol  = 1e-4;
        const double xMin     = hist->GetXaxis()->GetXmin();
        const double xMax     = hist->GetXaxis()->GetXmax();
        const double binW     = hist->GetBinWidth(1);

        // Effective low-edge floor: combine the user-supplied physical floor
        // (physicalLo) with the per-pT-bin floor Rmin = fit_minPt / pTBinLo
        // (when both are positive and EASYBJETS_RESPONSE_FIT_PTMIN is set).
        // Take the MAX of all three (xMin, physicalLo, Rmin) so we never go
        // below the most restrictive bound.
        double combinedFloor = xMin;
        if (physicalLo >= 0.0) combinedFloor = std::max(combinedFloor, physicalLo);
        const double fit_minPt = CurrentFitMinPt();
        if (fit_minPt > 0 && pTBinLo > 0) {
            const double rmin = fit_minPt / pTBinLo;
            combinedFloor = std::max(combinedFloor, rmin);
        }
        const double effLoFloor = combinedFloor;
        std::vector<GausFitResult> candidates;

        // Seed-window override - when EASYBJETS_SEED_MODE=raw, override the caller 
        // seedLo/seedHi with values computed from the raw histogram stats:
        //   seedLo = h.GetMean() - maxNSigma * h.GetStdDev()
        //   seedHi = h.GetMean() + maxNSigma * h.GetStdDev()
        // No truncation, lets the histogram drive the seed. 
        // Better for asymmetric distributions.
        if (CurrentSeedMode() == SeedMode::Raw) {
            const double rawMean = hist->GetMean();
            const double rawRMS  = hist->GetStdDev();
            if (rawRMS > 1e-9) {
                seedLo = rawMean - maxNSigma * rawRMS;
                seedHi = rawMean + maxNSigma * rawRMS;
                if (seedLo < xMin) seedLo = xMin;
                if (seedHi > xMax) seedHi = xMax;
            }
        }

        // --- Step 1: Seed fit ---
        // Clamp the seed window to the physical floor too.
        double seedLoClamped = std::max(effLoFloor, seedLo);
        TF1* seedFunc = new TF1(("seed_" + nameTag).c_str(), "gaus", seedLoClamped, seedHi);
        seedFunc->SetParameters(hist->GetMaximum(), hist->GetMean(), hist->GetStdDev());
        seedFunc->SetParLimits(2, 0.0, 1.0);
        hist->Fit(seedFunc, "RQ0");
        double mu0  = seedFunc->GetParameter(1);
        double sig0 = seedFunc->GetParameter(2);

        GausFitResult seed;
        seed.amp      = seedFunc->GetParameter(0);  seed.ampErr   = seedFunc->GetParError(0);
        seed.mean     = mu0;                         seed.meanErr  = seedFunc->GetParError(1);
        seed.sigma    = sig0;                        seed.sigmaErr = seedFunc->GetParError(2);
        seed.chi2     = seedFunc->GetChisquare();    seed.ndf      = seedFunc->GetNDF();
        seed.chi2ndf  = (seed.ndf > 0) ? seed.chi2 / seed.ndf : 999;
        seed.fitLo    = seedLoClamped;               seed.fitHi    = seedHi;
        seed.nSigma   = 0;
        seed.nSigmaLo = 0;
        seed.nSigmaHi = 0;
        seed.modelName = "gaus";
        seed.pars      = {seed.amp, seed.mean, seed.sigma};
        seed.valid    = (seed.ndf >= 1 && sig0 > 0 && std::isfinite(seed.chi2ndf));
        candidates.push_back(seed);
        delete seedFunc;

        // --- Step 2: Iterative window scan ---
        // Each candidate runs an inner loop that updates the window as (mu, sigma)
        // drift. 
        if (seed.valid && sig0 > 1e-6) {
            for (const auto& NP : nSigmaValues) {
                const double N_lo = NP.lo;
                const double N_hi = NP.hi;
                double mu  = mu0;
                double sig = sig0;
                double amp = seed.amp;
                GausFitResult converged;
                converged.valid = false;

                // Store last few iterates so we can
                // detect sigma flipping between two values (without converging).
                std::vector<GausFitResult> history;
                history.reserve(maxIter);
                int oscillationCount = 0;

                for (int iter = 0; iter < maxIter; ++iter) {
                    double lo = mu - N_lo * sig;
                    double hi = mu + N_hi * sig;
                    if (lo < effLoFloor) lo = effLoFloor;
                    if (hi > xMax)       hi = xMax;
                    if (hi - lo < 3 * binW) break;

                    // tag: encode both nSig values for unique TF1 names
                    TF1* func = new TF1(Form("iter_%s_%.0fL_%.0fH_i%d",
                                             nameTag.c_str(),
                                             N_lo * 10, N_hi * 10, iter),
                                        "gaus", lo, hi);
                    func->SetParameters(amp, mu, sig);
                    func->SetParLimits(2, 0.0, 1.0);
                    hist->Fit(func, "RQ0");

                    double mu_new  = func->GetParameter(1);
                    double sig_new = func->GetParameter(2);

                    if (mu_new < lo || mu_new > hi) {
                        if (history.empty()) {
                            // First iteration already failed the guard.
                            // Mark this candidate invalid and move on.
                            converged.valid = false;
                            delete func;
                            break;
                        }
                        // Restore the most recent successful iterate.
                        converged = history.back();
                        std::cerr << "[FitMode] WARN: mass guard triggered for '"
                                  << nameTag << "' iter " << iter
                                  << " (mu_fit=" << mu_new
                                  << " outside [" << lo << ", " << hi << "])"
                                  << "; restoring previous iterate.\n";
                        delete func;
                        break;
                    }

                    converged.amp      = func->GetParameter(0);  converged.ampErr   = func->GetParError(0);
                    converged.mean     = mu_new;                  converged.meanErr  = func->GetParError(1);
                    converged.sigma    = sig_new;                 converged.sigmaErr = func->GetParError(2);
                    converged.chi2     = func->GetChisquare();    converged.ndf      = func->GetNDF();
                    converged.chi2ndf  = (converged.ndf > 0) ? converged.chi2 / converged.ndf : 999;
                    converged.fitLo    = lo;                      converged.fitHi    = hi;
                    // set nSigma to lo when lo == hi else the average.
                    converged.nSigma   = (std::abs(N_lo - N_hi) < 1e-9) ? N_lo : 0.5 * (N_lo + N_hi);
                    converged.nSigmaLo = N_lo;
                    converged.nSigmaHi = N_hi;
                    converged.modelName = "gaus";
                    converged.pars      = {converged.amp, converged.mean, converged.sigma};
                    converged.valid    = (converged.ndf >= 1 && sig_new > 0 &&
                                         std::isfinite(converged.chi2ndf));
                    delete func;

                    if (!converged.valid) break;
                    history.push_back(converged);

                    bool muConv  = (std::abs(mu)  > 1e-12) ? std::abs((mu_new  - mu)  / mu)  < convTol : std::abs(mu_new  - mu)  < convTol;
                    bool sigConv = (std::abs(sig) > 1e-12) ? std::abs((sig_new - sig) / sig) < convTol : std::abs(sig_new - sig) < convTol;

                    // Sigma oscillating > 5% between consecutive
                    // iterations for 3+ iter in a row means the fit has
                    // hit a flip-flop between two local minima.
                    if (iter >= 5 && history.size() >= 2) {
                        const auto& prev = history[history.size()-2];
                        double relDelta = (prev.sigma > 1e-12)
                            ? std::abs((sig_new - prev.sigma) / prev.sigma)
                            : std::abs(sig_new - prev.sigma);
                        if (relDelta > 0.05) ++oscillationCount;
                        else                 oscillationCount = 0;
                    }

                    amp = converged.amp;
                    mu  = mu_new;
                    sig = sig_new;

                    if (muConv && sigConv) break;

                    if (oscillationCount >= 3) {
                        // Oscillating. Pick the iterate with lowest Chi^2/ndf
                        // (i.e., closest to the "right" Gaussian shape)
                        // from the history - this is better than just
                        // returning the last iterate, which is arbitrary.
                        double bestChi2ndf = 1e9;
                        for (const auto& h : history) {
                            if (h.valid && h.chi2ndf < bestChi2ndf) {
                                bestChi2ndf = h.chi2ndf;
                                converged = h;
                            }
                        }
                        break;
                    }
                }

                if (converged.valid) candidates.push_back(converged);
            }
        }

        // --- Step 2b: Alternative-model candidates ---
        // FitModelSet=AllModels: add Bukin + CrystalBall as extra candidates
        //   alongside the Gaussian grid. Selector picks the best overall.
        // FitModelSet=BukinOnly or CrystalBallOnly: fit only the requested
        //   model over the seed window, then clear the Gaussian candidates
        //   so only the alternative remains in the pool. (The seed Gaussian
        //   is still needed internally as starting values for the alt fit.)
        FitModelSet fm = CurrentFitModelSet();
        if ((fm == FitModelSet::AllModels ||
             fm == FitModelSet::BukinOnly ||
             fm == FitModelSet::CrystalBallOnly ||
             fm == FitModelSet::DoubleCBOnly ||
             fm == FitModelSet::DoubleGaussOnly) && !candidates.empty()) {

            // Find the best-by-penalty Gaussian so far, use it as seed.
            // For BukinOnly/CrystalBallOnly this falls back to the seed fit
            // from Step 1 if no N-sigma iterates survived.
            const GausFitResult* best_so_far = nullptr;
            double bestPen = 1e9;
            for (const auto& c : candidates) {
                if (!c.valid) continue;
                double pen = c.Penalty();
                if (pen < bestPen) { bestPen = pen; best_so_far = &c; }
            }
            if (best_so_far) {
                const bool want_bukin = (fm == FitModelSet::AllModels || fm == FitModelSet::BukinOnly);
                const bool want_cb    = (fm == FitModelSet::AllModels || fm == FitModelSet::CrystalBallOnly);
                const bool want_dscb  = (fm == FitModelSet::AllModels || fm == FitModelSet::DoubleCBOnly);
                const bool want_dgauss= (fm == FitModelSet::AllModels || fm == FitModelSet::DoubleGaussOnly);

                // Crucially we SCAN the alt model over EVERY Gaussian candidate
                // window (not just the single best), so scan_log/scan can select the
                // best alt-model window the same way it does for the Gaussian — and
                // so a fit that fails at one window still has others to succeed at.
                double loN = 3.0, hiN = 2.0;
                if (const char* e = std::getenv("EASYBJETS_FIT_MODEL_NSIGMA_LO")) { try { double v=std::stod(e); if(v>=1.0&&v<=8.0) loN=v; } catch(...){} }
                if (const char* e = std::getenv("EASYBJETS_FIT_MODEL_NSIGMA_HI")) { try { double v=std::stod(e); if(v>=1.0&&v<=8.0) hiN=v; } catch(...){} }
                const double axLo = hist->GetXaxis()->GetXmin();
                const double axHi = hist->GetXaxis()->GetXmax();

                std::vector<GausFitResult> gaussCands;
                for (const auto& c : candidates)
                    if (c.valid && c.modelName == "gaus") gaussCands.push_back(c);

                int wIdx = 0;
                for (const auto& g : gaussCands) {
                    double altLo = g.mean - loN * g.sigma;
                    double altHi = g.mean + hiN * g.sigma;
                    if (altLo < axLo) altLo = axLo;
                    if (altHi > axHi) altHi = axHi;
                    const std::string wtag = Form("%s_aw%d", nameTag.c_str(), wIdx++);

                    if (want_bukin) {
                        GausFitResult bukin = FitBukinOverWindow(
                            hist, wtag, altLo, altHi, g.mean, g.sigma, g.amp);
                        bukin.nSigma = g.nSigma; bukin.nSigmaLo = g.nSigmaLo; bukin.nSigmaHi = g.nSigmaHi;
                        if (bukin.valid) candidates.push_back(bukin);
                    }
                    if (want_cb) {
                        GausFitResult cb = FitCrystalBallOverWindow(
                            hist, wtag, altLo, altHi, g.mean, g.sigma, g.amp);
                        cb.nSigma = g.nSigma; cb.nSigmaLo = g.nSigmaLo; cb.nSigmaHi = g.nSigmaHi;
                        if (cb.valid) candidates.push_back(cb);
                    }
                    if (want_dscb) {
                        GausFitResult dscb = FitDoubleCBOverWindow(
                            hist, wtag, altLo, altHi, g.mean, g.sigma, g.amp);
                        dscb.nSigma = g.nSigma; dscb.nSigmaLo = g.nSigmaLo; dscb.nSigmaHi = g.nSigmaHi;
                        if (dscb.valid) candidates.push_back(dscb);
                    }
                    if (want_dgauss) {
                        GausFitResult dg = FitDoubleGaussianOverWindow(
                            hist, wtag, altLo, altHi, g.mean, g.sigma, g.amp);
                        dg.nSigma = g.nSigma; dg.nSigmaLo = g.nSigmaLo; dg.nSigmaHi = g.nSigmaHi;
                        if (dg.valid) candidates.push_back(dg);
                    }
                }

                // For the selectedOnly modes: strip out all Gaussian candidates so
                // only the alternative model is available for selection. We
                // keep only the candidates whose modelName matches the request.
                if (fm == FitModelSet::BukinOnly || fm == FitModelSet::CrystalBallOnly ||
                    fm == FitModelSet::DoubleCBOnly || fm == FitModelSet::DoubleGaussOnly) {
                    const std::string keep = (fm == FitModelSet::BukinOnly)       ? "bukin"
                                           : (fm == FitModelSet::CrystalBallOnly) ? "crystalball"
                                           : (fm == FitModelSet::DoubleCBOnly)    ? "dscb"
                                                                                  : "doublegauss";
                    std::vector<GausFitResult> filtered;
                    for (const auto& c : candidates)
                        if (c.modelName == keep) filtered.push_back(c);
                    if (!filtered.empty()) {
                        candidates = std::move(filtered);
                    } else {
                        // If the requested fit failed, fall back to the Gaussian scan
                        std::cerr << "[FitModel] WARN: requested '" << keep
                                  << "' fit failed on '" << nameTag
                                  << "'; falling back to the Gaussian scan\n";
                    }
                }
            }
        }

        // --- Step 3: Select best ---
        // Selection is FitMode mode-aware
        //            
        const double kMaxChi2NdfForSelection = 5.0;
        const double kMinFilledFraction      = 0.70;
        const int    kMinNDFLegacy7          = 2; 
        const double kMinChi2NdfLegacy7      = 0.01; 

        auto isSelectable = [&](const GausFitResult& c) {
            return c.valid
                && c.chi2ndf < kMaxChi2NdfForSelection
                && HasEnoughContentInWindow(hist, c.fitLo, c.fitHi, kMinFilledFraction);
        };

        int bestIdx = -1;
        double bestPenalty = 1e9;

        if (fitMode == FitMode::Legacy) {
            // The candidate pool is { seed, single +/-Nsigma iterate }. Always select
            // the iterate (index 1). If it failed (only seed present), fall back
            // to the seed and warn.
            if (candidates.size() >= 2 && candidates[1].valid) {
                bestIdx = 1;
            } else if (candidates.size() >= 1 && candidates[0].valid) {
                std::cerr << "[FitMode] WARN: legacy +/-" << CurrentLegacyNSigma()
                          << "sigma iterate failed for '" << nameTag
                          << "'; falling back to seed.\n";
                bestIdx = 0;
            }
        } else if (fitMode == FitMode::Legacy7Grid) {
            // Step 3a: among the 7-grid candidates that pass the floor,
            // pick the one with smallest |log(Chi^2/ndf)|.
            double bestLogDist = 1e9;
            for (size_t i = 0; i < candidates.size(); ++i) {
                const auto& c = candidates[i];
                if (!c.valid) continue;
                if (c.modelName != "gaus") continue;        
                if (i == 0) continue;                       
                if (c.ndf < kMinNDFLegacy7) continue;
                if (c.chi2ndf < kMinChi2NdfLegacy7) continue;
                double dist = (c.chi2ndf > 0) ? std::abs(std::log(c.chi2ndf)) : 1e9;
                if (dist < bestLogDist) { bestLogDist = dist; bestIdx = (int)i; }
            }
            // Step 3b: fallback - no candidate passed the floor. Pick the candidate
            // whose (ndf, Chi^2/ndf) is closest to the floor from above and warn.
            if (bestIdx < 0) {
                std::cerr << "[FitMode] WARN: no legacy7grid candidate passed the "
                          << "ndf>=" << kMinNDFLegacy7
                          << " AND Chi^2/ndf>=" << kMinChi2NdfLegacy7
                          << " floor for '" << nameTag
                          << "'. Selecting candidate closest to floor.\n";
                double bestFloorDev = 1e9;
                for (size_t i = 0; i < candidates.size(); ++i) {
                    const auto& c = candidates[i];
                    if (!c.valid) continue;
                    if (c.modelName != "gaus") continue;
                    if (i == 0) continue;
                    // Distance from floor: how far below ndf=6 OR Chi^2/ndf=0.3
                    double ndfShort = std::max(0.0, kMinNDFLegacy7 - (double)c.ndf);
                    double chi2Short = std::max(0.0, kMinChi2NdfLegacy7 - c.chi2ndf);
                    double dev = ndfShort + 10.0 * chi2Short;  
                    if (dev < bestFloorDev) { bestFloorDev = dev; bestIdx = (int)i; }
                }
            }
        } else if (fitMode == FitMode::ScanDense) {
            // Log-distance selection: among the dense grid candidates (symmetric
            // + asymmetric) passing the under-fit floor (ndf >= kMinNDFLegacy7 AND
            // Chi^2/ndf >= kMinChi2NdfLegacy7) and the window-content cut, pick the
            // one whose Chi^2/ndf is log-symmetrically closest to 1.
            // |log(Chi^2/ndf)| treats 0.5 and 2.0 as equally far from 1 (both 0.69),
            // unlike the asymmetric Penalty() this mode used previously.
            double bestLogDist = 1e9;
            for (size_t i = 0; i < candidates.size(); ++i) {
                const auto& c = candidates[i];
                if (!isSelectable(c)) continue;
                if (c.modelName != "gaus") continue;        
                if (i == 0) continue;                        
                if (c.ndf < kMinNDFLegacy7) continue;        
                if (c.chi2ndf < kMinChi2NdfLegacy7) continue;
                double dist = (c.chi2ndf > 0) ? std::abs(std::log(c.chi2ndf)) : 1e9;
                if (dist < bestLogDist) { bestLogDist = dist; bestIdx = (int)i; }
            }
            // Fallback: if no candidate passed the floor, drop the floor and
            // still use log-distance so the selection rule stays consistent.
            if (bestIdx < 0) {
                std::cerr << "[FitMode] WARN: no scan_dense candidate passed the "
                          << "ndf>=" << kMinNDFLegacy7
                          << " AND Chi^2/ndf>=" << kMinChi2NdfLegacy7
                          << " floor for '" << nameTag
                          << "'. Falling back to log-distance with no floor.\n";
                bestLogDist = 1e9;
                for (size_t i = 0; i < candidates.size(); ++i) {
                    const auto& c = candidates[i];
                    if (!c.valid) continue;
                    if (c.modelName != "gaus") continue;
                    if (i == 0) continue;
                    double dist = (c.chi2ndf > 0) ? std::abs(std::log(c.chi2ndf)) : 1e9;
                    if (dist < bestLogDist) { bestLogDist = dist; bestIdx = (int)i; }
                }
            }
        } else if (fitMode == FitMode::ScanLog) {
            // Log-distance selection: pick the candidate whose Chi^2/ndf is
            // log-symmetrically closest to 1, among those passing the
            // under-fit floor (ndf >= 6 AND Chi^2/ndf >= 0.3).
            // |log(Chi^2/ndf)| treats 0.5 and 2.0 as equally far from 1 (both = 0.69),
            // unlike the linear |Chi^2/ndf - 1| which would call 0.5 closer than 2.0.
            double bestLogDist = 1e9;
            for (size_t i = 0; i < candidates.size(); ++i) {
                const auto& c = candidates[i];
                if (!c.valid) continue;
                if (c.modelName != "gaus") continue;        
                if (i == 0) continue;                        
                if (c.ndf < kMinNDFLegacy7) continue;
                if (c.chi2ndf < kMinChi2NdfLegacy7) continue;
                double dist = (c.chi2ndf > 0) ? std::abs(std::log(c.chi2ndf)) : 1e9;
                if (dist < bestLogDist) { bestLogDist = dist; bestIdx = (int)i; }
            }
            // Fallback: if no candidate passes the floor, drop the floor and
            // still use log-distance (rather than penalty) so the rule remains
            // consistent. Warn so the user knows the floor wasn't satisfied.
            if (bestIdx < 0) {
                std::cerr << "[FitMode] WARN: no scan_log candidate passed the "
                          << "ndf>=" << kMinNDFLegacy7
                          << " AND Chi^2/ndf>=" << kMinChi2NdfLegacy7
                          << " floor for '" << nameTag
                          << "'. Falling back to log-distance with no floor.\n";
                bestLogDist = 1e9;
                for (size_t i = 0; i < candidates.size(); ++i) {
                    const auto& c = candidates[i];
                    if (!c.valid) continue;
                    if (c.modelName != "gaus") continue;
                    if (i == 0) continue;
                    double dist = (c.chi2ndf > 0) ? std::abs(std::log(c.chi2ndf)) : 1e9;
                    if (dist < bestLogDist) { bestLogDist = dist; bestIdx = (int)i; }
                }
            }
        } else {
            // FitMode::Scan: pick the candidate with the lowest penalty.
            for (size_t i = 0; i < candidates.size(); ++i) {
                if (!isSelectable(candidates[i])) continue;
                double pen = candidates[i].Penalty();
                if (pen < bestPenalty) { bestPenalty = pen; bestIdx = (int)i; }
            }
        }

        // Scan-mode fallbacks (the legacy/legacy7grid modes have their own
        // fallback paths above; do not mix). These activate only when we're in
        // FitMode::Scan and the strict selection above failed.
        if (fitMode == FitMode::Scan) {
            // Fallback 1: the window-content cut rejected everything -> drop just
            // that cut but still apply the Chi^2/ndf cap and the asymmetric penalty.
            if (bestIdx < 0) {
                for (size_t i = 0; i < candidates.size(); ++i) {
                    const auto& c = candidates[i];
                    if (!c.valid || c.chi2ndf >= kMaxChi2NdfForSelection) continue;
                    double pen = c.Penalty();
                    if (pen < bestPenalty) { bestPenalty = pen; bestIdx = (int)i; }
                }
            }

            // Fallback 2: even the Chi^2/ndf cap rejected everything -> revert to the
            // original "closest to 1" rule so we always return something non-empty.
            if (bestIdx < 0) {
                bestIdx = 0;
                double bestDev = candidates[0].valid ? candidates[0].AbsDevFrom1() : 999;
                for (size_t i = 1; i < candidates.size(); ++i) {
                    if (!candidates[i].valid) continue;
                    double dev = candidates[i].AbsDevFrom1();
                    if (dev < bestDev) { bestDev = dev; bestIdx = (int)i; }
                }
            }
        }

        // Final hard fallback for ANY mode: if still nothing valid, return seed.
        if (bestIdx < 0) {
            std::cerr << "[FitMode] WARN: no valid candidate found in mode '"
                      << FitModeName(fitMode) << "' for '" << nameTag
                      << "'; returning seed.\n";
            bestIdx = 0;
        }

        // --- Step 4: Print table ---
        std::cout << "\n  +---------------------- Iterative Gauss Fit: " << label
                  << " [mode=" << FitModeName(fitMode)
                  << ", seed=" << SeedModeName(CurrentSeedMode());
        if (effLoFloor > xMin + 1e-9) {
            std::cout << ", floor=" << Form("%.3f", effLoFloor);
        }
        std::cout << "] ------------------------------+" << std::endl;
        std::cout << "    │ Window           │    mu +/- dmu       │    sigma +/- dsigma       │ Chi^2/ndf        │ [lo, hi]          │" << std::endl;
        std::cout << "    +---------------------+---------------------+---------------------+------------------------+" << std::endl;
        for (size_t i = 0; i < candidates.size(); ++i) {
            const auto& c = candidates[i];
            std::string tag;
            if (i == 0) tag = Form("seed [%.2f,%.2f]", seedLo, seedHi);
            else if (c.modelName == "bukin")       tag = Form("bukin (asym)");
            else if (c.modelName == "crystalball") tag = Form("cball (tail)");
            else if (c.modelName == "dscb")        tag = Form("dscb (2tail)");
            else if (c.modelName == "doublegauss") tag = Form("dgauss (RMS)");
            else if (std::abs(c.nSigmaLo - c.nSigmaHi) > 1e-9)
                         tag = Form("mu -%.1f +%.1f sig", c.nSigmaLo, c.nSigmaHi);
            else        tag = Form("mu+/-%.1f sig conv", c.nSigma);

            std::string marker = ((int)i == bestIdx) ? " <- BEST" : "";
            std::cout << Form("  │ %-16s │ %.4f +/- %.4f │ %.4f +/- %.4f │ %6.2f / %-3d   │ [%.3f, %.3f]   │%s",
                              tag.c_str(),
                              c.mean, c.meanErr, c.sigma, c.sigmaErr,
                              c.chi2, c.ndf,
                              c.fitLo, c.fitHi,
                              marker.c_str()) << std::endl;
        }
        std::cout << "    +---------------------+---------------------+---------------------+------------------------+\n" << std::endl;

        // --- Step 5: Optional scan canvas ---
        if (!scanPdfPath.empty()) {
            TCanvas* scanCanvas = new TCanvas(("cscan_" + nameTag).c_str(), "Fit Window Scan", 800, 600);
            hist->SetLineColor(kBlack);
            hist->SetFillStyle(0);
            hist->Draw("HIST");
            hist->SetMaximum(hist->GetMaximum() * 1.4);
            hist->SetMinimum(0);

            int scanColors[] = {kGray+1, kBlue, kGreen+2, kOrange+1, kMagenta, kCyan+1};
            int nColors = sizeof(scanColors) / sizeof(scanColors[0]);

            TLegend* leg = new TLegend(0.62, 0.58, 0.92, 0.90);
            leg->SetBorderSize(0);
            leg->SetTextSize(0.022);
            leg->AddEntry(hist, "Histogram", "l");

            for (size_t i = 0; i < candidates.size(); ++i) {
                const auto& c = candidates[i];
                // Draw the ACTUAL fitted model for this candidate (DSCB/CB/Bukin
                // or Gaussian) - not a Gaussian reconstructed from (amp,mean,sigma).
                // So an alt-model run shows the curve that is actually selected.
                TF1* draw_func = MakeModelTF1(c, Form("draw_%zu_%s", i, nameTag.c_str()));
                int col = scanColors[i % nColors];
                int style = ((int)i == bestIdx) ? 1 : 2;
                int width = ((int)i == bestIdx) ? 3 : 1;
                draw_func->SetLineColor(col);
                draw_func->SetLineStyle(style);
                draw_func->SetLineWidth(width);
                draw_func->Draw("same");

                std::string lbl;
                if (i == 0) lbl = Form("Seed [%.2f,%.2f]", seedLo, seedHi);
                else if (std::abs(c.nSigmaLo - c.nSigmaHi) > 1e-9)
                            lbl = Form("#mu-%.1f/+%.1f#sigma [%.2f,%.2f]",
                                       c.nSigmaLo, c.nSigmaHi, c.fitLo, c.fitHi);
                else        lbl = Form("#mu#pm%.1f#sigma conv [%.2f,%.2f]", c.nSigma, c.fitLo, c.fitHi);
                // For alt-model candidates, prefix the model so the curve is
                // unambiguous (the window naming refers to the seed Gaussian it
                // was fit over). Short tags keep the legend readable.
                if (c.modelName == "bukin")            lbl = "Bukin " + lbl;
                else if (c.modelName == "crystalball") lbl = "CB " + lbl;
                else if (c.modelName == "dscb")        lbl = "DSCB " + lbl;
                else if (c.modelName == "doublegauss") lbl = "DG " + lbl;
                lbl += Form("  #chi^{2}/ndf=%.2f/%d", c.chi2, c.ndf);
                if ((int)i == bestIdx) lbl += "  #bf{BEST}";
                leg->AddEntry(draw_func, lbl.c_str(), "l");
            }
            leg->Draw();

            TLatex tx;
            tx.SetNDC(); tx.SetTextFont(42); tx.SetTextSize(0.030);
            tx.DrawLatex(0.15, 0.85, label.c_str());

            const auto& best = candidates[bestIdx];
            tx.SetTextSize(0.028);
            tx.DrawLatex(0.15, 0.80, Form("Best (%s): #mu=%.4f  #sigma=%.4f  #chi^{2}/ndf=%.2f/%d",
                                          ModelDisplayName(best.modelName).c_str(),
                                          best.mean, best.sigma, best.chi2, best.ndf));

            scanCanvas->SaveAs(scanPdfPath.c_str());
            delete scanCanvas;
        }

        return candidates[bestIdx];
    }

    static constexpr int MIN_BINS_FOR_FIT = 4;
    static constexpr int MIN_ENTRIES_SKIP = 20;

    inline TH1F* AdaptiveRebin(TH1F* h,
        int                target_per_bin = 4,
                                const std::string& name_suffix = "_rebinned")
    {
        if (!h) return nullptr;

        const double eff_entries = h->GetEffectiveEntries();
        const int    nbins_orig  = h->GetNbinsX();

        TH1F* hc = (TH1F*)h->Clone((std::string(h->GetName()) + name_suffix).c_str());
        hc->SetDirectory(nullptr);

        if (eff_entries < MIN_ENTRIES_SKIP) {
            hc->SetTitle((std::string(hc->GetTitle()) + " [SPARSE]").c_str());
            return hc;
        }

        int desired_bins = std::max(MIN_BINS_FOR_FIT,
                                    (int)std::floor(eff_entries / target_per_bin));

        int best_factor = 1;
        for (int f : {2, 4, 8, 16, 32}) {
            if (nbins_orig % f != 0) continue;
            if (nbins_orig / f < MIN_BINS_FOR_FIT) break;
            if (nbins_orig / f >= desired_bins) best_factor = f;
        }

        if (best_factor > 1) {
            hc->Rebin(best_factor);
        }

        int peak_lo  = hc->FindBin(0.8);
        int peak_hi  = hc->FindBin(1.2);
        int nonempty = 0;
        for (int b = peak_lo; b <= peak_hi; b++)
            if (hc->GetBinContent(b) > 0) nonempty++;

        if (nonempty < MIN_BINS_FOR_FIT) {
            hc->SetTitle((std::string(hc->GetTitle()) + " [SPARSE]").c_str());
        }

        return hc;
    }

    inline bool IsSparse(const TH1F* h) {
        return h && std::string(h->GetTitle()).find("[SPARSE]") != std::string::npos;
    }

    struct FitOrSkipResult {
        GausFitResult fit;
        bool          skipped  = false;
        int           nBinsUsed = 0;
        double        effEntries = 0;
        int           rebin_factor = 1;
    };

    inline std::pair<double,double> FindPeakSeed(
        TH1F*  h,
        double searchLo  = 0.8,
        double searchHi  = 1.2,
        double halfWidth = 0.08)
    {
        if (!h) return {0.9, 1.1};

        int binLo = h->FindBin(searchLo);
        int binHi = h->FindBin(searchHi);
        int peakBin = binLo;
        double peakContent = -1;
        for (int b = binLo; b <= binHi; b++) {
            if (h->GetBinContent(b) > peakContent) {
                peakContent  = h->GetBinContent(b);
                peakBin      = b;
            }
        }
        double peakX = h->GetBinCenter(peakBin);
        return { peakX - halfWidth, peakX + halfWidth };
    }

    inline FitOrSkipResult FitOrSkip(
        TH1F*              h,
        const std::string& name_tag,
        const std::string& label,
        const std::string& scan_pdf      = "",
        int                target_per_bin = 4,
        const std::string& warn_prefix   = "",
        double             seedLo        = 0.8,
        double             seedHi        = 1.2,
        double             maxNSigma     = 2.,
        bool               usePeakSeed  = false,
        double             peakHalfWidth = 0.08)
    {
        FitOrSkipResult out;
        if (!h) { out.skipped = true; return out; }

        out.effEntries = h->GetEffectiveEntries();

        TH1F* hc = AdaptiveRebin(h, target_per_bin, "_fit_" + name_tag);
        out.nBinsUsed    = hc->GetNbinsX();
        out.rebin_factor = h->GetNbinsX() / out.nBinsUsed;

        if (IsSparse(hc)) {
            out.skipped = true;
            std::cout << "  [SKIP] " << warn_prefix
                      << " — too few entries (" << (int)out.effEntries
                      << ") after rebinx" << out.rebin_factor
                      << " (" << out.nBinsUsed << " bins) — fit skipped\n";
            delete hc;
            return out;
        }

        if (out.rebin_factor > 1) {
            std::cout << "  [REBIN] " << warn_prefix
                      << " — " << h->GetNbinsX() << " -> " << out.nBinsUsed
                      << " bins (entries=" << (int)out.effEntries
                      << ", factor=" << out.rebin_factor << ")\n";
        }

        if (usePeakSeed) {
            auto [lo, hi] = FindPeakSeed(hc, 0.8, 1.2, peakHalfWidth);
            seedLo = lo;
            seedHi = hi;
            std::cout << "  [PEAK SEED] " << warn_prefix
                      << Form(" — peak found at %.3f -> seed [%.3f, %.3f]\n",
                              (lo+hi)/2.0, lo, hi);
        }

        out.fit = IterativeGausFit(hc, name_tag, label, scan_pdf, seedLo, seedHi, maxNSigma);
        delete hc;
        return out;
    }

    inline GausFitResult FitSingleWindow(
        TH1F*              h,
        const std::string& nameTag,
        double             searchLo     = 0.8, //0.85
        double             searchHi     = 1.2, //1.15
        double             fitHalfWidth = 0.15) //0.15
    {
        GausFitResult result;
        if (!h) return result;

        auto [peakLo, peakHi] = FindPeakSeed(h, searchLo, searchHi, fitHalfWidth);
        double peakX = (peakLo + peakHi) / 2.0;

        double xMin = h->GetXaxis()->GetXmin();
        double xMax = h->GetXaxis()->GetXmax();
        double lo = std::max(xMin, peakX - fitHalfWidth);
        double hi = std::min(xMax, peakX + fitHalfWidth);

        TF1* f = new TF1(("sw_" + nameTag).c_str(), "gaus", lo, hi);
        f->SetParameter(0, h->GetBinContent(h->FindBin(peakX)));
        f->SetParameter(1, peakX);
        f->SetParameter(2, (hi - lo) / 4.0);
        h->Fit(f, "RQ0");

        result.amp      = f->GetParameter(0);  result.ampErr   = f->GetParError(0);
        result.mean     = f->GetParameter(1);  result.meanErr  = f->GetParError(1);
        result.sigma    = std::abs(f->GetParameter(2));
        result.sigmaErr = f->GetParError(2);
        result.chi2     = f->GetChisquare();
        result.ndf      = f->GetNDF();
        result.chi2ndf  = (result.ndf > 0) ? result.chi2 / result.ndf : 999;
        result.fitLo    = lo;
        result.fitHi    = hi;
        result.nSigma   = 0;
        result.valid    = (result.ndf > 0 && result.sigma > 0
                           && std::isfinite(result.chi2ndf));
        delete f;
        return result;
    }

    inline FitOrSkipResult FitSingleWindowOrSkip(
        TH1F*              h,
        const std::string& name_tag,
        const std::string& warn_prefix   = "",
        int                target_per_bin = 4,
        double             searchLo      = 0.8,
        double             searchHi      = 1.2,
        double             fitHalfWidth  = 0.15)
    {
        FitOrSkipResult out;
        if (!h) { out.skipped = true; return out; }

        out.effEntries = h->GetEffectiveEntries();

        TH1F* hc = AdaptiveRebin(h, target_per_bin, "_sw_" + name_tag);
        out.nBinsUsed    = hc->GetNbinsX();
        out.rebin_factor = h->GetNbinsX() / out.nBinsUsed;

        if (IsSparse(hc)) {
            out.skipped = true;
            std::cout << "  [SKIP] " << warn_prefix
                      << " — too few entries (" << (int)out.effEntries
                      << ") after rebin×" << out.rebin_factor << " — fit skipped\n";
            delete hc;
            return out;
        }

        if (out.rebin_factor > 1)
            std::cout << "  [REBIN] " << warn_prefix
                      << " — " << h->GetNbinsX() << " -> " << out.nBinsUsed
                      << " bins (factor=" << out.rebin_factor << ")\n";

        out.fit = FitSingleWindow(hc, name_tag, searchLo, searchHi, fitHalfWidth);

        if (out.fit.valid && out.fit.chi2ndf > 3.0)
            std::cout << "  [WARN] " << warn_prefix
                      << Form(" — single-window fit chi2/ndf=%.2f (peak=%.3f, sigma=%.4f)\n",
                              out.fit.chi2ndf, out.fit.mean, out.fit.sigma);

        delete hc;
        return out;
    }

    // ============================================================================
    // b-JER data/MC ratio with correct error propagation
    // ============================================================================
    // Computes R_JER = sigma_insitu / sigma_insitu^MC with
    //   sigma_insitu    = sqrt(sigma_data^2 - sigma_truth^2)
    //   sigma_insitu^MC = sqrt(sigma_MC^2   - sigma_truth^2)
    // sigma_truth appears in both and is propagated as one variable.
    // valid is false if either difference of squares is negative.

    // ============================================================================
    struct JERRatio {
        double ratio    = 0;
        double ratioErr = 0;
        double sig_insitu      = 0, sig_insitu_err    = 0;
        double sig_MC_insitu   = 0, sig_MC_insitu_err = 0;
        bool   valid    = false;
    };

    inline JERRatio ComputeJERRatio(double sig_data,  double sig_data_err,
                                    double sig_mc,    double sig_mc_err,
                                    double sig_truth, double sig_truth_err)
    {
        JERRatio r;
        const double A = sig_data,  sA = sig_data_err;
        const double B = sig_mc,    sB = sig_mc_err;
        const double T = sig_truth, sT = sig_truth_err;

        const double N2 = A*A - T*T;   // sigma_insitu^2
        const double D2 = B*B - T*T;   // sigma^MC_insitu^2
        if (N2 <= 0 || D2 <= 0) return r;   

        const double N = std::sqrt(N2);
        const double D = std::sqrt(D2);

        r.sig_insitu        = N;
        r.sig_MC_insitu     = D;
        r.sig_insitu_err    = std::sqrt(std::pow(A*sA/N, 2) + std::pow(T*sT/N, 2));
        r.sig_MC_insitu_err = std::sqrt(std::pow(B*sB/D, 2) + std::pow(T*sT/D, 2));

        r.ratio = N / D;

        // Partial derivatives of R = N/D w.r.t. the three independent inputs:
        //   dR/dA =  A/(N*D)
        //   dR/dB = -N*B/D^3
        //   dR/dT =  T*(N^2 - D^2)/(N*D^3)
        const double dA =  A / (N * D);
        const double dB = -N * B / (D*D*D);
        const double dT =  T * (N2 - D2) / (N * D*D*D);

        r.ratioErr = std::sqrt(dA*dA*sA*sA + dB*dB*sB*sB + dT*dT*sT*sT);
        r.valid    = true;
        return r;
    }

    // ============================================================================
    // JES data/MC ratio (response means, simple uncorrelated error propagation)
    // ============================================================================
    // R_JES = mu_data / mu_MC
    // Independent inputs (no shared variable), so naive propagation is correct:
    //   (sigma_R / R)^2 = (sigma_mu_data/mu_data)^2 + (sigma_mu_MC/mu_MC)^2
    // ============================================================================
    struct JESRatio {
        double ratio    = 0;
        double ratioErr = 0;
        bool   valid    = false;
    };

    inline JESRatio ComputeJESRatio(double mu_data, double mu_data_err,
                                    double mu_mc,   double mu_mc_err)
    {
        JESRatio r;
        if (mu_mc == 0 || mu_data == 0) return r;
        r.ratio = mu_data / mu_mc;
        r.ratioErr = std::abs(r.ratio) * std::sqrt(
            std::pow(mu_data_err / mu_data, 2) +
            std::pow(mu_mc_err   / mu_mc,   2));
        r.valid = true;
        return r;
    }

}  // namespace FitUtils


// ============================================================================
// FitOutputUtils — shared fit-parameter I/O
// ============================================================================
namespace FitOutputUtils {

    // The fit_params_*.txt files read by the ratio macros.
    inline void SaveFitParameters(const char* filename, const char* fitType,
                                   double mean,  double meanErr,
                                   double sigma, double sigmaErr,
                                   double amp,   double ampErr,
                                   double chi2,  double ndf) {
        std::ofstream out(filename);
        if (!out.is_open()) {
            Log::error("SaveFitParameters", std::string("Cannot open: ") + filename);
            return;
        }
        double chi2ndf = (ndf > 0) ? chi2 / ndf : -1.0;
        out << "Fit Type: "  << fitType << "\n"
            << "Mean: "      << mean      << " +/- " << meanErr      << "\n"
            << "Sigma: "     << sigma     << " +/- " << sigmaErr     << "\n"
            << "Amplitude: " << amp       << " +/- " << ampErr       << "\n"
            << "Chi2/NDF: "  << chi2ndf   << " ("    << chi2
                             << "/"       << (int)ndf << ")\n";
        out.close();
        Log::info("SaveFitParameters", std::string("Saved: ") + filename);
    }

    inline void SaveSkippedSentinel(const char* filename, const char* fitType) {
        std::ofstream out(filename);
        if (!out.is_open()) return;
        out << "Fit Type: " << fitType << " [SKIPPED]\n"
            << "Mean: 0 +/- 0\n"
            << "Sigma: 0 +/- 0\n"
            << "Amplitude: 0 +/- 0\n"
            << "Chi2/NDF: -1 (0/0)\n"
            << "SKIPPED: 1\n";
        out.close();
    }

}  // namespace FitOutputUtils


#endif // COMMON_UTILS_H
