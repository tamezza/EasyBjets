/**
 * ============================================================================
 *  BinningConfig.h - histogram binning for the mass and response fits
 * ============================================================================
 * ===========================================
 * @author:    Thiziri AMEZZA
 * @email:     thiziri.amezza@cern.ch
 * @date:      10 March 2026 
 * ===========================================
 *  OVERVIEW
 *  --------
 *  The fit results depend noticeably on the binning, so instead of hardcoding
 *  one choice the bin count can be derived from the distribution itself.
 *  The mode is picked with EASYBJETS_BINNING (case-insensitive):
 *
 *   fixed        bin count from VariantConfig (default)
 *   fd           Freedman-Diaconis, 2 IQR N^-1/3
 *   scott        Scott, 3.49 sigma N^-1/3
 *   scottjes     Scott applied as a rebin factor of the 100-bin scratch
 *                histogram, like JES_BalanceFitter::OptimalRebin; keeps bin
 *                edges aligned between variants
 *   sturges      1 + log2(N)
 *   peakaware    pre-fit the peak, bin width = sigma / PEAKAWARE_K
 *   variable     fine bins under the peak, coarse in the sidebands
 *   fitfriendly  about EASYBJETS_BINNING_CORE_BINS bins (default 10) across
 *                the +-2 sigma core, sigma taken from the core only. Gives
 *                every variant the same number of bins inside its own peak,
 *                which is what we want when comparing fit results between them.
 * 
 * USAGE
 *  -----
 *    export EASYBJETS_BINNING=peakaware
 *    root 'Mass_Analysis.cxx("ALL")'
 *
 *  Unknown or empty values fall back to the compile-time default (fixed).
 *
 * The mode name is appended to the output folder so runs with different
 * binnings do not overwrite each other. Related knobs:
 *   EASYBJETS_BINNING_MAX_BINS         global cap on the bin count
 *   EASYBJETS_BINNING_MAX_BINS_TRUTH   same, for the truth-side histograms
 *   EASYBJETS_BINNING_CORE_BINS        target core bins for fitfriendly
 *   EASYBJETS_RESPONSE_BINNING_MODE    static (default) or dynamic, see below
 *   EASYBJETS_BINNING_TRUTH            how the truth histograms are binned
 *                                      in dynamic mode, see TruthPolicy
 *
 * For production we use fixed binning; the other modes were used to decide
 * on the fixed values and to check the binning systematic.
 * 
 * ==========================================================================
 * @version:  1.0
 * ========================================================================== 
 */

#ifndef EASYBJETS_BINNING_CONFIG_H
#define EASYBJETS_BINNING_CONFIG_H

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
#include <TH1F.h>
#include <TF1.h>

namespace BinningConfig {

    // -------- EDITABLE PARAMETERS -------------------------------------------

    // Fixed mode: uniform bin count per channel 
    constexpr int DB_FIXED_BINS  = 30;
    constexpr int MJB_FIXED_BINS = 25;

    // peakaware: K bins per sigma; the pre-fit sigma is clamped to a sane range.
    constexpr double PEAKAWARE_K = 5.0;
    constexpr double PEAKAWARE_MIN_SIGMA = 3.0;   
    constexpr double PEAKAWARE_MAX_SIGMA = 40.0; 
    constexpr int PEAKAWARE_PREFIT_BINS = 150;
    constexpr double PEAKAWARE_PREFIT_CENTER_GUESS = 90.0;  

    // ---- variable mode: fine bins in peak region, coarse outside ----
    // Peak region ~ [70, 110] GeV with 1.5 GeV bins; sidebands 5 GeV.
    inline std::vector<double> VariableEdgesDB() {
        std::vector<double> e;
        for (double m = 50.0; m < 70.0;  m += 5.0)  e.push_back(m);   // 50, 55, 60, 65
        for (double m = 70.0; m <= 110.0; m += 1.5) e.push_back(m);   // 70, 71.5, ..., 110
        for (double m = 115.0; m <= 200.0; m += 5.0) e.push_back(m);  // 115, 120, ..., 200
        return e;
    }

    inline std::vector<double> VariableEdgesMJB() {
        std::vector<double> e;
        for (double m = 50.0; m < 70.0;  m += 5.0)  e.push_back(m);
        for (double m = 70.0; m <= 110.0; m += 1.5) e.push_back(m);
        for (double m = 115.0; m <= 150.0; m += 5.0) e.push_back(m);
        return e;
    }

    // ------------------------------------------------------------------------
    // MODE SELECTION
    // ------------------------------------------------------------------------
    enum class Mode { Fixed, FD, Scott, ScottJES, Sturges, PeakAware, Variable, FitFriendly };

    #ifndef EASYBJETS_DEFAULT_BINNING
    #define EASYBJETS_DEFAULT_BINNING ::BinningConfig::Mode::Fixed
    #endif

    inline bool ParseMode(const std::string& s_in, Mode& out) {
        std::string s(s_in);
        for (auto& ch : s) ch = (char)std::tolower((unsigned char)ch);
        if (s == "fixed")     { out = Mode::Fixed;     return true; }
        if (s == "fd")        { out = Mode::FD;        return true; }
        if (s == "scott")     { out = Mode::Scott;     return true; }
        if (s == "scottjes" || s == "scott-jes" || s == "jes")
                              { out = Mode::ScottJES;  return true; }
        if (s == "sturges")   { out = Mode::Sturges;   return true; }
        if (s == "peakaware" || s == "peak-aware" || s == "peak")
                              { out = Mode::PeakAware; return true; }
        if (s == "variable" || s == "var")
                              { out = Mode::Variable;  return true; }
        if (s == "fitfriendly" || s == "fit-friendly" || s == "ff")
                              { out = Mode::FitFriendly; return true; }
        return false;
    }

    inline Mode CurrentMode() {
        static Mode cached = [](){
            const char* env = std::getenv("EASYBJETS_BINNING");
            if (env && *env) {
                Mode parsed;
                if (ParseMode(env, parsed)) return parsed;
                std::cerr << "[BinningConfig] WARN: unknown EASYBJETS_BINNING='"
                          << env << "', using default (fixed)\n";
            }
            return EASYBJETS_DEFAULT_BINNING;
        }();
        return cached;
    }

    // ===========================================================
    // Truth-histogram binning policy (EASYBJETS_BINNING_TRUTH)
    // ===========================================================
    // Controls how the truth-side response histograms (truth/ref, reco/truth) are
    // binned when EASYBJETS_RESPONSE_BINNING_MODE=dynamic. Ignored otherwise.
    //
    //   Match (default)    - same binning as the reco histogram (default)
    //                        EASYBJETS_BINNING_TRUTH=match (or unset).
    //
    //   Auto               - same algorithm as EASYBJETS_BINNING, run on the truth distriution, 
    //                        (e.g. pT_reco/pT_truth).
    //                        Bin edges DIFFER from reco. Each distribution
    //                        gets a binning optimal for its own width.
    //                        EASYBJETS_BINNING_TRUTH=auto.
    //
    //   Fixed              - Use cfg.GetNBinsTruth() (e.g. 50-bin),
    //                        EASYBJETS_BINNING_TRUTH=fixed.
    //
    //   Algorithm name     - a specific algorithm, e.g. EASYBJETS_BINNING_TRUTH=scott (etc.).
    //
    // Note: in static mode (EASYBJETS_RESPONSE_BINNING_MODE unset), this
    // policy is ignored — truth always uses cfg.GetNBinsTruth(). This env
    // var has effect only when dynamic mode is active.
    enum class TruthPolicy { Match, Auto, Fixed, Algorithm };

    // alg_out is only meaningful for TruthPolicy::Algorithm
    inline TruthPolicy CurrentTruthPolicy(Mode& alg_out) {
        static struct Cache {
            TruthPolicy policy;
            Mode        alg;
        } cached = [](){
            Cache c{TruthPolicy::Match, Mode::Fixed};
            const char* env = std::getenv("EASYBJETS_BINNING_TRUTH");
            if (!env || !*env) return c;   // default: match
            std::string s(env);
            for (auto& ch : s) ch = (char)std::tolower((unsigned char)ch);
            if (s == "match" || s == "shared" || s == "same")
                { c.policy = TruthPolicy::Match;  return c; }
            if (s == "auto")
                { c.policy = TruthPolicy::Auto;   return c; }
            if (s == "fixed" || s == "legacy" || s == "static")
                { c.policy = TruthPolicy::Fixed;  return c; }
            Mode parsed;
            if (ParseMode(s, parsed)) {
                c.policy = TruthPolicy::Algorithm;
                c.alg    = parsed;
                return c;
            }
            std::cerr << "[BinningConfig] WARN: unknown EASYBJETS_BINNING_TRUTH='"
                      << env << "' (valid: match|auto|fixed|<algorithm>), "
                      << "falling back to 'match'\n";
            return c;
        }();
        alg_out = cached.alg;
        return cached.policy;
    }

    inline const char* TruthPolicyName(TruthPolicy p) {
        switch (p) {
            case TruthPolicy::Match:     return "match";
            case TruthPolicy::Auto:      return "auto";
            case TruthPolicy::Fixed:     return "fixed";
            case TruthPolicy::Algorithm: return "algorithm";
        }
        return "unknown";
    }

    // Short tag for filename suffixes and log messages.
    inline std::string ModeTag(Mode m = CurrentMode()) {
        switch (m) {
            case Mode::Fixed:     return "fixed";
            case Mode::FD:        return "fd";
            case Mode::Scott:     return "scott";
            case Mode::ScottJES:  return "scottjes";
            case Mode::Sturges:   return "sturges";
            case Mode::PeakAware: return "peak";
            case Mode::Variable:  return "var";
            case Mode::FitFriendly: return "fitfriendly";
        }
        return "unknown";
    }

    // ------------------------------------------------------------------------
    // BINNING RESULT
    // ------------------------------------------------------------------------
    // If edges is empty the histogram is uniform with nbins bins, otherwise
    // the edges are used as they are.

    struct Result {
        int                 nbins   = 0;
        std::vector<double> edges;           // empty -> uniform
        std::string         tag;             // ModeTag() at construction time
        std::string         caption;         // printed on the plots
        double              sigma_peak = 0;  // peakaware only
    };

    // ------------------------------------------------------------------------
    // Helpers for effective sample size & quantile stats
    // ------------------------------------------------------------------------

    // EASYBJETS_BINNING_MAX_BINS lowers the cap passed by the caller;
    //  it never raises it. The env var acts as a global ceiling:
    //  if the call site already requested a smaller max_bins, 
    // the env var has no effect (we use std::min, not std::max).
    inline int RuntimeMaxBins(int default_max) {
        static const int v = []() {
            const char* e = std::getenv("EASYBJETS_BINNING_MAX_BINS");
            if (!e || !*e) return -1;
            try { int n = std::stoi(e); if (n > 1 && n < 10000) return n; }
            catch (...) {}
            return -1;
        }();
        return (v > 0) ? std::min(v, default_max) : default_max;
    }

    // Truth cap: EASYBJETS_BINNING_MAX_BINS_TRUTH if set, else the reco cap.
    // Resolution order when computing the truth cap:
    //   1. EASYBJETS_BINNING_MAX_BINS_TRUTH (if set, takes precedence)
    //   2. EASYBJETS_BINNING_MAX_BINS (falls back to the reco cap)
    //   3. The call-site default (no env override)
    //
    // This means: if you only set EASYBJETS_BINNING_MAX_BINS, BOTH reco and
    // truth use that cap. If you set BOTH env vars, they apply independently. 
    inline int RuntimeMaxBinsTruth(int default_max) {
        static const int v = []() {
            const char* e = std::getenv("EASYBJETS_BINNING_MAX_BINS_TRUTH");
            if (!e || !*e) return -1;
            try { int n = std::stoi(e); if (n > 1 && n < 10000) return n; }
            catch (...) {}
            return -1;
        }();
        if (v > 0) return std::min(v, default_max);
        return RuntimeMaxBins(default_max);
    }

    // Kish's effective N for weighted histograms. (sum w)^2 / sum w^2.
    inline double EffectiveN(TH1F* h) {
        double sumw  = h->GetSumOfWeights();
        double sumw2 = 0.0;
        if (h->GetSumw2N() > 0) {
            for (int i = 1; i <= h->GetNbinsX(); ++i)
                sumw2 += h->GetSumw2()->At(i);
        }
        if (sumw2 > 0 && sumw > 0) return (sumw * sumw) / sumw2;
        return static_cast<double>(h->GetEntries());
    }

    // ------------------------------------------------------------------------
    // Per-mode bin-count computations
    // ------------------------------------------------------------------------

    // Freedman–Diaconis: 2·IQR·N^{-1/3}
    inline int ComputeFD(TH1F* h, double xmin, double xmax, int min_bins, int max_bins) {
        double neff = EffectiveN(h);
        if (neff < 2) return min_bins;
        double q[3]; double p[3] = {0.25, 0.5, 0.75};
        h->GetQuantiles(3, q, p);
        double iqr = q[2] - q[0];
        if (iqr < 1e-6) return min_bins;
        double w = 2.0 * iqr * std::pow(neff, -1.0/3.0);
        int n = static_cast<int>(std::ceil((xmax - xmin) / w));
        return std::max(min_bins, std::min(max_bins, n));
    }

    // Scott's rule: 3.49 sigma N^-1/3
    inline int ComputeScott(TH1F* h, double xmin, double xmax, int min_bins, int max_bins) {
        double neff = EffectiveN(h);
        if (neff < 2) return min_bins;
        double sigma = h->GetRMS();
        if (sigma < 1e-6) return min_bins;
        double w = 3.49 * sigma * std::pow(neff, -1.0/3.0);
        int n = static_cast<int>(std::ceil((xmax - xmin) / w));
        return std::max(min_bins, std::min(max_bins, n));
    }

    inline double CoreRMS(TH1* h);   

    // Fit-friendly bin count: 
    // Bin width = 4 sigma_core / Nwin, i.e. Nwin bins inside +-2 sigma.
    inline int ComputeFitFriendly(TH1F* h, double xmin, double xmax, int min_bins, int max_bins) {
        if (!h) return min_bins;
        double neff = EffectiveN(h);
        if (neff < 2) return min_bins;
        double sigma = CoreRMS(h);
        if (sigma < 1e-6) return min_bins;
        static const double Nwin = []() {
            const char* e = std::getenv("EASYBJETS_BINNING_CORE_BINS");
            if (e && *e) { try { double v = std::stod(e); if (v >= 2.0 && v <= 200.0) return v; } catch (...) {} }
            return 10.0;
        }();
        double w = (4.0 * sigma) / Nwin;                 
        int n = static_cast<int>(std::ceil((xmax - xmin) / w));
        return std::max(min_bins, std::min(max_bins, n));
    }

    // RMS of the core only: recompute the RMS three times inside +-2.5 sigma
    // of the current mean. For a narrow peak with wide tails (the reco/truth
    // response for instance) the full RMS is roughly twice the core width.
    inline double CoreRMS(TH1* h) {
        if (!h) return 0.0;
        double mean = h->GetMean(), rms = h->GetRMS();
        if (rms <= 1e-6 || h->GetEntries() < 10) return rms;
        const int a0 = h->GetXaxis()->GetFirst(), a1 = h->GetXaxis()->GetLast();
        for (int it = 0; it < 3; ++it) {
            h->GetXaxis()->SetRangeUser(mean - 2.5 * rms, mean + 2.5 * rms);
            double m = h->GetMean(), r = h->GetRMS();
            if (r > 1e-6) { mean = m; rms = r; }
        }
        h->GetXaxis()->SetRange(a0, a1);   // restore the original axis range
        return rms;
    }

    // Fit-friendly bin count using the tail-robust CORE width (CoreRMS) instead of
    // the full RMS. Use for narrow peaked distributions (reco/truth) whose full RMS
    // would otherwise under-bin the core (giving fits with ndf=1-2). Same target:
    // ~EASYBJETS_BINNING_CORE_BINS bins across the +-2. sigma core.
    inline int ComputeFitFriendlyCore(TH1F* h, double xmin, double xmax, int min_bins, int max_bins) {
        if (!h) return min_bins;
        double neff = EffectiveN(h);
        if (neff < 2) return min_bins;
        double sigma = CoreRMS(h);
        if (sigma < 1e-6) return min_bins;
        static const double Nwin = []() {
            const char* e = std::getenv("EASYBJETS_BINNING_CORE_BINS");
            if (e && *e) { try { double v = std::stod(e); if (v >= 2.0 && v <= 200.0) return v; } catch (...) {} }
            return 10.0;
        }();
        double w = (4.0 * sigma) / Nwin;
        int n = static_cast<int>(std::ceil((xmax - xmin) / w));
        return std::max(min_bins, std::min(max_bins, n));
    }

    // Scott JES-style rebin (mirrors ATLAS JES_BalanceFitter::OptimalRebin).
    // Uses the same Scott formula (3.5·sigma·N^{-1/3}) but applied as a REBIN
    // factor on the scratch histogram, not as an arbitrary new bin count.
    inline int ComputeScottJES(TH1F* h, int min_bins, int max_bins) {
        double neff = EffectiveN(h);
        if (neff < 2) return min_bins;
        double sigma = h->GetRMS();
        if (sigma < 1e-6) return min_bins;

        double optWidth = 3.5 * sigma * std::pow(neff, -1.0/3.0);
        const int Nbins = h->GetNbinsX();
        const double range = h->GetBinLowEdge(Nbins + 1) - h->GetBinLowEdge(1);
        if (Nbins <= 1 || range <= 0) return min_bins;

        int rebin = 1;
        for (int i = 1; i < Nbins; ++i) {
            if (Nbins % i != 0) continue;          
            double binWidth = (range / Nbins) * i;
            if (binWidth < optWidth) rebin = i;    
            if (binWidth > optWidth) break;
        }

        int nbins = Nbins / rebin;
        return std::max(min_bins, std::min(max_bins, nbins));
    }

    // Sturges' rule: 1 + log2(N)
    inline int ComputeSturges(TH1F* h, int min_bins, int max_bins) {
        double neff = EffectiveN(h);
        if (neff < 2) return min_bins;
        int n = static_cast<int>(std::ceil(1.0 + std::log2(neff)));
        return std::max(min_bins, std::min(max_bins, n));
    }

    // Peak-aware: pre-fit a Gaussian on a coarsely-binned hist, then rebin
    // so bin_width = igma / K per bin. K is a tunable parameter.
    // If fit failure, returns a fallback FD bin count with a warning.
    inline int ComputePeakAware(TH1F* h, double xmin, double xmax,
                                 int min_bins, int max_bins,
                                 double& sigma_out) {
        double neff = EffectiveN(h);
        if (neff < 20) {
            std::cerr << "[BinningConfig] peakaware: too few events ("
                      << neff << "), falling back to FD\n";
            return ComputeFD(h, xmin, xmax, min_bins, max_bins);
        }

        TH1F* hc = (TH1F*)h->Clone("h_peakaware_prefit");
        hc->SetDirectory(nullptr);
        hc->Rebin(std::max(1, h->GetNbinsX() / PEAKAWARE_PREFIT_BINS));

        // Initial range: narrow window around prior peak guess. If the hist
        // doesn't reach the guess, use the middle of [xmin, xmax].
        double guess = PEAKAWARE_PREFIT_CENTER_GUESS;
        if (guess < xmin || guess > xmax) guess = 0.5 * (xmin + xmax);
        double lo = std::max(xmin, guess - 20.0);
        double hi = std::min(xmax, guess + 20.0);

        TF1 g("g_peakaware", "gaus", lo, hi);
        g.SetParameter(1, guess);
        g.SetParameter(2, 10.0);  // sigma seed
        int status = hc->Fit(&g, "QNR");
        double sigma = g.GetParameter(2);
        delete hc;

        if (status != 0 || sigma < PEAKAWARE_MIN_SIGMA || sigma > PEAKAWARE_MAX_SIGMA) {
            std::cerr << "[BinningConfig] peakaware: pre-fit failed (status="
                      << status << ", sigma=" << sigma << "), falling back to FD\n";
            return ComputeFD(h, xmin, xmax, min_bins, max_bins);
        }

        sigma_out = sigma;
        double w = sigma / PEAKAWARE_K;
        int n = static_cast<int>(std::ceil((xmax - xmin) / w));
        return std::max(min_bins, std::min(max_bins, n));
    }

    // ------------------------------------------------------------------------
    // MAIN ENTRY POINT
    // ------------------------------------------------------------------------
    // Computes the binning for a single pT bin, using the current mode.
    // `h` is a finely-binned MC scratch histogram filled with the pT slice.
    // `channel` is "DB" or "MJB" - used to pick per-channel defaults.
    // `min_bins` / `max_bins` clamp results for the data-driven modes.
    //
    // Returns a Result with either a uniform bin count (edges.empty()) or
    // a variable edge list (edges.size() > 1).
    inline Result Compute(TH1F* h, const std::string& channel,
                           double xmin, double xmax,
                           int min_bins, int max_bins) {
        Mode m = CurrentMode();
        Result r;
        r.tag = ModeTag(m);

        max_bins = RuntimeMaxBins(max_bins);

        switch (m) {
            case Mode::Fixed: {
                int variant_bins = max_bins;
                if (variant_bins <= 0 || variant_bins > 200) {
                    variant_bins = (channel == "MJB") ? MJB_FIXED_BINS : DB_FIXED_BINS;
                    r.nbins   = variant_bins;
                    r.caption = "binning: fixed [channel-default], "
                              + std::to_string(r.nbins) + " bins";
                } else {
                    r.nbins   = variant_bins;
                    r.caption = "binning: fixed [VariantConfig], "
                              + std::to_string(r.nbins) + " bins";
                }
                break;
            }
            case Mode::FD: {
                r.nbins = ComputeFD(h, xmin, xmax, min_bins, max_bins);
                r.caption = "binning: FD, " + std::to_string(r.nbins) + " bins";
                break;
            }
            case Mode::Scott: {
                r.nbins = ComputeScott(h, xmin, xmax, min_bins, max_bins);
                r.caption = "binning: Scott, " + std::to_string(r.nbins) + " bins";
                break;
            }
            case Mode::ScottJES: {
                r.nbins = ComputeScottJES(h, min_bins, max_bins);
                r.caption = "binning: Scott-JES, " + std::to_string(r.nbins) + " bins";
                break;
            }
            case Mode::Sturges: {
                r.nbins = ComputeSturges(h, min_bins, max_bins);
                r.caption = "binning: Sturges, " + std::to_string(r.nbins) + " bins";
                break;
            }
            case Mode::PeakAware: {
                double sigma = 0;
                r.nbins = ComputePeakAware(h, xmin, xmax, min_bins, max_bins, sigma);
                r.sigma_peak = sigma;
                if (sigma > 0) {
                    char buf[128];
                    std::snprintf(buf, sizeof(buf),
                                  "binning: peakaware, %d bins (#sigma=%.1f)",
                                  r.nbins, sigma);
                    r.caption = buf;
                } else {
                    r.caption = "binning: peakaware (FD fallback), "
                              + std::to_string(r.nbins) + " bins";
                }
                break;
            }
            case Mode::Variable: {
                r.edges = (channel == "MJB") ? VariableEdgesMJB() : VariableEdgesDB();
                // Clip to [xmin, xmax] in case the channel's range moves.
                while (!r.edges.empty() && r.edges.front() < xmin)
                    r.edges.erase(r.edges.begin());
                while (!r.edges.empty() && r.edges.back()  > xmax)
                    r.edges.pop_back();
                if (r.edges.size() < 2 || r.edges.front() > xmin || r.edges.back() < xmax) {
                    // Ensure the edge list spans the full range.
                    if (r.edges.empty() || r.edges.front() > xmin)
                        r.edges.insert(r.edges.begin(), xmin);
                    if (r.edges.back() < xmax)
                        r.edges.push_back(xmax);
                }
                r.nbins = (int)r.edges.size() - 1;
                r.caption = "binning: variable, " + std::to_string(r.nbins) + " bins";
                break;
            }
            case Mode::FitFriendly: {
                r.nbins = ComputeFitFriendly(h, xmin, xmax, min_bins, max_bins);
                r.caption = "binning: fit-friendly, " + std::to_string(r.nbins) + " bins";
                break;
            }
        }

        // Log line so runs are self-documenting in the terminal output.
        std::cout << "[BinningConfig] channel=" << channel
                  << " mode=" << r.tag
                  << " nbins=" << r.nbins
                  << " range=[" << xmin << "," << xmax << "]";
        if (r.sigma_peak > 0) std::cout << " sigma_peak=" << r.sigma_peak;
        std::cout << "\n";
        return r;
    }

    inline TH1F* MakeHist(const std::string& name, const std::string& title,
                          const Result& r, double xmin, double xmax) {
        if (r.edges.empty()) {
            return new TH1F(name.c_str(), title.c_str(), r.nbins, xmin, xmax);
        }
        return new TH1F(name.c_str(), title.c_str(),
                         (int)r.edges.size() - 1, r.edges.data());
    }

    // ========================================================================
    // RESPONSE HISTOGRAM BINNING 
    // ========================================================================
    // Parallel to Compute() but for response histograms (R_DB in [0, 2])
    // rather than mass (50-200 GeV). Shares the same data-driven algorithms
    // (FD, Scott, ScottJES, Sturges, PeakAware) - only the channel defaults
    // differ.
    //
    // Controlled by EASYBJETS_RESPONSE_BINNING_MODE env var:
    //    "static" (default) — code keeps using cfg.nBinsR (existing behavior)
    //    "dynamic"          — use this function to pick binning at runtime,
    //                          mode chosen by EASYBJETS_BINNING (same as mass)
    //
    // Only used when EASYBJETS_RESPONSE_BINNING_MODE=dynamic; otherwise the
    // macros use cfg.nBinsR / cfg.nBinsTruth. A tag_suffix containing
    // "truth" switches to the truth bin cap.
    // ========================================================================

    constexpr int DB_FIXED_BINS_R  = 40;
    constexpr int MJB_FIXED_BINS_R = 40;

    inline std::vector<double> VariableEdgesResponseDB() {
        std::vector<double> e;
        for (double r = 0.0;  r < 0.80;  r += 0.05) e.push_back(r);   // coarse tail
        for (double r = 0.80; r <= 1.20; r += 0.02) e.push_back(r);   // fine core
        for (double r = 1.25; r <= 2.00; r += 0.05) e.push_back(r);   // coarse tail
        return e;
    }
    inline std::vector<double> VariableEdgesResponseMJB() {
        return VariableEdgesResponseDB();  
    }

    inline bool ResponseDynamicEnabled() {
        static bool cached = [](){
            const char* env = std::getenv("EASYBJETS_RESPONSE_BINNING_MODE");
            if (!env || !*env) return false;
            std::string s(env);
            for (auto& c : s) c = (char)std::tolower((unsigned char)c);
            if (s == "dynamic" || s == "on" || s == "1" || s == "true" || s == "yes")
                return true;
            if (s == "static" || s == "off" || s == "0" || s == "false" || s == "no")
                return false;
            std::cerr << "[BinningConfig] WARN: unknown EASYBJETS_RESPONSE_BINNING_MODE='"
                      << env << "', using static\n";
            return false;
        }();
        return cached;
    }

    inline Result ComputeResponseImpl(TH1F* h, const std::string& channel,
                                       double xmin, double xmax,
                                       int min_bins, int max_bins,
                                       Mode m, const std::string& tag_suffix) {
        Result r;
        r.tag = ModeTag(m);
        if (!tag_suffix.empty()) r.tag += tag_suffix;


        const bool is_truth_pass =
            tag_suffix.find("truth") != std::string::npos;
        if (is_truth_pass) {
            max_bins = RuntimeMaxBinsTruth(max_bins);
        } else {
            max_bins = RuntimeMaxBins(max_bins);
        }

        switch (m) {
            case Mode::Fixed: {
                r.nbins = (channel == "MJB") ? MJB_FIXED_BINS_R : DB_FIXED_BINS_R;
                r.caption = "binning: fixed, " + std::to_string(r.nbins) + " bins";
                break;
            }
            case Mode::FD: {
                r.nbins = ComputeFD(h, xmin, xmax, min_bins, max_bins);
                r.caption = "binning: FD, " + std::to_string(r.nbins) + " bins";
                break;
            }
            case Mode::Scott: {
                r.nbins = ComputeScott(h, xmin, xmax, min_bins, max_bins);
                r.caption = "binning: Scott, " + std::to_string(r.nbins) + " bins";
                break;
            }
            case Mode::ScottJES: {
                r.nbins = ComputeScottJES(h, min_bins, max_bins);
                r.caption = "binning: Scott-JES, " + std::to_string(r.nbins) + " bins";
                break;
            }
            case Mode::Sturges: {
                r.nbins = ComputeSturges(h, min_bins, max_bins);
                r.caption = "binning: Sturges, " + std::to_string(r.nbins) + " bins";
                break;
            }
            case Mode::PeakAware: {
                double sigma = 0;
                double neff = EffectiveN(h);
                if (neff < 20) {
                    std::cerr << "[BinningConfig-Resp] peakaware: too few events ("
                              << neff << "), FD fallback\n";
                    r.nbins = ComputeFD(h, xmin, xmax, min_bins, max_bins);
                    r.caption = "binning: peakaware-FD-fb, " + std::to_string(r.nbins) + " bins";
                    break;
                }
                TF1 g("g_pa_resp", "gaus", 0.85, 1.15);
                g.SetParameter(1, 1.0);
                g.SetParameter(2, 0.1);
                int status = h->Fit(&g, "QNR");
                sigma = g.GetParameter(2);
                if (status != 0 || sigma < 0.02 || sigma > 0.30) {
                    std::cerr << "[BinningConfig-Resp] peakaware: pre-fit failed, "
                                 "FD fallback (sigma=" << sigma << ")\n";
                    r.nbins = ComputeFD(h, xmin, xmax, min_bins, max_bins);
                    r.caption = "binning: peakaware-FD-fb, " + std::to_string(r.nbins) + " bins";
                    break;
                }
                r.sigma_peak = sigma;
                double w = sigma / PEAKAWARE_K;
                int n = static_cast<int>(std::ceil((xmax - xmin) / w));
                r.nbins = std::max(min_bins, std::min(max_bins, n));
                char buf[128];
                std::snprintf(buf, sizeof(buf),
                              "binning: peakaware, %d bins (#sigma=%.3f)",
                              r.nbins, sigma);
                r.caption = buf;
                break;
            }
            case Mode::Variable: {
                r.edges = (channel == "MJB") ? VariableEdgesResponseMJB()
                                              : VariableEdgesResponseDB();
                while (!r.edges.empty() && r.edges.front() < xmin) r.edges.erase(r.edges.begin());
                while (!r.edges.empty() && r.edges.back()  > xmax) r.edges.pop_back();
                if (r.edges.empty() || r.edges.front() > xmin)
                    r.edges.insert(r.edges.begin(), xmin);
                if (r.edges.back() < xmax)
                    r.edges.push_back(xmax);
                r.nbins = (int)r.edges.size() - 1;
                r.caption = "binning: variable, " + std::to_string(r.nbins) + " bins";
                break;
            }
            case Mode::FitFriendly: {
                r.nbins = ComputeFitFriendly(h, xmin, xmax, min_bins, max_bins);
                r.caption = "binning: fit-friendly, " + std::to_string(r.nbins) + " bins";
                break;
            }
        }

        std::cout << "[BinningConfig-Resp] channel=" << channel
                  << " mode=" << r.tag
                  << " nbins=" << r.nbins
                  << " range=[" << xmin << "," << xmax << "]";
        if (r.sigma_peak > 0) std::cout << " sigma_peak=" << r.sigma_peak;
        std::cout << "\n";
        return r;
    }

    inline Result ComputeResponse(TH1F* h, const std::string& channel,
                                   double xmin, double xmax,
                                   int min_bins, int max_bins) {
        return ComputeResponseImpl(h, channel, xmin, xmax,
                                    min_bins, max_bins,
                                    CurrentMode(), "");
    }

    // Same with an explicit mode, used for the truth histograms when their
    // policy asks for a different algorithm than the reco one.
    inline Result ComputeResponseWithMode(TH1F* h, const std::string& channel,
                                           double xmin, double xmax,
                                           int min_bins, int max_bins,
                                           Mode mode_override,
                                           const std::string& tag_suffix = "") {
        return ComputeResponseImpl(h, channel, xmin, xmax,
                                    min_bins, max_bins,
                                    mode_override, tag_suffix);
    }

}  // namespace BinningConfig

#endif  // EASYBJETS_BINNING_CONFIG_H
