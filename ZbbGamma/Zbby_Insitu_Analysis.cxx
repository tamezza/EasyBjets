/**
 * =============================================================================================
 *  Zbby_Insitu_Analysis.cxx - slimming and event selection for Z(->bb)+gamma
 *  in-situ studies using large-R (bb)-tagged jets - Direct Balance response measurements
 * =============================================================================================
 * ===========================================
 * @author:    Thiziri AMEZZA
 * @email:     thiziri.amezza@cern.ch
 * @date:      01 April 2026
 * ===========================================
 *
 * Consolidated analysis script for:
 *   - nominal jet type   (Baseline)
 *   - bJER10v00 jet type (regressed model v00)
 *   - bJER10v01 jet type (regressed model v01)
 *   - MuonCorr jet type  (nominal jet type, muon-corrected dataset)
 * 
 * Reads a raw ntuple, applies the selection below, and writes a slim ntuple
 * with the selected events plus a few new branches (selected_ljet_index,
 * R_DB, R_truth_DB, weight, overlap dR). A second tree, CutMassNominal,
 * holds the events that also pass 50 < m_J < 350 GeV. The cutflow and the
 * flavour composition of the selected jet (MC) are stored as histograms and
 * the cutflow is also written as text under CutFlow/.
 *
 *   root 'Zbby_Insitu_Analysis.cxx("<file>", isData, "<jetType>", "<wp>", isMuonCorr)'
 *
 * Output goes to ../DataSets/ZbbGamma/NtupleSlim_Run{2,3}/ depending on
 * EASYBJETS_RUN, see Insitu_Analysis.h. The Run_DB_*.sh scripts loop over
 * files, jet types and working points.
 *
 * Selection (cutflow bins):
 *    0  all events
 *    1  pT(gamma) > 150 GeV
 *    2  |eta(gamma)| < 1.37
 *    3  leading large-R jet pT > 200 GeV
 *    4  dR(gamma, J) > 0.2
 *    5  |dphi(gamma, J)| > 2.8
 *    6  no small-R jet failing JVT (pT < 60, |eta| < 2.4, jvt <= 0.59)
 *    7  small-R jets within dR < 0.4 of the photon are ignored    (no event veto)
 *    8  small-R jets within dR < 1.4 of the large-R jet are ignored (no event veto)
 *    9  no remaining small-R jet with pT > max(15 GeV, 0.1 pT_ref)
 *   10  selected jet passes the working point
 *   11  50 < m_J < 350 GeV   (CutMassNominal tree only)
 *
* Direct balance: the photon recoils against the (bb)-tagged large-R jet and
 *   R_DB = pT(J) / pT_ref,   pT_ref = pT(gamma) |cos dphi(gamma, J)|
 * measures the jet response. For MC the truth response R_truth_DB uses the
 * matched truth jet pT instead of the reco one.
 *
 * @param fn           Input ROOT file path
 * @param isData       Boolean flag indicating if processing real data (true) or MC simulation (false)
 * @param jetTypeStr   String specifying which jet variables to use: "nominal", "bJER10v00", or "bJER10v01"
 * @param workingPoint Dxbb working point: "0p25", "0p30", "0p37", "0p46",
 *                     "0p58", "0p74", "0p94", "1p25", or "1p55"
 * @param isMuonCorr   true to label output as muon-corrected dataset
 * 
 * =============================================================================================
 * @version: 1.0
 * =============================================================================================
 */

#include "../common/Insitu_Analysis.h"  // shared: JetType, DeltaPhi, WP normalisation, output paths
#include <fstream>
#include <limits>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>


void Zbby_Insitu_Analysis(
    Str fn           = "/data/atlas/tamezza/bJets_Insitu_Analysis/bJets_Samples/Run2/Zbby/GN2x_JetCalibxNtuple_ZbbGamma/data.root",
    bool isData      = true,
    std::string jetTypeStr   = "nominal",
    std::string workingPoint = "0p46",
    bool isMuonCorr  = false)
{
    //==========================================================================
    // Setup
    //==========================================================================

    JetType jetType       = parseJetType(jetTypeStr);
    std::string jetTypeSuffix = getJetTypeSuffix(jetType);

    if (!NormalizeWorkingPoint(workingPoint)) return;

    printTime();

    // Get the main tree
    auto t = getTree(fn, "nominal");
    if (!t) {
        std::cerr << "Error: Cannot get tree from file " << fn << std::endl;
        return;
    }
    int Nall = t->GetEntries();

    // Build output paths (era-aware; set EASYBJETS_RUN to pick Run2/Run3)
    InsituPaths paths = MakeInsituPaths(fn, "ZbbGamma", "DB", jetTypeSuffix, workingPoint, isMuonCorr);
    CreateInsituDirs(paths);
    Str fileName      = paths.baseName;
    Str fn_out        = paths.fn_out;
    Str CutFlowFolder = paths.cutFlowFolder;

    std::cout << "\nInput file:  " << fn  << "\n"
              << "Output file: " << fn_out << "\n"
              << "Jet type:    " << jetTypeSuffix
              << (isMuonCorr ? " (MuonCorr dataset)" : "") << "\n"
              << "Working point: " << workingPoint << "\n\n";

    fileName.ReplaceAll(".root", "");

    //==========================================================================
    // Initialize Histograms and Counters
    //==========================================================================

    // Cutflow tracking
    std::map<std::string, float> cutflow;
    std::map<std::string, int>   eventCounter;
    for (int i = 0; i <= 11; ++i) {
        cutflow["Cut_"       + std::to_string(i)] = 0;
        eventCounter["Cut_" + std::to_string(i)] = 0;
    }

    // Create histograms 
    auto h_cutflow = std::make_unique<TH1D>(
        "h_cutflow", "Cutflow;Cuts;Events", 12, 0, 12);

    // CUT ORDER: 
    const std::vector<std::string> cutNames = {
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
        "D_{xbb} " + workingPoint + "% WP",        
        "M^{J}"                                  
    };

    // Set histogram labels
    for (size_t i = 0; i < cutNames.size(); ++i)
        h_cutflow->GetXaxis()->SetBinLabel(i + 1, cutNames[i].c_str());

    // Background composition histograms (MC only)
    const std::vector<std::string> jetCountLabels = {
        "bb (#geq 2 b-had.)",
        "bc (1 b-had., #geq 1 c-had.)",
        "b (1 b-had., 0 c-had.)",
        "cc (0 b-had., #geq 2 c-had.)",
        "c (0 b-had., 1 c-had.)",
        "light (0 b-had., 0 c-had.)"
    };

    auto h_jet_counts = std::make_unique<TH1D>(
        "h_jet_counts", "Jet Flavour Counts;Category;Events", 6, 0, 6);
    auto h_background_composition = isData ? nullptr : std::make_unique<TH1D>(
        "h_background_composition", "Background Composition (post-b-tag);Category;Events", 6, 0, 6);
    auto h_bgk_comp_Mcut = isData ? nullptr : std::make_unique<TH1D>(
        "h_bgk_comp_Mcut", "Background Composition (post-mass cut);Category;Events", 6, 0, 6);

    if (!isData) {
        for (size_t i = 0; i < jetCountLabels.size(); ++i) {
            h_jet_counts->GetXaxis()->SetBinLabel(i + 1, jetCountLabels[i].c_str());
            h_background_composition->GetXaxis()->SetBinLabel(i + 1, jetCountLabels[i].c_str());
            h_bgk_comp_Mcut->GetXaxis()->SetBinLabel(i + 1, jetCountLabels[i].c_str());
        }
    }

    //==========================================================================
    // Variable Declarations
    //==========================================================================

    // Event-level quantities
    double weight       = 0;
    double total_weight = 0;
    int    Npass        = 0;

    // Photon kinematics
    ROOT::RVecF *ph_pt  = nullptr;
    ROOT::RVecF *ph_eta = nullptr;
    ROOT::RVecF *ph_phi = nullptr;
    ROOT::RVecF *ph_e   = nullptr;

    // Large-R jet kinematics - all variants
    ROOT::RVecF *ljet_pt            = nullptr;
    ROOT::RVecF *ljet_eta           = nullptr;
    ROOT::RVecF *ljet_phi           = nullptr;
    ROOT::RVecF *ljet_e             = nullptr;
    ROOT::RVecF *ljet_m             = nullptr;
    ROOT::RVecF *ljet_bJR10v00_pt   = nullptr;
    ROOT::RVecF *ljet_bJR10v00_mass = nullptr;
    ROOT::RVecF *ljet_bJR10v01_pt   = nullptr;
    ROOT::RVecF *ljet_bJR10v01_mass = nullptr;

    // pT and mass of the active jet type. The lambdas read the pointer at
    // call time, so they see the vector ROOT filled for the current entry.
    std::function<ROOT::RVecF*()> getCurrentPt;
    std::function<ROOT::RVecF*()> getCurrentMass;

    // Small-R jet kinematics
    ROOT::RVecF *jet_pt  = nullptr;
    ROOT::RVecF *jet_eta = nullptr;
    ROOT::RVecF *jet_phi = nullptr;
    ROOT::RVecF *jet_e   = nullptr;
    ROOT::RVecF *jet_m   = nullptr;
    ROOT::RVecF *jet_jvt = nullptr;

    // Working point flags - (bb)-tagging
    ROOT::RVecF *dhbb = nullptr;          // GN2x (bb)-tagging score
    ROOT::RVecB *pass_FlatMassQCDEff_0p25 = nullptr;
    ROOT::RVecB *pass_FlatMassQCDEff_0p3  = nullptr;
    ROOT::RVecB *pass_FlatMassQCDEff_0p37 = nullptr;
    ROOT::RVecB *pass_FlatMassQCDEff_0p46 = nullptr;
    ROOT::RVecB *pass_FlatMassQCDEff_0p58 = nullptr;
    ROOT::RVecB *pass_FlatMassQCDEff_0p74 = nullptr;
    ROOT::RVecB *pass_FlatMassQCDEff_0p94 = nullptr;
    ROOT::RVecB *pass_FlatMassQCDEff_1p25 = nullptr;
    ROOT::RVecB *pass_FlatMassQCDEff_1p55 = nullptr;

    // Truth matching branches (MC only - needed for b-JER computation)
    ROOT::RVecI *ljet_truth_label        = nullptr;
    ROOT::RVecI *jet_parton_truth_label  = nullptr;
    //ROOT::RVecF *jet_eff_jvt             = nullptr;
    float jet_eff_jvt = 0;
    ROOT::RVecF *matched_tjet_pt         = nullptr;
    ROOT::RVecF *matched_tjet_m          = nullptr;
    ROOT::RVecF *matched_tjet_eta        = nullptr;
    ROOT::RVecF *matched_tjet_phi        = nullptr;
    ROOT::RVecF *matched_tjet_dR         = nullptr;

    // Background composition branches (MC only)
    ROOT::RVecI *ljet_n_bhadrons = nullptr; // Number of b-hadrons ghost-associated to large-R jet
    ROOT::RVecI *ljet_n_chadrons = nullptr; // Number of c-hadrons ghost-associated to large-R jet

    // Angular and balance variables
    int    selected_ljet_index = -1;

    double delta_phi          = 0;
    double delta_eta_ph_ljet  = 0;
    double delta_R_ph_ljet    = 0;
    double delta_phi_ph_jet   = 0;
    double delta_eta_ph_jet   = 0;
    double delta_R_ph_jet     = 0;
    double delta_phi_ljet_jet = 0;
    double delta_eta_ljet_jet = 0;
    double delta_R_ljet_jet   = 0;
    double delta_eta_overlap  = 0;
    double delta_phi_overlap  = 0;
    double delta_R_overlap    = 0;
    double pt_ref_DB          = 0;
    double R_DB               = 0;
    double R_truth_DB         = 0;

    //==========================================================================
    // Set Up Input Tree Branches
    //==========================================================================

    // Photon
    t->SetBranchAddress("ph_pt",  &ph_pt);
    t->SetBranchAddress("ph_eta", &ph_eta);
    t->SetBranchAddress("ph_phi", &ph_phi);
    t->SetBranchAddress("ph_e",   &ph_e);

    // Common large-R jet branches
    t->SetBranchAddress("ljet_pt",  &ljet_pt);
    t->SetBranchAddress("ljet_eta", &ljet_eta);
    t->SetBranchAddress("ljet_phi", &ljet_phi);
    t->SetBranchAddress("ljet_e",   &ljet_e);
    t->SetBranchAddress("ljet_m",   &ljet_m);

    // Jet-type-specific branches
    switch (jetType) {
        case JetType::BJER10V00:
            if (t->SetBranchAddress("ljet_bJR10v00_pt",   &ljet_bJR10v00_pt)   < 0 ||
                t->SetBranchAddress("ljet_bJR10v00_mass", &ljet_bJR10v00_mass) < 0) {
                std::cerr << "Error: Cannot find bJER10v00 branches in tree\n";
                return;
            }
            getCurrentPt   = [&]() { return ljet_bJR10v00_pt; };
            getCurrentMass = [&]() { return ljet_bJR10v00_mass; };
            std::cout << "Using bJER10v00 jet variables\n";
            break;
        case JetType::BJER10V01:
            if (t->SetBranchAddress("ljet_bJR10v01_pt",   &ljet_bJR10v01_pt)   < 0 ||
                t->SetBranchAddress("ljet_bJR10v01_mass", &ljet_bJR10v01_mass) < 0) {
                std::cerr << "Error: Cannot find bJER10v01 branches in tree\n";
                return;
            }
            getCurrentPt   = [&]() { return ljet_bJR10v01_pt; };
            getCurrentMass = [&]() { return ljet_bJR10v01_mass; };
            std::cout << "Using bJER10v01 jet variables\n";
            break;
        default: // NOMINAL
            getCurrentPt   = [&]() { return ljet_pt; };
            getCurrentMass = [&]() { return ljet_m; };
            std::cout << "Using nominal jet variables\n";
            break;
    }

    // Small-R jets
    t->SetBranchAddress("jet_pt",  &jet_pt);
    t->SetBranchAddress("jet_eta", &jet_eta);
    t->SetBranchAddress("jet_phi", &jet_phi);
    t->SetBranchAddress("jet_e",   &jet_e);
    t->SetBranchAddress("jet_m",   &jet_m);
    t->SetBranchAddress("jet_jvt", &jet_jvt);

    // (bb)-tagging
    t->SetBranchAddress("dhbb", &dhbb);
    if (t->SetBranchAddress("pass_FlatMassQCDEff_0p25", &pass_FlatMassQCDEff_0p25) < 0 ||
        t->SetBranchAddress("pass_FlatMassQCDEff_0p3",  &pass_FlatMassQCDEff_0p3)  < 0 ||
        t->SetBranchAddress("pass_FlatMassQCDEff_0p37", &pass_FlatMassQCDEff_0p37) < 0 ||
        t->SetBranchAddress("pass_FlatMassQCDEff_0p46", &pass_FlatMassQCDEff_0p46) < 0 ||
        t->SetBranchAddress("pass_FlatMassQCDEff_0p58", &pass_FlatMassQCDEff_0p58) < 0 ||
        t->SetBranchAddress("pass_FlatMassQCDEff_0p74", &pass_FlatMassQCDEff_0p74) < 0 ||
        t->SetBranchAddress("pass_FlatMassQCDEff_0p94", &pass_FlatMassQCDEff_0p94) < 0 ||
        t->SetBranchAddress("pass_FlatMassQCDEff_1p25", &pass_FlatMassQCDEff_1p25) < 0 ||
        t->SetBranchAddress("pass_FlatMassQCDEff_1p55", &pass_FlatMassQCDEff_1p55) < 0) {
        std::cerr << "Error: Cannot find working point branches in tree\n";
        return;
    }

    // Event weight
    t->SetBranchAddress("total_weight", &total_weight);

    // Truth matching branches - MC only (needed for b-JER computation)
    if (!isData) {
        if (t->SetBranchAddress("ljet_truth_label",       &ljet_truth_label)       < 0 ||
            t->SetBranchAddress("jet_parton_truth_label", &jet_parton_truth_label) < 0 ||
            t->SetBranchAddress("jet_eff_jvt",            &jet_eff_jvt)            < 0 ||
            t->SetBranchAddress("matched_tjet_pt",        &matched_tjet_pt)        < 0 ||
            t->SetBranchAddress("matched_tjet_m",         &matched_tjet_m)         < 0 ||
            t->SetBranchAddress("matched_tjet_eta",       &matched_tjet_eta)       < 0 ||
            t->SetBranchAddress("matched_tjet_phi",       &matched_tjet_phi)       < 0 ||
            t->SetBranchAddress("matched_tjet_dR",        &matched_tjet_dR)        < 0) {
            std::cerr << "Error: Cannot find truth matching branches in tree\n";
            return;
        }
        if (t->SetBranchAddress("ljet_n_bhadrons", &ljet_n_bhadrons) < 0 ||
            t->SetBranchAddress("ljet_n_chadrons", &ljet_n_chadrons) < 0) {
            std::cerr << "Error: Cannot find ljet_n_bhadrons/ljet_n_chadrons branches in tree\n";
            return;
        }
        std::cout << "Truth matching branches loaded (MC mode)\n";
    }

    //==========================================================================
    // Working Point Selection (done once at setup)
    //==========================================================================
    std::function<ROOT::RVecB*()> getActiveWPVector;

    if      (workingPoint == "0p25") getActiveWPVector = [&]() { return pass_FlatMassQCDEff_0p25; };
    else if (workingPoint == "0p30") getActiveWPVector = [&]() { return pass_FlatMassQCDEff_0p3; };
    else if (workingPoint == "0p37") getActiveWPVector = [&]() { return pass_FlatMassQCDEff_0p37; };
    else if (workingPoint == "0p46") getActiveWPVector = [&]() { return pass_FlatMassQCDEff_0p46; };
    else if (workingPoint == "0p58") getActiveWPVector = [&]() { return pass_FlatMassQCDEff_0p58; };
    else if (workingPoint == "0p74") getActiveWPVector = [&]() { return pass_FlatMassQCDEff_0p74; };
    else if (workingPoint == "0p94") getActiveWPVector = [&]() { return pass_FlatMassQCDEff_0p94; };
    else if (workingPoint == "1p25") getActiveWPVector = [&]() { return pass_FlatMassQCDEff_1p25; };
    else if (workingPoint == "1p55") getActiveWPVector = [&]() { return pass_FlatMassQCDEff_1p55; };
    else {
        std::cerr << "Error: Unknown working point '" << workingPoint << "'\n";
        return;
    }

    // Lambda captures by reference - always gets current pointer value after GetEntry()
    auto getWorkingPointFlag = [&](size_t jetIndex) -> bool {
        ROOT::RVecB* wpVec = getActiveWPVector();
        return (wpVec && jetIndex < wpVec->size()) ? wpVec->at(jetIndex) : false;
    };

    //==========================================================================
    // Output Trees
    //==========================================================================

    // Create output file and trees
    TFile *of = new TFile(fn_out, "RECREATE");
    if (!of || of->IsZombie()) {
        std::cerr << "Error: Cannot create output file " << fn_out << "\n";
        return;
    }
    of->cd();
    auto *ot = t->CloneTree(0);
    
    ot->SetDirectory(of);
    ot->Branch("selected_ljet_index", &selected_ljet_index, "selected_ljet_index/I");
    ot->Branch("R_DB",              &R_DB,              "R_DB/D");
    ot->Branch("R_truth_DB",        &R_truth_DB,        "R_truth_DB/D");
    ot->Branch("weight",            &weight,            "weight/D");
    ot->Branch("delta_eta_overlap", &delta_eta_overlap, "delta_eta_overlap/D");
    ot->Branch("delta_phi_overlap", &delta_phi_overlap, "delta_phi_overlap/D");
    ot->Branch("delta_R_overlap",   &delta_R_overlap,   "delta_R_overlap/D");

    // Secondary tree for events also passing the mass window
    TTree *CutMassNominal = t->CloneTree(0);
    CutMassNominal->SetName("CutMassNominal");
    CutMassNominal->SetDirectory(of);
    CutMassNominal->Branch("selected_ljet_index", &selected_ljet_index, "selected_ljet_index/I");
    CutMassNominal->Branch("R_DB",       &R_DB,       "R_DB/D");
    CutMassNominal->Branch("R_truth_DB", &R_truth_DB, "R_truth_DB/D");
    CutMassNominal->Branch("weight", &weight, "weight/D");

    //==========================================================================
    // Event Loop
    //==========================================================================

    std::cout << std::fixed << std::setprecision(1)
              << "Will loop over " << 1e-6 * Nall << "M events\n";

    for (int i = 0; i < Nall; ++i) {
        // Progress reporting
        if (i && i % 500000 == 0) {
            std::cout << std::setw(6) << 1e-6 * i << "M / "
                      << std::setw(6) << 1e-6 * Nall << "M evts.  "
                      << std::setprecision(1) << 100.0 * i / Nall << "%\n";
            printTime(false);
        }

        //----------------------------------------------------------------------
        // Reset per-event output variables before GetEntry()
        //----------------------------------------------------------------------
        selected_ljet_index = -1;
        delta_phi           = 0;
        delta_eta_ph_ljet   = 0;
        delta_R_ph_ljet     = 0;
        delta_phi_ph_jet    = 0;
        delta_eta_ph_jet    = 0;
        delta_R_ph_jet      = 0;
        delta_phi_ljet_jet  = 0;
        delta_eta_ljet_jet  = 0;
        delta_R_ljet_jet    = 0;
        delta_eta_overlap   = 0;
        delta_phi_overlap   = 0;
        delta_R_overlap     = 0;
        pt_ref_DB           = 0;
        R_DB                = 0;
        R_truth_DB          = 0;

        t->GetEntry(i);
        weight = isData ? 1.0 : total_weight;

        //----------------------------------------------------------------------
        // Cut 0: Initial count
        //----------------------------------------------------------------------
        cutflow["Cut_0"] += weight;
        eventCounter["Cut_0"]++;

        //----------------------------------------------------------------------
        // Cut 1: Photon pT > 150 GeV
        //----------------------------------------------------------------------
        if (ph_pt->size() < 1 || ph_pt->at(0) < 150) continue;
        cutflow["Cut_1"] += weight;
        eventCounter["Cut_1"]++;

        //----------------------------------------------------------------------
        // Cut 2: Photon |eta| < 1.37
        //----------------------------------------------------------------------
        if (ph_eta->size() < 1 || std::abs(ph_eta->at(0)) >= 1.37) continue;
        cutflow["Cut_2"] += weight;
        eventCounter["Cut_2"]++;

        //----------------------------------------------------------------------
        // Cut 3: Select highest-pT large-R jet with pT > 200 GeV
        //        ((bb)-tagging comes later at Cut 10)
        //----------------------------------------------------------------------
        if (getCurrentPt()->size() < 1) continue;

        // Find highest-pT large-R jet (no b-tag requirement yet)
        int index_selected = -1;
        double max_pt = -1.0;
        for (size_t j = 0; j < getCurrentPt()->size(); ++j) {
            double pt_j = getCurrentPt()->at(j);
            if (pt_j > max_pt) {
                max_pt = pt_j;
                index_selected = static_cast<int>(j);
            }
        }

        if (index_selected < 0 || max_pt < 200) continue;
        cutflow["Cut_3"] += weight;
        eventCounter["Cut_3"]++;

        selected_ljet_index = index_selected;

        // Ensure we have valid eta/phi for the selected large-R jet
        if (ljet_eta->size() <= (size_t)index_selected ||
            ljet_phi->size() <= (size_t)index_selected) continue;
        if (ph_phi->size() < 1) continue;

        // Calculate photon – selected large-R jet angular separation
        delta_phi         = DeltaPhi(ph_phi->at(0), ljet_phi->at(index_selected));
        delta_eta_ph_ljet = std::abs(ph_eta->at(0) - ljet_eta->at(index_selected));
        delta_R_ph_ljet   = std::sqrt(delta_phi * delta_phi +
                                      delta_eta_ph_ljet * delta_eta_ph_ljet);

        //----------------------------------------------------------------------
        // Cut 4: dR(gamma, J) > 0.2  -  photon–large-R jet overlap removal
        //----------------------------------------------------------------------
        if (delta_R_ph_ljet <= 0.2) continue;
        cutflow["Cut_4"] += weight;
        eventCounter["Cut_4"]++;

        //----------------------------------------------------------------------
        // Cut 5: |dphi(gamma, J)| > 2.8  -  back-to-back topology
        //----------------------------------------------------------------------
        if (std::abs(delta_phi) < 2.8) continue;
        cutflow["Cut_5"] += weight;
        eventCounter["Cut_5"]++;

        // Reference pT for Direct Balance: projected photon pT along the selected jet axis
        pt_ref_DB = ph_pt->at(0) * std::abs(std::cos(delta_phi));

        //----------------------------------------------------------------------
        // Cuts 6–9: Small-R jet selection
        //----------------------------------------------------------------------
        bool passJetCut6 = true;   // JVT pile-up removal
        bool passJetCut9 = true;   // Extra-radiation veto

        delta_eta_overlap = 0.0;
        delta_phi_overlap = 0.0;
        delta_R_overlap   = -1.0;

        double min_delta_R_overlap = std::numeric_limits<double>::max();

        if (jet_pt->size() > 0) {
            for (size_t j = 0; j < jet_pt->size(); ++j) {

                // dR between small-R jet and selected large-R jet
                double deta_overlap = std::abs(jet_eta->at(j) - ljet_eta->at(index_selected));
                double dphi_overlap = DeltaPhi(jet_phi->at(j), ljet_phi->at(index_selected));
                double dR_overlap   = std::sqrt(deta_overlap * deta_overlap +
                                                dphi_overlap * dphi_overlap);

                if (dR_overlap < min_delta_R_overlap) {
                    min_delta_R_overlap = dR_overlap;
                    delta_eta_overlap   = deta_overlap;
                    delta_phi_overlap   = dphi_overlap;
                    delta_R_overlap     = dR_overlap;
                }

                // Cut 6: JVT - reject pile-up jets in tracker acceptance
                if (jet_pt->at(j) < 60 &&
                    jet_jvt->at(j) <= 0.59 &&
                    std::abs(jet_eta->at(j)) < 2.4) {
                    passJetCut6 = false;
                    break;
                }

                // Cut 7: skip photon-ghost jets - overlap removal, not an event veto
                delta_phi_ph_jet = DeltaPhi(ph_phi->at(0), jet_phi->at(j));
                delta_eta_ph_jet = std::abs(ph_eta->at(0) - jet_eta->at(j));
                delta_R_ph_jet   = std::sqrt(delta_phi_ph_jet * delta_phi_ph_jet +
                                             delta_eta_ph_jet * delta_eta_ph_jet);
                if (delta_R_ph_jet <= 0.4) continue;  // exclude jet, keep event

                // Cut 8: skip jets inside large-R cone - overlap removal, not an event veto
                delta_phi_ljet_jet = DeltaPhi(ljet_phi->at(index_selected), jet_phi->at(j));
                delta_eta_ljet_jet = std::abs(ljet_eta->at(index_selected) - jet_eta->at(j));
                delta_R_ljet_jet   = std::sqrt(delta_phi_ljet_jet * delta_phi_ljet_jet +
                                               delta_eta_ljet_jet * delta_eta_ljet_jet);
                if (delta_R_ljet_jet <= 1.4) continue;  // exclude jet, keep event

                // Cut 9: extra-radiation veto - any surviving jet above threshold kills event
                if (jet_pt->at(j) >= std::max(15.0, 0.1 * pt_ref_DB)) {
                    passJetCut9 = false;
                    break;
                }
            }

            // Apply event-level decisions
            if (!passJetCut6) continue;
            cutflow["Cut_6"] += weight; eventCounter["Cut_6"]++;

            // Cuts 7 & 8 are jet-level overlap removals - all events reaching here pass them
            cutflow["Cut_7"] += weight; eventCounter["Cut_7"]++;
            cutflow["Cut_8"] += weight; eventCounter["Cut_8"]++;

            if (!passJetCut9) continue;
            cutflow["Cut_9"] += weight; eventCounter["Cut_9"]++;

        } else {
            // No small-R jets: event passes cuts 6–9 automatically
            for (int c = 6; c <= 9; ++c) {
                std::string key = "Cut_" + std::to_string(c);
                cutflow[key]      += weight;
                eventCounter[key]++;
            }
        }

        //----------------------------------------------------------------------
        // Cut 10: D_xbb (bb)-tagging working point on the selected jet
        //         Check if the SELECTED jet is (bb)-tagged at the current WP
        //----------------------------------------------------------------------
        if (!getWorkingPointFlag(index_selected)) continue;
        cutflow["Cut_10"] += weight;
        eventCounter["Cut_10"]++;

        // Background composition for the selected (bb)-tagged large-R jet (MC only)
        if (!isData &&
            (size_t)index_selected < ljet_n_bhadrons->size() &&
            (size_t)index_selected < ljet_n_chadrons->size()) {
            int nb = ljet_n_bhadrons->at(index_selected);
            int nc = ljet_n_chadrons->at(index_selected);
            int jet_bin = (nb >= 2)                ? 0 :  // bb
                          (nb == 1 && nc >= 1)     ? 1 :  // bc
                          (nb == 1 && nc == 0)     ? 2 :  // b
                          (nb == 0 && nc >= 2)     ? 3 :  // cc
                          (nb == 0 && nc == 1)     ? 4 :  // c
                                                     5;   // light
            h_jet_counts->Fill(jet_bin, weight);
            h_background_composition->Fill(jet_bin, weight);
        }

        //----------------------------------------------------------------------
        // Calculate R_DB (Direct Balance response)
        //----------------------------------------------------------------------
        R_DB = (pt_ref_DB > 0) ? getCurrentPt()->at(index_selected) / pt_ref_DB : 0;

        // Compute truth response for b-JER: p_T^{truth} / p_T^{ref,db}
        R_truth_DB = 0;
        if (!isData && matched_tjet_pt &&
            (int)matched_tjet_pt->size() > index_selected && pt_ref_DB > 0) {
            R_truth_DB = matched_tjet_pt->at(index_selected) / pt_ref_DB;
        }

        // Save event to main output tree
        ot->Fill();
        ++Npass;

        //----------------------------------------------------------------------
        // Cut 11: Mass window  50 < M(J) < 350 GeV
        //----------------------------------------------------------------------
        if (getCurrentMass()->size() <= (size_t)index_selected) continue;
        double ljet_mass_val = getCurrentMass()->at(index_selected);
        if (ljet_mass_val > 50 && ljet_mass_val < 350) {
            CutMassNominal->Fill();
            cutflow["Cut_11"] += weight;
            eventCounter["Cut_11"]++;
            if (!isData &&
                (size_t)index_selected < ljet_n_bhadrons->size() &&
                (size_t)index_selected < ljet_n_chadrons->size()) {
                int nb = ljet_n_bhadrons->at(index_selected);
                int nc = ljet_n_chadrons->at(index_selected);
                int jet_bin = (nb >= 2)                ? 0 :
                              (nb == 1 && nc >= 1)     ? 1 :
                              (nb == 1 && nc == 0)     ? 2 :
                              (nb == 0 && nc >= 2)     ? 3 :
                              (nb == 0 && nc == 1)     ? 4 :
                                                         5;
                h_bgk_comp_Mcut->Fill(jet_bin, weight);
            }
        }
    }

    //==========================================================================
    // Write Cutflow
    //==========================================================================

    {
        Str cutflowFileName = CutFlowFolder + Str(fileName) + "_cutflow_results.txt";
        std::ofstream outfile(cutflowFileName.Data());
        if (!outfile.is_open()) {
            std::cerr << "Warning: Cannot create cutflow file " << cutflowFileName << "\n";
        } else {
            outfile << "Cutflow results - jet type: " << jetTypeSuffix
                    << (isMuonCorr ? " (MuonCorr)" : "")
                    << ", working point: " << workingPoint << "\n";

            for (int i = 0; i <= 11; ++i) {
                std::string key = "Cut_" + std::to_string(i);
                float absEff = (i == 0) ? 1.0f
                    : cutflow[key] / cutflow["Cut_0"];
                float relEff = 1.0f;
                if (i > 0) {
                    std::string prevKey = "Cut_" + std::to_string(i - 1);
                    relEff = (cutflow[prevKey] > 0)
                        ? cutflow[key] / cutflow[prevKey]
                        : 0.0f;
                }
                auto line = key + " " + cutNames[i] + ": "
                    + std::to_string(cutflow[key])
                    + " (Events: " + std::to_string(eventCounter[key])
                    + ", Abs Eff: " + std::to_string(absEff * 100) + "%"
                    + ", Rel Eff: " + std::to_string(relEff * 100) + "%)\n";

                outfile    << line;
                std::cout  << line;
            }
            outfile.close();
        }
    }

    //==========================================================================
    // Write Histograms and Cutflow Plot
    //==========================================================================

    for (int i = 0; i <= 11; ++i)
        h_cutflow->SetBinContent(i + 1, cutflow["Cut_" + std::to_string(i)]);

    {
        TCanvas c1("c1", "Cutflow Histogram", 800, 600);
        c1.SetLogy();
        h_cutflow->Draw();

        TLatex latex;
        latex.SetNDC();
        latex.SetTextSize(0.03);
        latex.DrawLatex(0.15, 0.92,
            "Z(#rightarrowbb)+#gamma - Photon (#gamma), Large-R Jet (J) & Small-R Jet (j)");
        c1.SaveAs((CutFlowFolder + Str(fileName) + "_cutflow.png").Data());
    }

    //==========================================================================
    // Write Histograms to Output File
    //==========================================================================

    of->cd();
    h_cutflow->Write();
    if (!isData) {
        h_jet_counts->Write();
        h_background_composition->Write();
        h_bgk_comp_Mcut->Write();
    }

    //==========================================================================
    // Cleanup and Final Statistics
    //==========================================================================

    // Print final statistics
    std::cout << std::fixed << std::setprecision(1)
              << "\nProcessing completed - jet type: " << jetTypeSuffix
              << (isMuonCorr ? " (MuonCorr)" : "")
              << ", working point: " << workingPoint << "\n"
              << "Events in:  " << Nall  << "\n"
              << "Events out: " << Npass << "  ("
              << std::setprecision(2) << 100.0 * Npass / Nall << "%)\n";

    printTime(false);

    // Write and close output file
    of->cd();
    ot->Write();
    CutMassNominal->Write();
    of->Write();

    std::cout << "Output written to: " << fn_out << "\n";

    // Clean up - remove from ROOT's global list before deleting
    gROOT->GetListOfFiles()->Remove(of);
    of->Close();
    delete of;
    of = nullptr;  

    printTime();
}


//==============================================================================
// Convenience wrapper functions
//==============================================================================

void Zbby_Nominal_Slim(
    Str fn = "../../JetCalibNtuples/ZbbGamma_bJER2/data.root",
    bool isData = true, std::string workingPoint = "0p46")
{
    Zbby_Insitu_Analysis(fn, isData, "nominal", workingPoint, false);
}

void Zbby_bJER10v00_Slim(
    Str fn = "../../JetCalibNtuples/ZbbGamma_bJER2/data.root",
    bool isData = true, std::string workingPoint = "0p46")
{
    Zbby_Insitu_Analysis(fn, isData, "bJER10v00", workingPoint, false);
}

void Zbby_bJER10v01_Slim(
    Str fn = "../../JetCalibNtuples/ZbbGamma_bJER2/data.root",
    bool isData = true, std::string workingPoint = "0p46")
{
    Zbby_Insitu_Analysis(fn, isData, "bJER10v01", workingPoint, false);
}

void Zbby_MuonCorr_Slim(
    Str fn = "../../JetCalibNtuples/ZbbGamma_muon_corr/data.root",
    bool isData = true, std::string workingPoint = "0p46")
{
    Zbby_Insitu_Analysis(fn, isData, "nominal", workingPoint, true);
}

// Run all three jet types for a single working point
void Zbby_All_JetTypes_SingleWP(
    Str fn = "../../JetCalibNtuples/ZbbGamma_bJER2/data.root",
    bool isData = true, std::string workingPoint = "0p46")
{
    std::cout << "\n=== Running all jet types (WP: " << workingPoint << ") ===\n";
    Zbby_Insitu_Analysis(fn, isData, "nominal",   workingPoint, false);
    Zbby_Insitu_Analysis(fn, isData, "bJER10v00", workingPoint, false);
    Zbby_Insitu_Analysis(fn, isData, "bJER10v01", workingPoint, false);
    std::cout << "\n=== All jet type analyses completed ===\n";
}

// Run all jet types for all working points
void Zbby_All_Analyses(
    Str fn = "../../JetCalibNtuples/ZbbGamma_bJER2/data.root",
    bool isData = true)
{
    const std::vector<std::string> workingPoints =
        {"0p25","0p30","0p37","0p46","0p58","0p74","0p94","1p25","1p55"};
    const std::vector<std::string> jetTypes =
        {"nominal","bJER10v00","bJER10v01"};

    int total = workingPoints.size() * jetTypes.size(), n = 0;
    std::cout << "\n=== Running all jet types and working points ===\n";

    for (const auto& wp : workingPoints) {
        for (const auto& jt : jetTypes) {
            std::cout << "\n--- Analysis " << ++n << "/" << total
                      << ": " << jt << " (WP: " << wp << ") ---\n";
            Zbby_Insitu_Analysis(fn, isData, jt, wp, false);
        }
    }
    std::cout << "\n=== All analyses completed ===\n";
}

// Run selected working points only
void Zbby_Selected_WorkingPoints(
    Str fn = "../../JetCalibNtuples/ZbbGamma_bJER2/data.root",
    bool isData = true,
    std::vector<std::string> selectedWPs = {"0p37","0p46","0p58"})
{
    const std::vector<std::string> jetTypes =
        {"nominal","bJER10v00","bJER10v01"};

    int total = selectedWPs.size() * jetTypes.size(), n = 0;
    std::cout << "\n=== Running selected working points ===\n";

    for (const auto& wp : selectedWPs) {
        for (const auto& jt : jetTypes) {
            std::cout << "\n--- Analysis " << ++n << "/" << total
                      << ": " << jt << " (WP: " << wp << ") ---\n";
            Zbby_Insitu_Analysis(fn, isData, jt, wp, false);
        }
    }
    std::cout << "\n=== Selected WP analyses completed ===\n";
}