/**
 * =============================================================================================
 *  ZbbJets_Insitu_Analysis.cxx - slimming and event selection for Z(->bb)+jets
 *  in-situ studies using large-R jets for Multijets Balance response measurements
 * =============================================================================================
 * ===========================================
 * @author:    Thiziri AMEZZA
 * @email:     thiziri.amezza@cern.ch
 * @date:      21 January 2025
 * ===========================================
 *
 * Consolidated analysis script for:
 *   - nominal jet type   (Baseline)
 *   - bJER10v00 jet type (regressed model v00)
 *   - bJER10v01 jet type (regressed model v01)
 *   - MuonCorr jet type  (nominal jet type, muon-corrected dataset)
 *
 * Reads a raw ntuple, applies the selection below and writes a slim ntuple
 * with the selected events plus the balance variables (index_tagged, R_MJB,
 * R_truth_MJB, pt_recoil, alpha, weight, ...). Trees named Pre<Cut> hold the
 * events after each intermediate cut, for cutflow checks. The cutflow and the
 * flavour composition of the tagged jet (MC) are stored as histograms, the
 * cutflow also as text under CutFlow/.
 *
 *   root 'ZbbJets_Insitu_Analysis.cxx("<file>", isData, "<jetType>", "<wp>", isMuonCorr)'
 *
 * Output goes to ../DataSets/ZbbJets/NtupleSlim_Run{2,3}/ depending on
 * EASYBJETS_RUN, see Insitu_Analysis.h. The Run_MJB_*.sh scripts loop over
 * files, jet types and working points.
 *
 * Selection (cutflow bins). MC events with total_weight > 100 are dropped
 * before the first bin, and events with more than two large-R jets after it.
 *    0  all events
 *    1  a large-R jet passes the working point (the highest-pT one is used)
 *    2  pT(J) > 450 GeV
 *    3  |eta(J)| < 1.2
 *    4  at least one small-R jet with pT > 25 GeV, |eta| < 2.8, dR(j, J) > 1
 *    5  none of those fails JVT (pT < 60, |eta| < 2.4, jvt <= 0.59)
 *    6  |alpha - pi| < 0.3, alpha = dphi(J, recoil)
 *    7  dR(j, J) > 1.5 for every selected small-R jet
 * 
 * Multijet balance: the (bb)-tagged large-R jet recoils against the system of
 * small-R jets and
 *   R_MJB = pT(J) / pT_recoil
 * measures the jet response, pT_recoil being the vector sum of the selected
 * small-R jets (the single jet pT when there is only one). For MC the truth
 * response R_truth_MJB uses the matched truth jet pT.
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

#include "../common/Insitu_Analysis.h"  
#include <fstream>
#include <limits>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>


void ZbbJets_Insitu_Analysis(
    Str fn = "/data/atlas/tamezza/BJetCalib_Run3/JetCalibProcess/Zbbj/Zbbj_xbbCalib_19Jan2026/multijets_801165.root", 
    bool isData = false, 
    std::string jetTypeStr = "nominal",
    std::string workingPoint = "0p46",
    bool isMuonCorr  = false) 
{
    //==========================================================================
    // Setup
    //==========================================================================

    JetType jetType = parseJetType(jetTypeStr);
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
    InsituPaths paths = MakeInsituPaths(fn, "ZbbJets", "MJB", jetTypeSuffix, workingPoint, isMuonCorr);
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
    for (int i = 0; i <= 7; ++i) {
        cutflow["Cut_" + std::to_string(i)] = 0;
        eventCounter["Cut_" + std::to_string(i)] = 0;
    }

    // Create histograms 
    auto h_cutflow = std::make_unique<TH1D>("h_cutflow", "Event Selection;Cut;Events", 8, 0, 8);
    auto h_n_large_r_jets = std::make_unique<TH1D>("h_n_large_r_jets", 
        "Number of Large-R Jets;N_{large-R jets};Events", 7, 0, 7);
    auto h_n_recoil_jets = std::make_unique<TH1D>("h_n_recoil_jets", 
        "Number of Small-R Jets;N_{recoil jets};Events", 7, 0, 7);

    // Post-selection histograms
    auto h_n_ljets = std::make_unique<TH1D>("h_n_ljets", 
        "Final Large-R Jet Multiplicity;N_{large-R jets};Events", 10, 0, 10);
    auto h_n_jets = std::make_unique<TH1D>("h_n_jets", 
        "Final Small-R Jet Multiplicity;N_{small-R jets};Events", 10, 0, 10);

    // Define labels for cuts and categories
    std::vector<std::string> cutNames = {
        "Initial Events",                               
        "D_{xbb} " + workingPoint + "% WP",            
        "p_{T}^{J} > 450",                             
        "|#eta^{J}| < 1.2",                             
        "Small-R p_{T}^{J} > 25  & #eta^{j} < 2.8",     
        "JVT > 0.59",                                   
        "|#alpha - #pi| < 0.3",                         
        "#Delta R (J1, MJB) > 1.5",                     
    };

    // Set histogram labels
    for (size_t i = 0; i < cutNames.size(); ++i) {
        h_cutflow->GetXaxis()->SetBinLabel(i + 1, cutNames[i].c_str());
    }

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

    if (!isData) {
        for (size_t i = 0; i < jetCountLabels.size(); ++i) {
            h_jet_counts->GetXaxis()->SetBinLabel(i + 1, jetCountLabels[i].c_str());
            h_background_composition->GetXaxis()->SetBinLabel(i + 1, jetCountLabels[i].c_str());
        }
    }

    //==========================================================================
    // Variable Declarations
    //==========================================================================
    
    // Event-level quantities
    double weight = 0; 
    double total_weight = 0;
    
    // Jet kinematics - declare all possible variables
    ROOT::RVecF *ljet_pt = nullptr, *ljet_eta = nullptr, *ljet_phi = nullptr;
    ROOT::RVecF *ljet_e = nullptr, *ljet_m = nullptr;
    ROOT::RVecF *ljet_bJR10v00_pt = nullptr, *ljet_bJR10v00_mass = nullptr;
    ROOT::RVecF *ljet_bJR10v01_pt = nullptr, *ljet_bJR10v01_mass = nullptr;
    ROOT::RVecF *jet_pt = nullptr, *jet_eta = nullptr, *jet_phi = nullptr;
    ROOT::RVecF *jet_e = nullptr, *jet_m = nullptr, *jet_jvt = nullptr;
    
    // pT and mass of the active jet type. The lambdas read the pointer at
    // call time, so they see the vector ROOT filled for the current entry.
    std::function<ROOT::RVecF*()> getCurrentPt;
    std::function<ROOT::RVecF*()> getCurrentMass;

    // Working point flags - (bb)-tagging
    ROOT::RVecF *dhbb = nullptr;          // GN2x (bb)-tagging score
    ROOT::RVecB *pass_FlatMassQCDEff_0p25 = nullptr;
    ROOT::RVecB *pass_FlatMassQCDEff_0p3 = nullptr;
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
    ROOT::RVecF *jet_eff_jvt             = nullptr;  
    ROOT::RVecF *matched_tjet_pt         = nullptr;  
    ROOT::RVecF *matched_tjet_m          = nullptr;  
    ROOT::RVecF *matched_tjet_eta        = nullptr;  
    ROOT::RVecF *matched_tjet_phi        = nullptr;  
    ROOT::RVecF *matched_tjet_dR         = nullptr;  

    // Background composition branches (MC only)
    ROOT::RVecI *ljet_n_bhadrons = nullptr;  // Number of b-hadrons ghost-associated to large-R jet
    ROOT::RVecI *ljet_n_chadrons = nullptr;  // Number of c-hadrons ghost-associated to large-R jet

    // Event topology
    int n_small_r_jets = 0;
    int n_large_r_jets = 0;
    
    // Angular variables
    double delta_phi_leading = -999, delta_phi_subleading = -999, delta_phi_subsubleading = -999;
    double delta_eta_overlap = 0, delta_phi_overlap = 0, delta_R_overlap = 0;
    double delta_eta_jet_jet = 0, delta_phi_jet_jet = 0, delta_R_jet_jet = 0;
    double delta_eta_jet_ljet = 0, delta_phi_jet_ljet = 0, delta_R_jet_ljet = 0;
    

    // Recoil system
    double pt_recoil_x = 0, pt_recoil_y = 0, pt_recoil = 0, eta_recoil = 0, phi_recoil = 0;
    
    // Jet ordering
    int leading_idx = -1, second_leading_idx = -1, third_leading_idx = -1;
    double leading_pt = -1, second_leading_pt = -1, third_leading_pt = -1;

    int index_tagged = -1;       
    int notag_index = -1;       
    int selected_ljet_index = -1;
    
    // Balance observables
    double pt_ref_MJB = 0;  
    double R_MJB = 0;       
    double Asym = 0;       

    double R_truth_MJB = -1.0;
    double pt_truth_jet = -1.0;
    double dR_truth_match = -1.0;

    double alpha = 0;
    double min_delta_R = 0;

    //==========================================================================
    // Set Up Input Tree Branches
    //==========================================================================
    
    // Common large-R jet branches
    t->SetBranchAddress("ljet_pt", &ljet_pt);
    t->SetBranchAddress("ljet_eta", &ljet_eta);
    t->SetBranchAddress("ljet_phi", &ljet_phi);
    t->SetBranchAddress("ljet_e", &ljet_e);
    t->SetBranchAddress("ljet_m", &ljet_m);
    
    // Set up branches based on jet type
    switch (jetType) {
        case JetType::BJER10V00:
            if (t->SetBranchAddress("ljet_bJR10v00_pt", &ljet_bJR10v00_pt) < 0 ||
                t->SetBranchAddress("ljet_bJR10v00_mass", &ljet_bJR10v00_mass) < 0) {
                std::cerr << "Error: Cannot find bJER10v00 branches in tree" << std::endl;
                return;
            }
            getCurrentPt   = [&]() { return ljet_bJR10v00_pt; };
            getCurrentMass = [&]() { return ljet_bJR10v00_mass; };
            std::cout << "Using bJER10v00 jet variables" << std::endl;
            break;
        case JetType::BJER10V01:
            if (t->SetBranchAddress("ljet_bJR10v01_pt", &ljet_bJR10v01_pt) < 0 ||
                t->SetBranchAddress("ljet_bJR10v01_mass", &ljet_bJR10v01_mass) < 0) {
                std::cerr << "Error: Cannot find bJER10v01 branches in tree" << std::endl;
                return;
            }
            getCurrentPt   = [&]() { return ljet_bJR10v01_pt; };
            getCurrentMass = [&]() { return ljet_bJR10v01_mass; };
            std::cout << "Using bJER10v01 jet variables" << std::endl;
            break;
        default: 
            getCurrentPt   = [&]() { return ljet_pt; };
            getCurrentMass = [&]() { return ljet_m; };
            std::cout << "Using nominal jet variables" << std::endl;
            break;
    }
    
    // Small-R jet branches
    t->SetBranchAddress("jet_pt", &jet_pt);
    t->SetBranchAddress("jet_eta", &jet_eta);
    t->SetBranchAddress("jet_phi", &jet_phi);
    t->SetBranchAddress("jet_e", &jet_e);
    t->SetBranchAddress("jet_m", &jet_m);
    t->SetBranchAddress("jet_jvt", &jet_jvt);
    
    // (bb)-tagging branches
    t->SetBranchAddress("dhbb", &dhbb);
    if (t->SetBranchAddress("pass_FlatMassQCDEff_0p25", &pass_FlatMassQCDEff_0p25) < 0 ||
        t->SetBranchAddress("pass_FlatMassQCDEff_0p3", &pass_FlatMassQCDEff_0p3) < 0 ||
        t->SetBranchAddress("pass_FlatMassQCDEff_0p37", &pass_FlatMassQCDEff_0p37) < 0 ||
        t->SetBranchAddress("pass_FlatMassQCDEff_0p46", &pass_FlatMassQCDEff_0p46) < 0 ||
        t->SetBranchAddress("pass_FlatMassQCDEff_0p58", &pass_FlatMassQCDEff_0p58) < 0 ||
        t->SetBranchAddress("pass_FlatMassQCDEff_0p74", &pass_FlatMassQCDEff_0p74) < 0 ||
        t->SetBranchAddress("pass_FlatMassQCDEff_0p94", &pass_FlatMassQCDEff_0p94) < 0 ||
        t->SetBranchAddress("pass_FlatMassQCDEff_1p25", &pass_FlatMassQCDEff_1p25) < 0 ||
        t->SetBranchAddress("pass_FlatMassQCDEff_1p55", &pass_FlatMassQCDEff_1p55) < 0) {
        std::cerr << "Error: Cannot find working point branches in tree" << std::endl;
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
            std::cerr << "Error: Cannot find truth matching branches in tree" << std::endl;
            return;
        }
        if (t->SetBranchAddress("ljet_n_bhadrons", &ljet_n_bhadrons) < 0 ||
            t->SetBranchAddress("ljet_n_chadrons", &ljet_n_chadrons) < 0) {
            std::cerr << "Error: Cannot find ljet_n_bhadrons/ljet_n_chadrons branches in tree" << std::endl;
            return;
        }
        std::cout << "Truth matching branches loaded (MC mode)" << std::endl;
    }

    //==========================================================================
    // Set Up Output Trees
    //==========================================================================
    
    // Create output file and trees
    TFile *of = new TFile(fn_out, "RECREATE");
    if (!of || of->IsZombie()) {
        std::cerr << "Error: Cannot create output file " << fn_out << std::endl;
        return;
    }
    of->cd();  
    auto ot = t->CloneTree(0);  

    // Create intermediate trees for cutflow validation
    std::vector<std::unique_ptr<TTree>> analysis_trees;
    const std::vector<std::string> tree_names = {
        "PreNominal", "PreCutNominal", "PreDxbb", "PreLargeRpT",
        "PreLargeReta", "PreSmallRpT", "PreJVT", "PreDeltaRvetoSmall",
        "PreVetoPtRSmall", "PreDeltaPhiJj", "PreDeltaRJj"
    };

    for (const auto& name : tree_names) {
        auto tree = std::unique_ptr<TTree>(t->CloneTree(0));
        tree->SetName(name.c_str());
        tree->SetDirectory(of);  
        analysis_trees.push_back(std::move(tree));
    }

    // Helper function to set up common branches for all trees
    auto setupBranches = [&](TTree* tree) {
        tree->Branch("index_tagged", &index_tagged, "index_tagged/I");
        tree->Branch("notag_index", &notag_index, "notag_index/I");
        tree->Branch("selected_ljet_index", &selected_ljet_index, "selected_ljet_index/I");
        tree->Branch("weight", &weight, "weight/D");
        tree->Branch("R_MJB", &R_MJB, "R_MJB/D");
        tree->Branch("delta_eta_overlap", &delta_eta_overlap, "delta_eta_overlap/D");
        tree->Branch("delta_phi_overlap", &delta_phi_overlap, "delta_phi_overlap/D");
        tree->Branch("delta_R_overlap", &delta_R_overlap, "delta_R_overlap/D");
        tree->Branch("pt_recoil", &pt_recoil, "pt_recoil/D");
        tree->Branch("eta_recoil", &eta_recoil, "eta_recoil/D");
        tree->Branch("phi_recoil", &phi_recoil, "phi_recoil/D");
        tree->Branch("alpha", &alpha, "alpha/D");
        tree->Branch("min_delta_R", &min_delta_R, "min_delta_R/D");
        tree->Branch("delta_phi_jet_ljet", &delta_phi_jet_ljet, "delta_phi_jet_ljet/D");
        tree->Branch("delta_eta_jet_ljet", &delta_eta_jet_ljet, "delta_eta_jet_ljet/D");
        tree->Branch("delta_R_jet_ljet", &delta_R_jet_ljet, "delta_R_jet_ljet/D");
        tree->Branch("delta_phi_leading", &delta_phi_leading, "delta_phi_leading/D");
        tree->Branch("delta_phi_subleading", &delta_phi_subleading, "delta_phi_subleading/D");
        tree->Branch("delta_phi_subsubleading", &delta_phi_subsubleading, "delta_phi_subsubleading/D");
        tree->Branch("delta_eta_jet_jet", &delta_eta_jet_jet, "delta_eta_jet_jet/D");
        tree->Branch("delta_phi_jet_jet", &delta_phi_jet_jet, "delta_phi_jet_jet/D");
        tree->Branch("delta_R_jet_jet", &delta_R_jet_jet, "delta_R_jet_jet/D");
        tree->Branch("Asym", &Asym, "Asym/D");
        tree->Branch("R_truth_MJB",     &R_truth_MJB,     "R_truth_MJB/D");
        tree->Branch("pt_truth_jet",    &pt_truth_jet,    "pt_truth_jet/D");
        tree->Branch("dR_truth_match",  &dR_truth_match,  "dR_truth_match/D");
    };

    // Set up branches for all trees
    for (auto& tree : analysis_trees) {
        setupBranches(tree.get());
    }
    setupBranches(ot);

    bool found_jets_with_eta_gt_2p8 = false; 

    std::cout << "Will loop over " << std::fixed << std::setprecision(1) 
              << 1e-6 * Nall << "M events\n";

    //==========================================================================
    // Event Loop
    //==========================================================================
    
    for (int i = 0; i < Nall; ++i) {
        if (i && i % 500000 == 0) {
            std::cout << std::fixed << std::setprecision(1)
                      << std::setw(6) << 1e-6 * i << "M / "
                      << std::setw(6) << 1e-6 * Nall << "M evts.  "
                      << std::setprecision(1) << 100.0 * i / Nall << "%\n";
            printTime(false);
        }

        t->GetEntry(i);
        
        // Initialize event weight
        weight = isData ? 1 : total_weight;
        if (weight > 100) continue; 
        cutflow["Cut_0"] += weight;
        eventCounter["Cut_0"]++;

        selected_ljet_index = -1;

        //======================================================================
        // Large-R Jet Selection
        //======================================================================
        
        // Reset event variables
        pt_recoil_x = pt_recoil_y = eta_recoil = phi_recoil = 0;
        n_small_r_jets = 0;
        leading_idx = second_leading_idx = third_leading_idx = -1;
        leading_pt = second_leading_pt = third_leading_pt = -1;
        delta_phi_leading = delta_phi_subleading = delta_phi_subsubleading = -999;

        R_truth_MJB    = -1.0;
        pt_truth_jet   = -1.0;
        dR_truth_match = -1.0;

        if (!getCurrentPt()) {
            std::cerr << "Error: getCurrentPt() returned null!" << std::endl;
            continue;
        }

        n_large_r_jets = getCurrentPt()->size();
        h_n_large_r_jets->Fill(n_large_r_jets);

        // Allow at most 2 large-R jets
        if (getCurrentPt()->size() > 2) continue;  

        // Find (bb)-tagged large-R jet
        bool is_tagged = false;
        size_t index_tagged_local = 0;  
        int notag_index_local = -1;    

        // Function to get the appropriate working point flag.
        auto getWorkingPointFlag = [&](size_t jetIndex) -> bool {
            ROOT::RVecB* wpVec =
                (workingPoint == "0p25") ? pass_FlatMassQCDEff_0p25 :
                (workingPoint == "0p30") ? pass_FlatMassQCDEff_0p3  :
                (workingPoint == "0p37") ? pass_FlatMassQCDEff_0p37 :
                (workingPoint == "0p46") ? pass_FlatMassQCDEff_0p46 :
                (workingPoint == "0p58") ? pass_FlatMassQCDEff_0p58 :
                (workingPoint == "0p74") ? pass_FlatMassQCDEff_0p74 :
                (workingPoint == "0p94") ? pass_FlatMassQCDEff_0p94 :
                (workingPoint == "1p25") ? pass_FlatMassQCDEff_1p25 :
                (workingPoint == "1p55") ? pass_FlatMassQCDEff_1p55 :
                                           pass_FlatMassQCDEff_0p46; 
            return (wpVec && jetIndex < wpVec->size()) ? wpVec->at(jetIndex) : false;
        };

        double highest_tagged_pt = -1;
        double highest_untagged_pt = -1;
        for (size_t j = 0; j < getCurrentPt()->size(); ++j) {
            if (getWorkingPointFlag(j)) {  // This jet is (bb)-tagged
                if (getCurrentPt()->at(j) > highest_tagged_pt) {
                    highest_tagged_pt = getCurrentPt()->at(j);
                    index_tagged_local = j;
                    is_tagged = true;
                }
            } else {  // This jet is not (bb)-tagged
                if (getCurrentPt()->at(j) > highest_untagged_pt) {
                    highest_untagged_pt = getCurrentPt()->at(j);
                    notag_index_local = j;
                }
            }
        }
        if (!is_tagged) continue;
        
        index_tagged = static_cast<int>(index_tagged_local);
        notag_index = static_cast<int>(notag_index_local);
        selected_ljet_index = index_tagged;

        // Background composition for the selected large-R jet (MC only)
        if (!isData) {
            int nb = ljet_n_bhadrons->at(index_tagged);
            int nc = ljet_n_chadrons->at(index_tagged);
            int jet_bin = (nb >= 2)                ? 0 :  // bb
                          (nb == 1 && nc >= 1)     ? 1 :  // bc
                          (nb == 1 && nc == 0)     ? 2 :  // b
                          (nb == 0 && nc >= 2)     ? 3 :  // cc
                          (nb == 0 && nc == 1)     ? 4 :  // c
                                                     5;   // light
            h_jet_counts->Fill(jet_bin, weight);
            h_background_composition->Fill(jet_bin, weight);
        }
        
        cutflow["Cut_1"] += weight;
        eventCounter["Cut_1"]++;
        analysis_trees[2]->Fill();  // PreDxbb

        // Apply kinematic cuts on tagged large-R jet
        if (getCurrentPt()->at(index_tagged) <= 450) continue;
        cutflow["Cut_2"] += weight;
        eventCounter["Cut_2"]++;
        analysis_trees[3]->Fill();  // PreLargeRpT

        if (std::abs(ljet_eta->at(index_tagged)) >= 1.2) continue;
        cutflow["Cut_3"] += weight;
        eventCounter["Cut_3"]++;
        analysis_trees[4]->Fill();  // PreLargeReta

        //======================================================================
        // Small-R Jet Selection
        //======================================================================
        
        std::vector<int> selected_jet_indices;
        for (size_t j = 0; j < jet_pt->size(); ++j) {
            
            // Calculate delta_R between this small-R jet and the (bb)-tagged large-R jet
            delta_eta_overlap = std::abs(jet_eta->at(j) - ljet_eta->at(index_tagged));
            delta_phi_overlap = DeltaPhi(jet_phi->at(j), ljet_phi->at(index_tagged));
            delta_R_overlap = std::sqrt(delta_eta_overlap * delta_eta_overlap + 
                                 delta_phi_overlap * delta_phi_overlap);
            
            // Basic kinematic selection
            if (jet_pt->at(j) > 25 && std::abs(jet_eta->at(j)) < 2.8 && delta_R_overlap > 1.) {
                // Calculate recoil system components
                pt_recoil_x += jet_pt->at(j) * std::cos(jet_phi->at(j));
                pt_recoil_y += jet_pt->at(j) * std::sin(jet_phi->at(j));
                ++n_small_r_jets;
                selected_jet_indices.push_back(j);

                // Determine pT ordering of jets
                if (jet_pt->at(j) > leading_pt) {
                    third_leading_pt = second_leading_pt;
                    third_leading_idx = second_leading_idx;
                    second_leading_pt = leading_pt;
                    second_leading_idx = leading_idx;
                    leading_pt = jet_pt->at(j);
                    leading_idx = j;
                } else if (jet_pt->at(j) > second_leading_pt) {
                    third_leading_pt = second_leading_pt;
                    third_leading_idx = second_leading_idx;
                    second_leading_pt = jet_pt->at(j);
                    second_leading_idx = j;
                } else if (jet_pt->at(j) > third_leading_pt) {
                    third_leading_pt = jet_pt->at(j);
                    third_leading_idx = j;
                }
            } else {
                // jet failed - check if it was specifically the eta cut
                if (jet_pt->at(j) > 25 && std::abs(jet_eta->at(j)) >= 2.8) {
                    found_jets_with_eta_gt_2p8 = true;
                }
            }
        }
        
        // Check if we have at least one selected jet before proceeding
        if (leading_idx < 0) continue;  // No selected jets remain after cuts

        cutflow["Cut_4"] += weight;
        eventCounter["Cut_4"]++;
        analysis_trees[5]->Fill();  // PreSmallRpT

        h_n_recoil_jets->Fill(n_small_r_jets);

        // Calculate total recoil pT
        pt_recoil = std::sqrt(pt_recoil_x * pt_recoil_x + pt_recoil_y * pt_recoil_y);

        //====================================================================================================
        // Calculate Angular Separations
        //====================================================================================================
        
        // Calculate Delta_Phi between each Small-R jet and the (bb)-tagged large-R jet
        if (leading_idx != -1) {
            delta_phi_leading = DeltaPhi(
                jet_phi->at(leading_idx), ljet_phi->at(index_tagged));
        }

        if (second_leading_idx != -1) {
            delta_phi_subleading = DeltaPhi(
                jet_phi->at(second_leading_idx), ljet_phi->at(index_tagged));
        }

        if (third_leading_idx != -1) {
            delta_phi_subsubleading = DeltaPhi(
                jet_phi->at(third_leading_idx), ljet_phi->at(index_tagged));
        }

        // Calculate separation between leading small-R jet and large-R jet
        delta_eta_jet_ljet = std::abs(jet_eta->at(leading_idx) - ljet_eta->at(index_tagged));
        delta_phi_jet_ljet = DeltaPhi(
            jet_phi->at(leading_idx), ljet_phi->at(index_tagged));
        
        delta_R_jet_ljet = std::sqrt(delta_eta_jet_ljet * delta_eta_jet_ljet + 
                               delta_phi_jet_ljet * delta_phi_jet_ljet);

        //====================================================================================================
        // JVT Selection
        //====================================================================================================
        
        // JVT applies only to low-pT jets within tracker acceptance (pT < 60 GeV, |eta| < 2.4)
        /*if (jet_pt->at(leading_idx) < 60 && 
            (*jet_jvt)[leading_idx] <= 0.59 && 
            std::abs(jet_eta->at(leading_idx)) < 2.4) {
            continue;  // Fail JVT
        }*/
        bool passJVT = true;
        for (int idx : selected_jet_indices) {
            if (jet_pt->at(idx) < 60 &&
                jet_jvt->at(idx) <= 0.59 &&
                std::abs(jet_eta->at(idx)) < 2.4) {
                passJVT = false;
                break;
            }
        }
        if (!passJVT) continue;
        cutflow["Cut_5"] += weight;
        eventCounter["Cut_5"]++;
        analysis_trees[6]->Fill();  // PreJVT

        //====================================================================================================
        // Multijets Balance Selection: Multijets Balance
        //====================================================================================================

        if (n_small_r_jets >= 1) {

            // Calculate recoil system properties based on number of jets
            if (n_small_r_jets == 1) {
                // For single jet: recoil pT is just the jet pT
                pt_recoil = jet_pt->at(leading_idx);
                eta_recoil = jet_eta->at(leading_idx);
                phi_recoil = jet_phi->at(leading_idx);
            } else {
                // For multiple jets: recoil is vector sum of all selected jets
                pt_recoil = std::sqrt(pt_recoil_x * pt_recoil_x + pt_recoil_y * pt_recoil_y);
                
                // Calculate eta and phi of recoil system
                float pz_recoil = 0.0;
                for (size_t j = 0; j < selected_jet_indices.size(); ++j) {
                    int idx = selected_jet_indices[j];
                    float jet_pz = jet_pt->at(idx) * std::sinh(jet_eta->at(idx));
                    pz_recoil += jet_pz;
                }
                
                // Calculate recoil system eta and phi
                float p_recoil = std::sqrt(pt_recoil_x * pt_recoil_x + pt_recoil_y * pt_recoil_y + pz_recoil * pz_recoil);
                eta_recoil = 0.5 * std::log((p_recoil + pz_recoil) / (p_recoil - pz_recoil));
                phi_recoil = std::atan2(pt_recoil_y, pt_recoil_x);
            }
            
            // alpha is the geometric angle between the b-jet and recoil system.
            // DeltaPhi returns a signed value in (-pi, pi], so we take std::abs to map it to [0, pi]
            // keeps events where alpha is within 0.3 rad (~17°) of π.
            alpha = std::abs(DeltaPhi(ljet_phi->at(index_tagged), phi_recoil));

            // Require alpha close to pi (back-to-back configuration)
            if (std::abs(alpha - PI) >= 0.3) continue;
            cutflow["Cut_6"] += weight;
            eventCounter["Cut_6"]++;
            analysis_trees[1]->Fill();  // PreCutNominal

            // Find minimum separation between any small-R jet and the tagged large-R jet
            min_delta_R = std::numeric_limits<double>::max();
            for (size_t j = 0; j < selected_jet_indices.size(); ++j) {
                int idx = selected_jet_indices[j];
                double delta_eta = jet_eta->at(idx) - ljet_eta->at(index_tagged);
                double delta_phi = DeltaPhi(jet_phi->at(idx), ljet_phi->at(index_tagged));
                double delta_R_jet_ljet_temp = std::sqrt(delta_eta * delta_eta + delta_phi * delta_phi);
                min_delta_R = std::min(min_delta_R, delta_R_jet_ljet_temp);
            }

            // Require minimum separation between all small-R jets and tagged large-R jet
            if (min_delta_R <= 1.5) continue;
            cutflow["Cut_7"] += weight;
            eventCounter["Cut_7"]++;
            analysis_trees[10]->Fill();  // PreDeltaRJj (after all cuts)

            // Calculate pt_ref_MJB based on number of small-R jets
            if (n_small_r_jets == 1) {
                // For single jet: use projected momentum
                pt_ref_MJB = jet_pt->at(leading_idx);
            } else {
                // For multiple jets: use total recoil system pT
                pt_ref_MJB = pt_recoil;
            }
        }

        //======================================================================
        // Calculate The Balance 
        //======================================================================
        
        // Calculate response (ratio of tagged jet pT to reference pT)
        R_MJB = (pt_ref_MJB != 0) ? getCurrentPt()->at(index_tagged) / pt_ref_MJB : 0;

        if (!isData && index_tagged >= 0 &&
            index_tagged < (int)matched_tjet_pt->size() &&
            matched_tjet_pt->at(index_tagged) > 0) {
            pt_truth_jet   = matched_tjet_pt->at(index_tagged);
            dR_truth_match = matched_tjet_dR->at(index_tagged);
            R_truth_MJB    = (pt_ref_MJB > 0) ? pt_truth_jet / pt_ref_MJB : -1.0;
        } else {
            pt_truth_jet = -1.0;
            dR_truth_match = -1.0;
            R_truth_MJB = -1.0;
        }

        // Fill final state histograms
        h_n_ljets->Fill(n_large_r_jets);
        h_n_jets->Fill(n_small_r_jets);

        ot->Fill();
    }

    //==========================================================================
    // Write Output
    //==========================================================================
    
    // Print final statistics and save results
    if (found_jets_with_eta_gt_2p8) {
        std::cout << "Warning: some jets with |eta| > 2.8 were removed by the eta cut." << std::endl;
    } else {
        std::cout << "No jets with |eta| > 2.8 found - eta cut had no effect." << std::endl;
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
            for (int i = 0; i <= 7; ++i) {
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
                outfile   << line;
                std::cout << line;
            }
            outfile.close();
        }
    }

    // Fill cutflow histogram
    for (int i = 0; i <= 7; ++i) {
        h_cutflow->SetBinContent(i + 1, cutflow["Cut_" + std::to_string(i)]);
    }

    //==========================================================================
    // Write Histograms to Output File
    //==========================================================================
    
    of->cd();
    
    // Write all histograms to output file
    h_cutflow->Write();
    h_n_large_r_jets->Write();
    h_n_recoil_jets->Write();
    h_n_ljets->Write();
    h_n_jets->Write();
    if (!isData) {
        h_jet_counts->Write();
        h_background_composition->Write();
    }

    //==========================================================================
    // Cleanup and Final Statistics
    //==========================================================================
    
    // Print final statistics
    std::cout << "\nProcessing completed for jet type: " << jetTypeSuffix
              << ", working point: " << workingPoint << "\n"
              << "Total events processed: " << Nall << "\n"
              << std::fixed << std::setprecision(1)
              << "Events passing all cuts: " << cutflow["Cut_7"] << "\n"
              << std::setprecision(2)
              << "Final efficiency: " << 100.0 * cutflow["Cut_7"] / cutflow["Cut_0"] << "%\n";

    printTime(false);
    
    ot->Show(0);
    of->Write();
    std::cout << "Output written to: " << fn_out << "\n";
    for (auto& tree : analysis_trees) tree.release();

    // Clean up
    gROOT->GetListOfFiles()->Remove(of);

    of->Close();
    delete of;  
    of = nullptr;  

    
    printTime();
}

//==============================================================================
// Convenience wrapper functions
//==============================================================================

void MJB_Nominal_Analysis(Str fn = "/data/atlas/tamezza/BJetCalib_Run3/JetCalibProcess/Zbbj/Zbbj_xbbCalib_19Jan2026/multijets_801165.root", bool isData = false, std::string workingPoint = "0p46") {
    ZbbJets_Insitu_Analysis(fn, isData, "nominal", workingPoint);
}

void MJB_bJER10v00_Analysis(Str fn = "/data/atlas/tamezza/BJetCalib_Run3/JetCalibProcess/Zbbj/Zbbj_xbbCalib_19Jan2026/multijets_801165.root", bool isData = false, std::string workingPoint = "0p46") {
    ZbbJets_Insitu_Analysis(fn, isData, "bJER10v00", workingPoint);
}

void MJB_bJER10v01_Analysis(Str fn = "/data/atlas/tamezza/BJetCalib_Run3/JetCalibProcess/Zbbj/Zbbj_xbbCalib_19Jan2026/multijets_801165.root", bool isData = false, std::string workingPoint = "0p46") {
    ZbbJets_Insitu_Analysis(fn, isData, "bJER10v01", workingPoint);
}

// Run all three jet types for a single working point
void MJB_All_JetTypes_SingleWP(Str fn = "/data/atlas/tamezza/BJetCalib_Run3/JetCalibProcess/Zbbj/Zbbj_xbbCalib_19Jan2026/multijets_801165.root", bool isData = false, std::string workingPoint = "0p46") {
    std::cout << "\n=== Running Analysis for All Jet Types (WP: " << workingPoint << ") ===\n" << std::endl;
    
    std::cout << "\n--- Analysis 1/3: Nominal Jets (WP: " << workingPoint << ") ---" << std::endl;
    ZbbJets_Insitu_Analysis(fn, isData, "nominal", workingPoint);
    
    std::cout << "\n--- Analysis 2/3: bJER10v00 Jets (WP: " << workingPoint << ") ---" << std::endl;
    ZbbJets_Insitu_Analysis(fn, isData, "bJER10v00", workingPoint);
    
    std::cout << "\n--- Analysis 3/3: bJER10v01 Jets (WP: " << workingPoint << ") ---" << std::endl;
    ZbbJets_Insitu_Analysis(fn, isData, "bJER10v01", workingPoint);
    
    std::cout << "\n=== All Jet Type Analyses Completed for WP: " << workingPoint << " ===\n" << std::endl;
}

// Run all jet types for all working points
void MJB_All_Analyses(Str fn = "/data/atlas/tamezza/BJetCalib_Run3/JetCalibProcess/Zbbj/Zbbj_xbbCalib_19Jan2026/multijets_801165.root", bool isData = false) {
    std::cout << "\n=== Running Analysis for All Jet Types and Working Points ===\n" << std::endl;
    
    std::vector<std::string> workingPoints = {"0p25", "0p30", "0p37", "0p46", "0p58", "0p74", "0p94", "1p25", "1p55"};
    std::vector<std::string> jetTypes = {"nominal", "bJER10v00", "bJER10v01"};
    
    int totalAnalyses = workingPoints.size() * jetTypes.size();
    int currentAnalysis = 0;
    
    for (const auto& wp : workingPoints) {
        for (const auto& jetType : jetTypes) {
            currentAnalysis++;
            std::cout << "\n--- Analysis " << currentAnalysis << "/" << totalAnalyses 
                     << ": " << jetType << " Jets (WP: " << wp << ") ---" << std::endl;
            ZbbJets_Insitu_Analysis(fn, isData, jetType, wp);
        }
    }
    
    std::cout << "\n=== All Analyses Completed ===\n" << std::endl;
}

// Run selected working points only
void MJB_Selected_WorkingPoints(Str fn = "/data/atlas/tamezza/BJetCalib_Run3/JetCalibProcess/Zbbj/Zbbj_xbbCalib_19Jan2026/multijets_801165.root", 
                               bool isData = false, 
                               std::vector<std::string> selectedWPs = {"0p37", "0p46", "0p58"}) {
    std::cout << "\n=== Running Analysis for Selected Working Points ===\n" << std::endl;
    
    std::vector<std::string> jetTypes = {"nominal", "bJER10v00", "bJER10v01"};
    
    int totalAnalyses = selectedWPs.size() * jetTypes.size();
    int currentAnalysis = 0;
    
    for (const auto& wp : selectedWPs) {
        for (const auto& jetType : jetTypes) {
            currentAnalysis++;
            std::cout << "\n--- Analysis " << currentAnalysis << "/" << totalAnalyses 
                     << ": " << jetType << " Jets (WP: " << wp << ") ---" << std::endl;
            ZbbJets_Insitu_Analysis(fn, isData, jetType, wp);
        }
    }
    
    std::cout << "\n=== Selected Working Point Analyses Completed ===\n" << std::endl;
}
