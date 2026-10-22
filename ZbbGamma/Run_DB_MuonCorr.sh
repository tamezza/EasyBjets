#!/bin/bash

# Run era - used to select the correct input ntuple subdir and output plot subdir - drives RunEra.h: run2|run3
export EASYBJETS_RUN=run3

# Files
MC_FILES=(
    "/data/atlas/tamezza/bJets_Insitu_Analysis/bJets_Samples/Run3/Zbby/GN2x_JetCalibxNtuple_Zbby_run3_muonCorr/Zbby.root"
    "/data/atlas/tamezza/bJets_Insitu_Analysis/bJets_Samples/Run3/Zbby/GN2x_JetCalibxNtuple_Zbby_run3_muonCorr/gamma_jets_pythia.root"
    "/data/atlas/tamezza/bJets_Insitu_Analysis/bJets_Samples/Run3/Zbby/GN2x_JetCalibxNtuple_Zbby_run3_muonCorr/tty.root"
    "/data/atlas/tamezza/bJets_Insitu_Analysis/bJets_Samples/Run3/Zbby/GN2x_JetCalibxNtuple_Zbby_run3_muonCorr/Wqqy.root"
    "/data/atlas/tamezza/bJets_Insitu_Analysis/bJets_Samples/Run3/Zbby/GN2x_JetCalibxNtuple_Zbby_run3_muonCorr/Zqqy.root"
    "/data/atlas/tamezza/bJets_Insitu_Analysis/bJets_Samples/Run3/Zbby/GN2x_JetCalibxNtuple_Zbby_run3_muonCorr/Zqqy_MG.root"

)

DATA_FILES=(
    "/data/atlas/tamezza/bJets_Insitu_Analysis/bJets_Samples/Run3/Zbby/GN2x_JetCalibxNtuple_Zbby_run3_muonCorr/data.root"
)

#WORKING_POINTS=("0p25" "0p30" "0p37" "0p58" "0p74" "0p94")
#WORKING_POINTS=("0p46")
WORKING_POINTS=("0p25" "0p30" "0p37" "0p46" "0p58" "0p74" "0p94" "1p25" "1p55")

# === MuonCorr MC ===
echo "=== Processing MuonCorr MC ==="
for file in "${MC_FILES[@]}"; do
    for wp in "${WORKING_POINTS[@]}"; do
        echo "Processing: $(basename $file) (WP: $wp, MuonCorr)"
        root -l -b -q "Zbby_Insitu_Analysis.cxx(\"$file\", false, \"nominal\", \"$wp\", true)"
    done
done

# === MuonCorr Data ===
echo "=== Processing MuonCorr Data ==="
for file in "${DATA_FILES[@]}"; do
    for wp in "${WORKING_POINTS[@]}"; do
        echo "Processing Data: $(basename $file) (WP: $wp, MuonCorr)"
        root -l -b -q "Zbby_Insitu_Analysis.cxx(\"$file\", true, \"nominal\", \"$wp\", true)"
    done
done

echo "=== MuonCorr processing complete ==="
