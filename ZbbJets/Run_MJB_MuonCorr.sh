#!/bin/bash

# Run era - used to select the correct input ntuple subdir and output plot subdir - drives RunEra.h: run2|run3
export EASYBJETS_RUN=run3

# Working points
WORKING_POINTS=("0p46")
#WORKING_POINTS=("0p25" "0p30" "0p37" "0p46" "0p58" "0p74" "0p94")

run_root_analysis() {
    root -l -b -q "$1" 2>&1
    local exit_code=$?
    if [ $exit_code -eq 139 ]; then
        echo "WARNING: ROOT segfaulted (139) on: $1 — output may be truncated."
        return 0
    fi
    if [ $exit_code -ne 0 ]; then
        echo "ERROR: ROOT crashed with exit code $exit_code"
        return $exit_code
    fi
    return 0
}

# Files (from the _muonCorr directory)
MC_FILES=(
    "/data/atlas/tamezza/bJets_Insitu_Analysis/bJets_Samples/Run3/Zbbj/GN2x_JetCalibxNtuple_ZbbJets_run3_muonCorr/multijets_801165.root"
    "/data/atlas/tamezza/bJets_Insitu_Analysis/bJets_Samples/Run3/Zbbj/GN2x_JetCalibxNtuple_ZbbJets_run3_muonCorr/multijets_801166.root"
    "/data/atlas/tamezza/bJets_Insitu_Analysis/bJets_Samples/Run3/Zbbj/GN2x_JetCalibxNtuple_ZbbJets_run3_muonCorr/multijets_801167.root"
    "/data/atlas/tamezza/bJets_Insitu_Analysis/bJets_Samples/Run3/Zbbj/GN2x_JetCalibxNtuple_ZbbJets_run3_muonCorr/multijets_801168.root"
    "/data/atlas/tamezza/bJets_Insitu_Analysis/bJets_Samples/Run3/Zbbj/GN2x_JetCalibxNtuple_ZbbJets_run3_muonCorr/multijets_801169.root"
    "/data/atlas/tamezza/bJets_Insitu_Analysis/bJets_Samples/Run3/Zbbj/GN2x_JetCalibxNtuple_ZbbJets_run3_muonCorr/multijets_801170.root"
    "/data/atlas/tamezza/bJets_Insitu_Analysis/bJets_Samples/Run3/Zbbj/GN2x_JetCalibxNtuple_ZbbJets_run3_muonCorr/multijets_801171.root"
    "/data/atlas/tamezza/bJets_Insitu_Analysis/bJets_Samples/Run3/Zbbj/GN2x_JetCalibxNtuple_ZbbJets_run3_muonCorr/multijets_801172.root"
    "/data/atlas/tamezza/bJets_Insitu_Analysis/bJets_Samples/Run3/Zbbj/GN2x_JetCalibxNtuple_ZbbJets_run3_muonCorr/multijets_801173.root"
    "/data/atlas/tamezza/bJets_Insitu_Analysis/bJets_Samples/Run3/Zbbj/GN2x_JetCalibxNtuple_ZbbJets_run3_muonCorr/ZbbJets.root"
    #"/data/atlas/tamezza/bJets_Insitu_Analysis/bJets_Samples/Run2/Zbbj/GN2x_JetCalibxNtuple_ZbbJets_muonCorr/top.root"
    #"/data/atlas/tamezza/bJets_Insitu_Analysis/bJets_Samples/Run2/Zbbj/GN2x_JetCalibxNtuple_ZbbJets_muonCorr/WqqJets.root"
    #"/data/atlas/tamezza/bJets_Insitu_Analysis/bJets_Samples/Run2/Zbbj/GN2x_JetCalibxNtuple_ZbbJets_muonCorr/ZqqJets.root"
)

DATA_FILES=(
    "/data/atlas/tamezza/bJets_Insitu_Analysis/bJets_Samples/Run3/Zbbj/GN2x_JetCalibxNtuple_ZbbJets_run3_muonCorr/data.root"
)

# === MuonCorr MC ===
echo "=== Processing MJB MuonCorr MC for Working Points: ${WORKING_POINTS[@]} ==="
for file in "${MC_FILES[@]}"; do
    if [[ ! -f "$file" ]]; then
        echo "SKIP (missing): $(basename $file)"
        continue
    fi
    for wp in "${WORKING_POINTS[@]}"; do
        echo "Processing: $(basename $file) (WP: $wp, MuonCorr)"
        run_root_analysis "ZbbJets_Insitu_Analysis.cxx(\"$file\", false, \"nominal\", \"$wp\", true)"
    done
done

# === MuonCorr Data ===
echo "=== Processing MJB MuonCorr Data for Working Points: ${WORKING_POINTS[@]} ==="
for file in "${DATA_FILES[@]}"; do
    if [[ ! -f "$file" ]]; then
        echo "SKIP (missing): $(basename $file)"
        continue
    fi
    for wp in "${WORKING_POINTS[@]}"; do
        echo "Processing Data: $(basename $file) (WP: $wp, MuonCorr)"
        run_root_analysis "ZbbJets_Insitu_Analysis.cxx(\"$file\", true, \"nominal\", \"$wp\", true)"
    done
done

echo "=== MJB MuonCorr processing complete ==="
