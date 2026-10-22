#!/bin/bash

# Run era - used to select the correct input ntuple subdir and output plot subdir - drives RunEra.h: run2|run3
export EASYBJETS_RUN=run2

run_root_analysis() {
    root -l -b -q "$1" 2>&1
    local exit_code=$?
    if [ $exit_code -ne 0 ] && [ $exit_code -ne 139 ]; then
        echo "ERROR: ROOT crashed with exit code $exit_code"
        return $exit_code
    fi
    return 0
}

# Files
MC_FILES=(
    "/data/atlas/tamezza/bJets_Insitu_Analysis/bJets_Samples/Run2/Zbbj/GN2x_JetCalibxNtuple_ZbbJets/ZbbJets.root"
    "/data/atlas/tamezza/bJets_Insitu_Analysis/bJets_Samples/Run2/Zbbj/GN2x_JetCalibxNtuple_ZbbJets/multijets_364702.root"
    "/data/atlas/tamezza/bJets_Insitu_Analysis/bJets_Samples/Run2/Zbbj/GN2x_JetCalibxNtuple_ZbbJets/multijets_364703.root"
    "/data/atlas/tamezza/bJets_Insitu_Analysis/bJets_Samples/Run2/Zbbj/GN2x_JetCalibxNtuple_ZbbJets/multijets_364704.root"
    "/data/atlas/tamezza/bJets_Insitu_Analysis/bJets_Samples/Run2/Zbbj/GN2x_JetCalibxNtuple_ZbbJets/multijets_364705.root"
    "/data/atlas/tamezza/bJets_Insitu_Analysis/bJets_Samples/Run2/Zbbj/GN2x_JetCalibxNtuple_ZbbJets/multijets_364706.root"
    "/data/atlas/tamezza/bJets_Insitu_Analysis/bJets_Samples/Run2/Zbbj/GN2x_JetCalibxNtuple_ZbbJets/multijets_364707.root"
    "/data/atlas/tamezza/bJets_Insitu_Analysis/bJets_Samples/Run2/Zbbj/GN2x_JetCalibxNtuple_ZbbJets/multijets_364708.root"
    "/data/atlas/tamezza/bJets_Insitu_Analysis/bJets_Samples/Run2/Zbbj/GN2x_JetCalibxNtuple_ZbbJets/multijets_364709.root"
    "/data/atlas/tamezza/bJets_Insitu_Analysis/bJets_Samples/Run2/Zbbj/GN2x_JetCalibxNtuple_ZbbJets/multijets_364710.root"
    "/data/atlas/tamezza/bJets_Insitu_Analysis/bJets_Samples/Run2/Zbbj/GN2x_JetCalibxNtuple_ZbbJets/multijets_364711.root"
    "/data/atlas/tamezza/bJets_Insitu_Analysis/bJets_Samples/Run2/Zbbj/GN2x_JetCalibxNtuple_ZbbJets/multijets_364712.root"
    "/data/atlas/tamezza/bJets_Insitu_Analysis/bJets_Samples/Run2/Zbbj/GN2x_JetCalibxNtuple_ZbbJets/multijets_sherpa.root"
    "/data/atlas/tamezza/bJets_Insitu_Analysis/bJets_Samples/Run2/Zbbj/GN2x_JetCalibxNtuple_ZbbJets/top.root"
    "/data/atlas/tamezza/bJets_Insitu_Analysis/bJets_Samples/Run2/Zbbj/GN2x_JetCalibxNtuple_ZbbJets/WqqJets.root"
    "/data/atlas/tamezza/bJets_Insitu_Analysis/bJets_Samples/Run2/Zbbj/GN2x_JetCalibxNtuple_ZbbJets/ZqqJets.root"
)

DATA_FILES=(
    "/data/atlas/tamezza/bJets_Insitu_Analysis/bJets_Samples/Run2/Zbbj/GN2x_JetCalibxNtuple_ZbbJets/data.root"
)

# Working points and jet types
#WORKING_POINTS=("0p25" "0p30" "0p37" "0p58" "0p74" "0p94" "0p74")
WORKING_POINTS=("0p46")
JET_TYPES=("nominal" "bJER10v00" "bJER10v01")
#JET_TYPES=("nominal")

# Process MC files
"=== Processing MC Files for Selected Working Points: ${WORKING_POINTS[@]} ==="
for file in "${MC_FILES[@]}"; do
    for wp in "${WORKING_POINTS[@]}"; do
        for jet_type in "${JET_TYPES[@]}"; do
            echo "Processing: $(basename $file) (WP: $wp, Jet: $jet_type)"
            # root -l -b -q "ZbbJets_Insitu_Analysis.cxx(\"$file\", false, \"$jet_type\", \"$wp\")"
            run_root_analysis "ZbbJets_Insitu_Analysis.cxx(\"$file\", false, \"$jet_type\", \"$wp\")"
        done
    done
done

# Process DATA files
echo "=== Processing Data Files for Selected Working Points: ${WORKING_POINTS[@]} ==="
for file in "${DATA_FILES[@]}"; do
    for wp in "${WORKING_POINTS[@]}"; do
        for jet_type in "${JET_TYPES[@]}"; do
            echo "Processing Data: $file (WP: $wp, Jet: $jet_type)"
            root -l -b -q "ZbbJets_Insitu_Analysis.cxx(\"$file\", true, \"$jet_type\", \"$wp\")"
        done
    done
done

echo "=== All analyses completed ==="
