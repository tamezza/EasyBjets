#!/bin/bash

# Run era - used to select the correct input ntuple subdir and output plot subdir - drives RunEra.h: run2|run3
export EASYBJETS_RUN=run3

# Files
MC_FILES=(
    "/data/atlas/tamezza/bJets_Insitu_Analysis/bJets_Samples/Run3/Zbby/GN2x_JetCalibxNtuple_Zbby_run3/Zbby.root"
    "/data/atlas/tamezza/bJets_Insitu_Analysis/bJets_Samples/Run3/Zbby/GN2x_JetCalibxNtuple_Zbby_run3/gamma_jets_pythia.root"
    "/data/atlas/tamezza/bJets_Insitu_Analysis/bJets_Samples/Run3/Zbby/GN2x_JetCalibxNtuple_Zbby_run3/tty.root"
    "/data/atlas/tamezza/bJets_Insitu_Analysis/bJets_Samples/Run3/Zbby/GN2x_JetCalibxNtuple_Zbby_run3/Wqqy.root"
    "/data/atlas/tamezza/bJets_Insitu_Analysis/bJets_Samples/Run3/Zbby/GN2x_JetCalibxNtuple_Zbby_run3/Zqqy.root"
    #"/data/atlas/tamezza/bJets_Insitu_Analysis/bJets_Samples/Run3/Zbby/GN2x_JetCalibxNtuple_Zbby_run3/Zqqy_MG.root"

)

DATA_FILES=(
    "/data/atlas/tamezza/bJets_Insitu_Analysis/bJets_Samples/Run3/Zbby/GN2x_JetCalibxNtuple_Zbby_run3/data.root"
)

# Working points and jet types
WORKING_POINTS=("0p46")
#WORKING_POINTS=("0p25" "0p30" "0p37" "0p46" "0p58" "0p74" "0p94" "1p25" "1p55")
JET_TYPES=("nominal" "bJER10v00" "bJER10v01")

# Process MC files
echo "=== Processing MC Files for Selected Working Points: ${WORKING_POINTS[@]} ==="
for file in "${MC_FILES[@]}"; do
    for wp in "${WORKING_POINTS[@]}"; do
        for jet_type in "${JET_TYPES[@]}"; do
            echo "Processing: $(basename $file) (WP: $wp, Jet: $jet_type)"
            root -l -b -q "Zbby_Insitu_Analysis.cxx(\"$file\", false, \"$jet_type\", \"$wp\")"
        done
    done
done

# Process DATA files
echo "=== Processing Data Files for Selected Working Points: ${WORKING_POINTS[@]} ==="
for file in "${DATA_FILES[@]}"; do
    for wp in "${WORKING_POINTS[@]}"; do
        for jet_type in "${JET_TYPES[@]}"; do
            echo "Processing Data: $file (WP: $wp, Jet: $jet_type)"
            root -l -b -q "Zbby_Insitu_Analysis.cxx(\"$file\", true, \"$jet_type\", \"$wp\")"
        done
    done
done

echo "=== All analyses completed ==="
