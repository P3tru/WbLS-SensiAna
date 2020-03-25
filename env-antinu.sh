#!/bin/bash
echo "Setting up ROOT 6"
source /data/snoplus/home/zsoldos/.local/root-6.06.08/build/bin/thisroot.sh
#export DISPLAY=:0
echo "Setting up Geant4 10.01.p03"
export G4BUILD=/data/snoplus/home/zsoldos/.local/geant4.10.01.p03-install
source ${G4BUILD}/bin/geant4.sh; source ${G4BUILD}/share/Geant4-10.1.3/geant4make/geant4make.sh

export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/data/snoplus/home/zsoldos/.local/boost-1.71.0/lib:/data/snoplus/shared/code/hdf5/lib
source env.sh
