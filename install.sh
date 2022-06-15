#!/bin/bash


# create the executables for the primer generation script

cd synthesis_modules/primer_generator/src
make
cd -

# create the executables for the kmer counting script DSK

cd sequencing_modules/
git clone --recursive https://github.com/GATB/dsk.git
cd dsk
sh INSTALL
mkdir -p tmp
cd --

echo "installation completed"
