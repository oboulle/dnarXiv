#!/bin/bash

#clone the git subprojects

git submodule init
git submodule update

# create the executables for the primer generation script

cd synthesis_modules/primer_generator/src
make
cd -

cd sequencing_modules/
git clone --recursive https://github.com/GATB/dsk.git
cd dsk
sh INSTALL
mkdir -p tmp
cd --

echo "installation completed"
