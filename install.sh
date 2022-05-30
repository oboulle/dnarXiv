#!/bin/bash

#clone the git subprojects

git submodule init
git submodule update

# create the executables for the primer generation script

cd synthesis_modules/primer_generator/src
make

echo "installation completed"
