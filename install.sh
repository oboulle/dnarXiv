#!/bin/bash

#clone the git subprojects

git submodule init
git submodule update

cd synthesis_modules/primer_generator/src
make

echo "installation completed"
