#!/bin/bash

git submodule init && git submodule update #clone the git subprojects

sudo apt-get update
sudo apt-get install xmlstarlet

#install julia
if [ ! -f /usr/local/bin/julia ]; then
	wget https://julialang-s3.julialang.org/bin/linux/x64/1.7/julia-1.7.2-linux-x86_64.tar.gz
	tar -xvzf julia-1.7.2-linux-x86_64.tar.gz
	sudo cp -r julia-1.7.2 /opt/
	sudo ln -s /opt/julia-1.7.2/bin/julia /usr/local/bin/julia
	rm -rf julia-*
fi
