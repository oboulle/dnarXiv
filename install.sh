#!/bin/bash

git submodule init && git submodule update #clone the git subprojects

sudo apt-get update
sudo apt-get install xmlstarlet
