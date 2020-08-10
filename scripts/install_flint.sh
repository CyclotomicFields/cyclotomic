#!/bin/bash -x

flint_version=2.6.2

sudo apt-get update
sudo apt-get install -y libgmp-dev libmpfr-dev
mkdir -p build
cd build
wget http://flintlib.org/flint-${flint_version}.tar.gz
tar -xvf flint-${flint_version}.tar.gz
cd flint-${flint_version}
./configure
make -j12
sudo make install
