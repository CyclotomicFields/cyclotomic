#!/bin/bash -x

antic_version=0.2.2

sudo apt install -y libgmp-dev libmpfr-dev libflint-dev
mkdir build
cd build
wget https://github.com/wbhart/antic/archive/${antic_version}.tar.gz
tar -xvf ${antic_version}.tar.gz
cd antic-${antic_version}
./configure
make
sudo make install
