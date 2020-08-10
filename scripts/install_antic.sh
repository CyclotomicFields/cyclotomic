#!/bin/bash -x

antic_version=0.2.2

mkdir -p build
cd build
wget https://github.com/wbhart/antic/archive/${antic_version}.tar.gz
tar -xvf ${antic_version}.tar.gz
cd antic-${antic_version}
./configure
make -j12
sudo make install
