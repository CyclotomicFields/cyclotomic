#!/bin/bash -x

antic_version=0.2.2

apt install -y libgmp-dev libmpfr-dev libflint-2.6.2 libflint-dev
mkdir build
cd build
wget https://github.com/wbhart/antic/archive/${antic_version}.tar.gz
cd antic-${antic_version}
make
make install
