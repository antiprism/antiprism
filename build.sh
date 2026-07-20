#!/bin/sh
rm -rf $PWD/install
rm -rf $PWD/build
mkdir -p $PWD/build
cd $PWD/build
cmake -G "Unix Makefiles" -DCMAKE_INSTALL_PREFIX=../install ..
cmake --build . --target install -- -j4
cd ..

