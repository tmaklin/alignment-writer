#!/bin/bash

## Build script for compiling alignment-writer
## inside the Holy Build Box v3.0.2 container.

set -e

## Arg 1 should be the git tag to build
VER=$1
if [[ -z $VER ]]; then
  echo "Error: specify version"
  exit;
fi

## Install git (HBB git doesn't have https support)
yum -y install git
rm --force /hbb/bin/git
ln -s /usr/bin/git /hbb/bin/git

## Activate Holy Build Box environment.
source /hbb_exe/activate
export LDFLAGS="-L/lib64 -static-libstdc++"
set -x

## Clone source tree
mkdir /io/tmp
cd /io/tmp
git clone https://github.com/tmaklin/alignment-writer
cd alignment-writer
git checkout ${VER}

## Configure source build
mkdir build
cd build
cmake -DCMAKE_CXX_FLAGS="-march=x86-64 -mtune=generic -m64" \
      -DCMAKE_C_FLAGS="-march=x86-64 -mtune=generic -m64" \
      -DCMAKE_WITH_FLTO=1 \
      -DCMAKE_WITH_NATIVE_INSTRUCTIONS=OFF \
      ..

## Compile binary
make VERBOSE=1 -j

## Gather distributables
target=alignment-writer-${VER}-$(gcc -v 2>&1 | grep "^Target" | cut -f2 -d':' | sed 's/[[:space:]]*//g')
path=/io/tmp/$target
mkdir $path
cp ../build/bin/alignment-writer $path/
cp ../README.md $path/
cp ../LICENSE $path/
cd /io/tmp
tar -zcvf $target.tar.gz $target
mv $target.tar.gz /io/
cd /io/
rm -rf tmp cache

