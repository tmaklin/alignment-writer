#!/bin/bash
## Build script for cross-compiling alignment-writer for macOS x86-64 or arm64.

set -exo pipefail

## Arg 1 should be the git tag to build
VER=$1
if [[ -z $VER ]]; then
  echo "Error: specify version"
  exit;
fi

## Arg 2 should be the target architecture
ARCH=$2
if [[ -z $ARCH ]]; then
  echo "Error: specify architecture (one of x86-64,arm64)"
  exit;
fi

## Install dependencies
apt update
apt install -y cmake git libomp5 libomp-dev curl

## Clone source tree
mkdir /io/tmp
cd /io/tmp
git clone https://github.com/tmaklin/alignment-writer
cd alignment-writer
git checkout ${VER}

## Configure source build
mkdir build
cd build
target_arch=""
if [ "$ARCH" = "x86-64" ]; then
    # compile x86_64
    cmake -DCMAKE_TOOLCHAIN_FILE="/io/$ARCH-toolchain_GNU.cmake" \
          -DCMAKE_C_FLAGS="-march=x86-64 -mtune=generic -m64 -fPIC -fPIE" \
          -DCMAKE_CXX_FLAGS="-march=x86-64 -mtune=generic -m64 -fPIC -fPIE" \
	  -DCMAKE_WITH_NATIVE_INSTRUCTIONS=OFF \
	  -DCMAKE_WITH_FLTO=1 \
          -DZLIB_LIBRARY="/osxcross/SDK/MacOSX13.0.sdk/usr/lib/libz.tbd" -DZLIB_INCLUDE_DIR="/osxcross/SDK/MacOSX13.0.sdk/usr/include" \
	  -DBZIP2_FOUND=0 \
	  -DLIBLZMA_FOUND=0 \
	  ..
    target_arch="x86_64-apple-darwin22"
elif [ "$ARCH" = "arm64" ]; then
    # compile aarch64
    cmake -DCMAKE_TOOLCHAIN_FILE="/io/$ARCH-toolchain_GNU.cmake" \
          -DCMAKE_C_FLAGS="-march=armv8-a -mtune=generic -m64 -fPIC -fPIE" \
          -DCMAKE_CXX_FLAGS="-march=armv8-a -mtune=generic -m64 -fPIC -fPIE" \
	  -DCMAKE_WITH_NATIVE_INSTRUCTIONS=OFF \
	  -DCMAKE_WITH_FLTO=1 \
          -DZLIB_LIBRARY="/osxcross/SDK/MacOSX13.0.sdk/usr/lib/libz.tbd" -DZLIB_INCLUDE_DIR="/osxcross/SDK/MacOSX13.0.sdk/usr/include" \
	  -DBZIP2_FOUND=0 \
	  -DLIBLZMA_FOUND=0 \
	  ..
    target_arch="aarch64-apple-darwin22"
fi

## Compile
make VERBOSE=1 -j

## Gather distributables
target=alignment-writer-${VER}-$target_arch
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

