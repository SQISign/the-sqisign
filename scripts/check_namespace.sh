#!/bin/bash

set -e

if [ ! -f "include/sqisign_namespace.h" ]; then
  echo "Please run script from the sqisign root directory"
  exit 1
fi

if [[ "$OSTYPE" == "darwin"* ]]; then
    sed -i '' 's|//#define DISABLE_NAMESPACING|#define DISABLE_NAMESPACING|' ./include/sqisign_namespace.h
else
    sed -i 's|//#define DISABLE_NAMESPACING|#define DISABLE_NAMESPACING|' ./include/sqisign_namespace.h
fi

mkdir -p build_broadwell && cd build_broadwell && cmake -DSQISIGN_BUILD_TYPE=broadwell .. && make -j8 && cd ..
mkdir -p build && cd build && cmake .. && make -j8
find . ../build_broadwell -name '*.a' -exec nm {} \; | grep '.c.o:\|T ' | scala -nc ../scripts/Namespace.scala > sqisign_namespace.h

if [[ "$OSTYPE" == "darwin"* ]]; then
    sed -i '' 's|#define DISABLE_NAMESPACING|//#define DISABLE_NAMESPACING|' ../include/sqisign_namespace.h
else
    sed -i 's|#define DISABLE_NAMESPACING|//#define DISABLE_NAMESPACING|' ../include/sqisign_namespace.h
fi

diff sqisign_namespace.h ../include/sqisign_namespace.h


# Check the exit code of diff
if [ $? -eq 0 ]; then
  echo "No change in namespace."
  exit 0
else
  echo "Namespace changed, please update."
  exit 1
fi
