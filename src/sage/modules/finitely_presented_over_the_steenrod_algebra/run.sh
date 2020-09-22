#!/bin/sh

# This is the location of my Sage source installation:
cd $HOME/sage/

# Run doctests on a single file:
./sage -t ./src/sage/modules/finitely_presented_over_the_steenrod_algebra/resolutions.py
#./sage -t ./src/sage/modules/finitely_presented_over_the_steenrod_algebra/module.py

# Run doctests on four different files in parallell:
./sage -tp 4 ./src/sage/modules/finitely_presented_over_the_steenrod_algebra/module.py ./src/sage/modules/finitely_presented_over_the_steenrod_algebra/element.py ./src/sage/modules/finitely_presented_over_the_steenrod_algebra/morphism.py ./src/sage/modules/finitely_presented_over_the_steenrod_algebra/homspace.py

# Build the reference manual for modules:
# ./sage --docbuild reference/modules html

