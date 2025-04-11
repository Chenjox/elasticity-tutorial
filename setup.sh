#!/bin/bash

set -u

cd ..

PROJECT_DEPENDENCIES="'dune-common dune-grid dune-istl dune-geometry dune-localfunctions dune-functions dune-typetree dune-uggrid dune-vtk dune-gmsh4 dune-curvedgeometry dune-curvedgrid dune-assembler dune-fufem dune-matrix-vector dune-solvers'"

DUNE_CORE_VERSION="releases/2.10"
DUNE_STAGING_VERSION="v2.10.0"

DUNE_ASSEMBLER_VERSION="a7ea5486bcfecf1dcd3451473713f772a9c8c37e"

git clone --branch $DUNE_CORE_VERSION https://gitlab.dune-project.org/core/dune-common.git
git clone --branch $DUNE_CORE_VERSION https://gitlab.dune-project.org/core/dune-grid.git
git clone --branch $DUNE_CORE_VERSION https://gitlab.dune-project.org/core/dune-istl.git
git clone --branch $DUNE_CORE_VERSION https://gitlab.dune-project.org/core/dune-geometry.git
git clone --branch $DUNE_CORE_VERSION https://gitlab.dune-project.org/core/dune-localfunctions.git
git clone --branch $DUNE_STAGING_VERSION https://gitlab.dune-project.org/staging/dune-functions.git
git clone --branch $DUNE_CORE_VERSION https://gitlab.dune-project.org/staging/dune-typetree.git
git clone --branch $DUNE_CORE_VERSION https://gitlab.dune-project.org/staging/dune-uggrid.git
git clone --branch $DUNE_CORE_VERSION https://gitlab.dune-project.org/extensions/dune-vtk.git
git clone --branch $DUNE_CORE_VERSION https://gitlab.dune-project.org/extensions/dune-gmsh4.git
git clone --branch $DUNE_CORE_VERSION https://gitlab.dune-project.org/extensions/dune-curvedgeometry.git
git clone --branch $DUNE_CORE_VERSION https://gitlab.dune-project.org/extensions/dune-curvedgrid.git
git clone --branch $DUNE_CORE_VERSION https://gitlab.dune-project.org/fufem/dune-fufem.git
git clone --branch $DUNE_CORE_VERSION https://gitlab.dune-project.org/fufem/dune-matrix-vector.git
git clone --branch $DUNE_CORE_VERSION https://git.imp.fu-berlin.de/agnumpde/dune-solvers.git

# Dune Assembler ist so neu, das hat nicht einmal tags
git clone https://gitlab.dune-project.org/staging/dune-assembler.git
cd dune-assembler
git reset --hard $DUNE_ASSEMBLER_VERSION
cd ..

ln -s dune-common/bin/dunecontrol dunecontrol
ln -s dune-common/bin/duneproject duneproject
