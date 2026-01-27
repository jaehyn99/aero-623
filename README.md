# aero-623

Hi! This is the AEROSP 623 repo for our group. I've included a bit of my source code from my research project, so hopefully we won't have to reinvent the wheel. The codes
are built on the Eigen library, which consists matrix/vector classes and linear algebra routines. You will need to clone their source code down from their repo
(git@github.com:PX4/eigen.git) and put the entire `Eigen` subdirectory under `external` before you can compile.

Build: make -j

Projection test with arbitrary point: ./bin/main.exe proj1  projects/Project-1/bladeupper.txt projects/Project-1/bladelower.txt

Refinement: ./bin/main.exe refine projects/Project-1/mesh_coarse.gri projects/Project-1/bladeupper.txt projects/Project-1/bladelower.txt mesh_refined.gri 2 3

curve 1 is upperblade
curve 5 is lowerblade