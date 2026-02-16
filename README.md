# aero-623

Hi! This is the AEROSP 623 repo for our group. I've included a bit of my source code from my research project, so hopefully we won't have to reinvent the wheel. The codes
are built on the Eigen library, which consists matrix/vector classes and linear algebra routines. You will need to clone their source code down from their repo
(git@github.com:PX4/eigen.git) and put the entire `Eigen` subdirectory under `external` before you can compile.

# Airfoil Geometry
curve 1 is upperblade
curve 5 is lowerblade

# In/Out
curve 3: In
curve 7: Out

# Periodic
curve 2: left bottom
curve 4: left top
curve 8: right bottom
curve 6: right top

# Build: in "/aero-623"
make clean
make -j

# Verification: in "/aero-623"
./bin/main.exe
# or
./bin/main.exe verify

# Refinement: in "/aero-623"
./bin/main.exe refine projects/Project-1/mesh_coarse.gri projects/Project-1/bladeupper.txt projects/Project-1/bladelower.txt projects/Project-1/mesh_refined_out.gri  (alpha-sclar for size function) (refinement iteration)

# Run the first order solver
./bin/main.exe --mode steady-local --prefix projects/Project-1/mesh_refined_2394
./bin/main.exe --mode steady-global --mesh projects/Project-1/mesh_coarse.gri