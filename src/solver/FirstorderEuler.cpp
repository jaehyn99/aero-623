"""
Data structure for first order euler method

Mesh: /home/jaehyn/CFD2/aero-623/projects/Project-1/mesh_refined_2394.gri
nNode nElemTot Dim
for i = 1:nNode
    x(i) y(i) [z(i)]
nBGroup
for i = 1:nBGroup
    nBFace(i) nf(i) [Title(i)]
    for j = 1:nBFace(i)
        NB(i,j,1) .. NB(i,j,nf(i))
for i = 1:nElemGroup
    nElem(i) Order(i) Basis(i)
    for j = 1:nElem(i)
        NE(i,j,1) .. NE(i,j,nn(i))
nPG "PeriodicGroup"
for i = 1:nPG
    nPGNode(i) Periodicity(i)
    for j = 1:nPGNode(i)
        NP(i,j,1) NP(i,j,2)

I2E: /home/jaehyn/CFD2/aero-623/projects/Project-1/mesh_refined_2394I2E.txt
for i = 1:nIntFaceTot
    elemL(i) faceL(i) elemR(i) faceR(i)

B2E: /home/jaehyn/CFD2/aero-623/projects/Project-1/mesh_refined_2394B2E.txt
for i = 1:nBFaceTot
    elemIndex(i) nBFace(i) nBGroup(i)

In: /home/jaehyn/CFD2/aero-623/projects/Project-1/mesh_refined_2394In.txt
for  i = 1:nIntFaceTot
    IN(x(i),y(i))

Bn: /home/jaehyn/CFD2/aero-623/projects/Project-1/mesh_refined_2394Bn.txt
for i = 1:nBFaceTot
    BN(x(i),y(i))

Area: /home/jaehyn/CFD2/aero-623/projects/Project-1/mesh_refined_2394Area.txt
for i = 1:nElemTot
    Area(i)

Periodic Edges:
/home/jaehyn/CFD2/aero-623/projects/Project-1/mesh_refined_2394periodicEdges.txt
"""

#include "FirstorderEuler.h"

FirstorderEuler::FirstorderEuler(const std::string& meshFile, const std::string& areaFile, const std::string& b2eFile, const std::string& bnFile, const std::string& i2eFile, const std::string& inFile, const std::string& periodicEdgesFile)
    : mesh(meshFile), area(areaFile), b2e(b2eFile), bn(bnFile), i2e(i2eFile), in(inFile), periodicEdges(periodicEdgesFile) {

    // Input: mesh_file='mesh_refined_2394.gri', Tfinal=2.0, CFL=0.8, g=9.8, save_every=0, out_prefix='sol', ...

    // Process mesh
    choose the mesh file and read in the mesh data (nodes, elements, boundary groups, etc.)
    // Process normals and length of edges
        // interior-edge normals and lengths (outward for left elem, inward for right)
        // boundary-edge normals and lengths (outward from the cell)
    // Process areas and centroids (triangle area (positive))
    
    // Initialize parameters (Input parameters can be implemented as command line arguments)
    gamma = 1.4 // Ratio of specific heats
    rho0 = 1 // inlet stagnation density
    a0 = 1 // inlet stagnation speed of sound
    p0 = 1 // inlet stagnation pressure
    alpha = 0.0 // flow angle in radians
    pout = 1 // outlet static pressure
    CFL = 0.8 // CFL number for time-stepping
    Tfinal = 2.0 // final time for unsteady simulations
    save_every = 0 // save solution every N iterations (0 for no saving)
    out_prefix = "sol" // prefix for output files
    M = 0.1 // initial guess for Mach number for initializing the flow field

    // Define conservative variables (Unknown is the cell average conservative state)
    U[] = [rho, rhou, rhov, rhoE] // Density, Momentum in x, Momentum in y, Energy
    
    // Initialize U everywhere from inlet stagnation conditions and guessed M =~ 0.1
    // Initialize the flow to a uniform state obtained using the inlet conditions and a guessed Mach number of 0.1. For uniform-refinement studies, you
    // can use the solution from a coarser mesh to initialize the state on a finer mesh.

    for i = 1:nElemTot
        u = M * a0 * cos(alpha) // x-velocity
        v = M * a0 * sin(alpha) // y-velocity
        T = T0 / (1 + 0.5 * (gamma - 1) * M^2) // static temperature from stagnation temperature and Mach number
        p = p0 * (T / T0)^(gamma / (gamma - 1)) // static pressure from stagnation pressure and temperature

        U(i, rho) = p0 / (R * T0)
        U(i, rhou) = U(i, rho) * u
        U(i, rhov) = U(i, rho) * v
        U(i, rhoE) = p/ (gamma - 1) + 0.5 * U(i, rho) * (u^2 + v^2) // total energy from static pressure and kinetic energy
    end for

    // Initialize Residual and Flux arrays
    Residual(i) = 0.0
    Flux(i, j, bn) = 0.0 // From flux function defined from collaborators
    Flux(i, j, in) = 0.0 // From flux function defined from collaborators

    // Apply boundary conditions (*Function outside with .cpp and .h files for each type of BC)
    In the edge loop, if an edge is boundary:
    compute U+ = UL
    build U-=(U+,bc data,n)
    call num_flux(U_plus, U_minus, n_hat)

    for i = 1:nBFaceTot
        if BN(i) == "inflow"
            // Apply inflow boundary condition using inflow_flux function defined from collaborators
            Flux(i, j, bn) = inflow_flux(U(i), Tt, pt, alpha)

            // inflow_flux function: (*Function outside with .cpp and .h files for inflow_flux)
            // function [Fb] = inflow_flux(UI, Tt, pt, n, alpha)
            // now implemented in InflowFlux.cpp and InflowFlux.h files, which calculates the boundary flux given the interior state (UI), total temperature (Tt), total pressure (pt), outward-pointing edge normal (n), and flow angle (alpha) in radians.

        else if BN(i) == "outflow"
            // Apply outflow boundary condition using outflow_flux function defined from collaborators
            Flux(i, j, bn) = outflow_flux(U(i), pout)
            
            / outflow_flux function: (*Function outside with .cpp and .h files for outflow_flux)
            // function [Fb] = outflow_flux(UI, pout)
            // Vb is calculated from the Reimann invariant J+(U+), interior entropy S+ = p+/(rho+)^gamma, and the given outlet static pressure pout. The boundary flux is then calculated using the boundary state (rho_b, Vb, p_b) and the outward normal vector n.
            // rhob = (pb/S+)^(1/gamma)
            // Vb = [ub, vb] = J+ - 2*c_b/(gamma-1) where c_b = sqrt(gamma*pb/rhob)
            // where Vb = Vb - (Vb dot n) * n
            // rhoEb = pb/(gamma-1) + 0.5*rhob*dot(Vb,Vb)
            // Fb = dot(F(Vb),n)

        else if BN(i) == "inviscid wall"
            // Apply wall boundary condition using wall_flux function defined from collaborators
            Flux(i, j, bn) = wall_flux(U(i))

            // wall_flux function: (*Function outside with .cpp and .h files for wall_flux)
            // wall_flux(U) = [0, pb*n_x, pb*n_y, 0]
            // where pb = (gamma - 1) * (rhoE - 0.5 *rho sqrt(ub^2 + vb^2)^2) is the static pressure from the interior state U, and n_x and n_y are the components of the outward normal vector for the boundary edge.
            // Vb = [ub, vb] = Vint - (Vint dot n) * n 

        else if BN(i) == "periodic"
            // Apply periodic boundary condition using periodic_flux function defined from collaborators
            Flux(i, j, bn) = periodic_flux(U(i), U(periodicEdges(i)))
    end for

    // Time-stepping loop (Options for global time stepping or local time stepping)
    t = 0.0; it = 0; Rhist = []

    // Print header for convergence history
    print(f"{'it':>6s} {'t':>10s} {'dt':>10s} {'||R||2':>12s}") 
    
    while (if steady: i not converged 
           else if unsteady: i for each iteration and t < Tfinal)
        
        // Update U using the calculated Residual and Flux arrays
        for i = 1:nElemTot
            for j = 1: edgePerElem
                // Boundary conditions imposing 
                // Precompute boundary states
                UL=U(iL)
                UR=U_ghost(UL, BC data, n) // U_ghost is the ghost state built from the interior state and the boundary condition data (e.g., for wall BC, U_ghost is the mirror state across the wall)
                Flux(i, j, bn) = flux funtion defined from collaborators for boundary edges
                
                // Interior flux calculation
                UL=U(iL)
                UR=U(iR)
                Flux(i, j, in) = flux funtion defined from collaborators for interior edges
                
                // Write one function: (*Function outside with .cpp and .h files for Residual and wave speed calculation)
                // compute_residual_and_wavespeed(mesh, U, t, bc, flux) -> (R, S)
                // where:
                //  residual vector per cell
                //  accumulated wave-speed measure for CFL (e.g., sum over edges of |𝑢⃗⋅𝑛⃗|+𝑐)
                // This is the heart of the solver and should be the only place that loops over edges.
                Residual(i) = Residual(i) + sum(Flux(ui, u_neighbor(i,e), normal(i,e)))* edge_length(i,e) for each edge e of element i
                P(i) = sum(edge_length(i,e) for each edge e of element i) // Perimeter of element i
                d(i) = 2 * Area(i) / P(i) // P(i) is the perimeter of element i
                s(i) = sum of edge(abs(s(i,e))*edge_length(i,e))/P(i)
                deltaT_global(i) = min(CFL * d(i)) / abs(avg(s(i)))
                deltaT_local(i) = (2 Area(i) * CFL) / (s(i) * P(i))
            end for
            if using global time stepping:
                U(i) = U(i) - deltaT_global(i) / Area(i) * Residual(i)
            else if using local time stepping:
                U(i) = U(i) - deltaT_local(i) / Area(i) * Residual(i)
        end for
        
        // Force update of U using the calculated Residual and time step size (local or global)
        // calculate x and y force coefficients on the blade (cx, cy) which are the forces (Fx, Fy) normalized by qout*c with c = 18.804 mm 

        // Print convergence history every 5 iterations or at the final time
        if it % 10 == 0 or t >= Tfinal:
        Rnorm = float(np.linalg.norm(R)); Rhist.append(Rnorm)
        print(f"{it:6d} {t:10.4e} {dt:10.4e} {Rnorm:12.4e}")

        t += dt; it += 1

        // Check for convergence (e.g., based on residual norms or changes in U)
        if (residualNorm < tolerance)
            converged = true
        end if

        // Save solution every save_every iterations (if save_every > 0), and at the final time for unsteady cases
        // if save_every and (it % save_every == 0):
        //     writeU(f"{out_prefix}_{it:07d}.dat", U)

    end while

// Post-processing (*python scripts for post-processing)
// From the converged field:
// 1. Compute Mach and entropy per cell for contour plots
// Compute blade on blade boundary edges using boundary pressure and the given normalization
// cp = (p - pout) / (qout) // pressure coefficient for blade on blade boundary edges
// where qout = 0.5 * gamma * pout *Mout^2
// where Mout^2 = 2/(gamma-1) * ((pout/p0)^((gamma-1)/gamma) - 1)

// 2. calculate x and y force coefficients on the blade (cx, cy) which are the forces (Fx, Fy) normalized by qout*c with c = 18.804 mm 
// Using ground truth on values of cx and cy from a simulation with a finest mesh.

}