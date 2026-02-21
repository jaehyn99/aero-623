import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import os

#-----------------------------------------------------------
def readgri(fname):
    with open(fname, 'r') as f:
        Nn, Ne, dim = [int(s) for s in f.readline().split()]
        if dim != 2:
            raise ValueError(f"Expected dim=2, got {dim}")

        V = np.array([[float(s) for s in f.readline().split()] for _ in range(Nn)])

        NB = int(f.readline())
        B = []; Bname = []
        for _ in range(NB):
            s = f.readline().split()
            Nb = int(s[0]); Bname.append(s[2])
            Bi = np.array([[int(t)-1 for t in f.readline().split()] for _ in range(Nb)])
            B.append(Bi)

        Ne0 = 0; E = []
        while Ne0 < Ne:
            s = f.readline().split()
            ne_blk = int(s[0])
            Ei = np.array([[int(t)-1 for t in f.readline().split()] for _ in range(ne_blk)])
            E = Ei if (Ne0 == 0) else np.concatenate((E, Ei), axis=0)
            Ne0 += ne_blk

    return {'V': V, 'E': E, 'B': B, 'Bname': Bname}

def read_xy_txt(path: str):
    arr = np.loadtxt(path, dtype=float)
    if arr.ndim == 1:
        arr = arr.reshape(1, -1)
    return arr[:, :2]

#-----------------------------------------------------------
def read_scalar_field_maybe_header(path):
    """
    Accepts either:
      A) N lines of floats
      B) first line is integer N, followed by N floats
    Returns a 1D float array.
    """
    vals = np.loadtxt(path, dtype=float, ndmin=1)
    if vals.size == 0:
        return vals

    # If first entry is an integer-like count and matches remaining length, drop it.
    n0 = int(round(vals[0]))
    if abs(vals[0] - n0) < 1e-12 and (vals.size - 1) == n0:
        return vals[1:]
    return vals

#-----------------------------------------------------------
def plot_wall_distance(Mesh, dist, fname, show_mesh=True, use_log=False, clim=None, plot_sizing=False):
    V = Mesh['V']; E = Mesh['E']

    if dist.shape[0] != V.shape[0]:
        raise ValueError(f"dist length {dist.shape[0]} != number of nodes {V.shape[0]}")

    if plot_sizing:
        if use_log:
            field = np.log10(np.maximum(field, 1e-16))
            title = "log10(size function h)"
        else:
            field = dist.copy()
            title = "Size function h"
    else:
        if use_log:
            field = np.log10(np.maximum(field, 1e-16))
            title = "log10(wall distance)"
        else:
            field = dist.copy()
            title = "Wall distance"
    

    fig = plt.figure(figsize=(12, 8))
    tpc = plt.tripcolor(V[:, 0], V[:, 1], E, field, shading='gouraud')

    if clim is not None:
        tpc.set_clim(clim[0], clim[1])

    plt.colorbar(tpc, shrink=0.85, label=title)

    if show_mesh:
        plt.triplot(V[:, 0], V[:, 1], E, 'k-', linewidth=0.2)

    plt.axis('equal')
    plt.tight_layout()
    plt.savefig(fname, dpi=300)
    plt.close(fig)

def plot_mesh_with_blades(mesh, blade_upper, blade_lower, out_png):
    V = mesh["V"]; E = mesh["E"]

    fig = plt.figure(figsize=(10, 10))
    plt.triplot(V[:,0], V[:,1], E, linewidth=0.3, color='black', alpha=1)
    up = read_xy_txt(blade_upper)
    lo = read_xy_txt(blade_lower)
    # plt.plot(up[:,0], up[:,1], linewidth=2.0)
    # plt.plot(lo[:,0], lo[:,1], linewidth=2.0)
    plt.legend(loc="best")
    plt.axis("equal")
    plt.tight_layout()
    plt.savefig(out_png, dpi=400)
    plt.close(fig)

#-----------------------------------------------------------
def main():

    base = os.getcwd()

    # gri_file  = base + "/mesh_refined.gri"
    # dist_file = gri_file.replace(".gri", ".walldist.txt")
    # out_png_dist   = base + "/wall_distance.png"
    # out_png_size = base + "/size_function.png"

    # Mesh = readgri(gri_file)
    # dist = read_scalar_field_maybe_header(dist_file)

    # print("Nnodes =", Mesh["V"].shape[0])
    # print("Ndist  =", dist.shape[0])
    # if dist.size:
    #     print("min/max=", dist.min(), dist.max())

    # plot_wall_distance(Mesh, dist, out_png_dist, show_mesh=True, use_log=False)

    # h = np.loadtxt(base + "/mesh_refined.hnode.txt")
    # plot_wall_distance(Mesh, h, out_png_size, show_mesh=True, use_log=False, plot_sizing=True)
    
    mesh_file = base + "/projects/Project-1/mesh_coarse.gri"  # available in sandbox
    mesh = readgri(str(mesh_file))

    out_png = base + "/mesh_coarse_overlay_blades.png"
    plot_mesh_with_blades(mesh, str(base + "/projects/Project-1/bladeupper.txt"), str(base + "/projects/Project-1/bladelower.txt"), str(out_png))


if __name__ == "__main__":
    main()
