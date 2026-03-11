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
    plt.legend(loc="best")
    plt.axis("equal")
    plt.tight_layout()
    plt.savefig(out_png, dpi=400)
    plt.close(fig)

# -----------------------------------------------------------
def read_curved_edges_txt(path: str):
    """
    Reads files of the form:

    Curve5 44 45
    x0 y0
    x1 y1
    x2 y2

    Curve5 45 46
    x0 y0
    x1 y1
    x2 y2
    ...

    Returns a list of dicts:
      [{"curve": str, "n1": int, "n2": int, "pts": np.ndarray}, ...]
    """
    edges = []
    with open(path, "r") as f:
        lines = [ln.strip() for ln in f if ln.strip()]

    i = 0
    while i < len(lines):
        header = lines[i].split()
        if len(header) != 3:
            raise ValueError(f"Bad header line in {path}: {lines[i]}")
        curve_name = header[0]
        n1 = int(header[1])
        n2 = int(header[2])
        i += 1

        pts = []
        while i < len(lines):
            parts = lines[i].split()
            if len(parts) == 3:
                # next header
                break
            if len(parts) != 2:
                raise ValueError(f"Bad point line in {path}: {lines[i]}")
            pts.append([float(parts[0]), float(parts[1])])
            i += 1

        edges.append({
            "curve": curve_name,
            "n1": n1,
            "n2": n2,
            "pts": np.array(pts, dtype=float)
        })

    return edges

# -----------------------------------------------------------
def plot_mesh_with_blades_and_lagrange_points(
    mesh,
    blade_upper,
    blade_lower,
    upper_curved_file,
    lower_curved_file,
    out_png
):
    V = mesh["V"]
    E = mesh["E"]

    up = read_xy_txt(blade_upper)
    lo = read_xy_txt(blade_lower)

    upper_edges = read_curved_edges_txt(upper_curved_file)
    lower_edges = read_curved_edges_txt(lower_curved_file)

    fig = plt.figure(figsize=(10, 10))

    # mesh
    plt.triplot(V[:, 0], V[:, 1], E, linewidth=0.3, color='black', alpha=0.6)

    # # original blade curves
    # plt.plot(up[:, 0], up[:, 1], '-', linewidth=1.5, label='upper blade curve')
    # plt.plot(lo[:, 0], lo[:, 1] + 18.0, '-', linewidth=1.5, label='lower blade curve')

    # curved edge polylines
    for k, edge in enumerate(upper_edges):
        pts = edge["pts"]
        plt.plot(
            pts[:,0], pts[:,1],
            's',
            markersize=2.5,
            markerfacecolor='none',
            markeredgecolor='black',
            linestyle='None',
            label='upper q-nodes' if k == 0 else None
        )

    for k, edge in enumerate(lower_edges):
        pts = edge["pts"]
        plt.plot(
            pts[:,0], pts[:,1],
            's',
            markersize=2.5,
            markerfacecolor='none',
            markeredgecolor='black',
            linestyle='None',
            label='lower q-nodes' if k == 0 else None
        )

    plt.text(
        0.98, 0.98,
        r"$Q = 1$",
        transform=plt.gca().transAxes,
        ha="right",
        va="top",
        fontsize=14,
        bbox=dict(facecolor="white", edgecolor="none", alpha=0.7)
    )
    
    plt.axis("equal")
    plt.legend(loc="best")
    plt.tight_layout()
    plt.savefig(out_png, dpi=400)
    plt.close(fig)


# -----------------------------------------------------------
def main():
    base = os.getcwd()

    mesh_file = base + "/mesh_coarse.gri"
    mesh = readgri(mesh_file)

    blade_upper = base + "/bladeupper.txt"
    blade_lower = base + "/bladelower.txt"

    upper_curved_file = base + "/upper_curved_edges_Q1_coarse.txt"
    lower_curved_file = base + "/lower_curved_edges_Q1_coarse.txt"
    # upper_curved_file = base + "/upper_curved_edges_Q2_coarse.txt"
    # lower_curved_file = base + "/lower_curved_edges_Q2_coarse.txt"
    # upper_curved_file = base + "/upper_curved_edges_Q3_coarse.txt"
    # lower_curved_file = base + "/lower_curved_edges_Q3_coarse.txt"

    out_png = base + "/mesh_coarse_overlay_blades_lagrange.png"

    plot_mesh_with_blades_and_lagrange_points(
        mesh,
        blade_upper,
        blade_lower,
        upper_curved_file,
        lower_curved_file,
        out_png
    )


if __name__ == "__main__":
    main()