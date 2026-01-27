# !/usr/bin/env python3
# Usage: 
# python pw_to_gri.py \
#   --nodes /home/jaehyn/CFD2/aero-623/projects/Project-1/nodes.txt \
#   --edges /home/jaehyn/CFD2/aero-623/projects/Project-1/edges.txt \
#   --cells /home/jaehyn/CFD2/aero-623/projects/Project-1/cells.txt \
#   -o /home/jaehyn/CFD2/aero-623/projects/Project-1/mesh.gri

import re
import argparse
from collections import defaultdict

def parse_nodes(path):
    coords = {}
    with open(path, "r") as f:
        lines = f.read().splitlines()

    i = 0
    while i < len(lines):
        m = re.match(r"\s*Node\s+(\d+)\s+--", lines[i])
        if m:
            nid = int(m.group(1))
            j = i + 1
            # Find the first line "x y z"
            while j < len(lines):
                parts = lines[j].split()
                if len(parts) >= 3:
                    try:
                        x, y, z = float(parts[0]), float(parts[1]), float(parts[2])
                        coords[nid] = (x, y, z)
                        break
                    except ValueError:
                        pass
                j += 1
            i = j
        i += 1
    if not coords:
        raise RuntimeError(f"No nodes parsed from {path}")
    return coords

def parse_edges(path):
    # returns list of (start, end, ownerCurve) with 1-based node ids
    edges = []
    with open(path, "r") as f:
        lines = f.read().splitlines()

    i = 0
    while i < len(lines):
        m = re.match(r"\s*Edge\s+(\d+)\s+--", lines[i])
        if m:
            start = end = owner = None
            j = i + 1
            while j < len(lines) and (start is None or owner is None):
                # start end length
                m2 = re.search(r"(\d+)\s+(\d+)\s+([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)", lines[j])
                if m2 and start is None:
                    start = int(m2.group(1))
                    end   = int(m2.group(2))
                m3 = re.search(r"Owner:\s*Curve\s+(\d+)", lines[j])
                if m3 and owner is None:
                    owner = int(m3.group(1))
                j += 1
            if start is not None and end is not None and owner is not None:
                edges.append((start, end, owner))
            i = j
        i += 1
    if not edges:
        raise RuntimeError(f"No edges parsed from {path}")
    return edges

def parse_tris(path):
    tris = []
    with open(path, "r") as f:
        for line in f:
            m = re.search(r"Nodal Connectivity:\s*(\d+)\s+(\d+)\s+(\d+)", line)
            if m:
                tris.append((int(m.group(1)), int(m.group(2)), int(m.group(3))))
    if not tris:
        raise RuntimeError(f"No triangles parsed from {path}")
    return tris

def write_gri(outpath, nodes, edges, tris, dim=2):
    # nodes: dict {nid: (x,y,z)}, nids are 1-based
    nmax = max(nodes.keys())
    # Make a dense array 1..nmax
    V = [(0.0, 0.0, 0.0)] * (nmax + 1)
    for nid, xyz in nodes.items():
        V[nid] = xyz

    # Group boundary edges by owner curve
    bcurves = defaultdict(list)
    for s, e, c in edges:
        bcurves[c].append((s, e))

    curve_ids = sorted(bcurves.keys())
    NB = len(curve_ids)

    with open(outpath, "w") as f:
        f.write(f"{nmax} {len(tris)} {dim}\n")
        # vertices: output x y (2D)
        for nid in range(1, nmax + 1):
            x, y, z = V[nid]
            if dim == 2:
                f.write(f"{x:.16e} {y:.16e}\n")
            else:
                f.write(f"{x:.16e} {y:.16e} {z:.16e}\n")

        # boundaries
        f.write(f"{NB}\n")
        for cid in curve_ids:
            bedges = bcurves[cid]
            # readgri() uses token[0] as Nb and token[2] as name
            f.write(f"{len(bedges)} 2 Curve{cid}\n")
            for s, e in bedges:
                f.write(f"{s} {e}\n")

        # elements - single block
        f.write(f"{len(tris)} 1 TriLagrange\n")
        for a, b, c in tris:
            f.write(f"{a} {b} {c}\n")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nodes", required=True)
    ap.add_argument("--edges", required=True)
    ap.add_argument("--cells", required=True)
    ap.add_argument("-o", "--out", required=True)
    ap.add_argument("--dim", type=int, default=2)
    args = ap.parse_args()

    nodes = parse_nodes(args.nodes)
    edges = parse_edges(args.edges)
    tris  = parse_tris(args.cells)

    write_gri(args.out, nodes, edges, tris, dim=args.dim)
    print(f"Wrote {args.out}")

if __name__ == "__main__":
    main()
