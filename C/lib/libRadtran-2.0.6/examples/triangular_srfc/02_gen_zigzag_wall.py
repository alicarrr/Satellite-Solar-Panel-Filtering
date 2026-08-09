#!/usr/bin/env python3

import xarray as xr
import numpy as np
import mesh_utils

def add_cube(verts, faces, material_ids, origin, dx, dy, dz, materials):
    new_verts = np.array(origin, dtype=np.float64) + np.array([
        [0 ,0 , 0], # 0 Ab
        [dx,0 , 0], # 1 Bb
        [0 ,dy, 0], # 2 Cb
        [dx,dy, 0], # 3 Db
        [0 ,0 ,dz], # 4 At
        [dx,0 ,dz], # 5 Bt
        [0 ,dy,dz], # 6 Ct
        [dx,dy,dz], # 7 Dt
        ], dtype=np.float64)

    new_faces = []
    new_faces.append([0,3,1])
    new_faces.append([0,2,3])
    Nground_faces = len(new_faces)
    new_faces.append([0,1,5])
    new_faces.append([0,5,4])
    new_faces.append([0,4,6])
    new_faces.append([0,6,2])
    new_faces.append([1,3,7])
    new_faces.append([1,7,5])
    new_faces.append([2,7,3])
    new_faces.append([2,6,7])
    Nside_faces = len(new_faces) - Nground_faces
    new_faces.append([4,5,7])
    new_faces.append([4,7,6])
    Ntop_faces = len(new_faces) - Nside_faces - Nground_faces
    new_faces = np.array(new_faces, dtype=np.int64)

    new_material_ids = []
    new_material_ids += [ materials.index('GroundSurface') for _ in range(Nground_faces) ]
    new_material_ids += [ materials.index('WallSurface') for _ in range(Nside_faces) ]
    new_material_ids += [ materials.index('RoofSurface') for _ in range(Ntop_faces) ]
    new_material_ids = np.array(new_material_ids, dtype=np.int64)

    new_faces += len(verts)
    verts = np.concatenate((verts, new_verts))
    faces = np.concatenate((faces, new_faces))
    material_ids = np.concatenate((material_ids, new_material_ids))
    return verts, faces, material_ids


def setup_wall(height, dx, dy, N, srfc_height, materials=None):
    """ Generate a zigzag wall from little cubes
        a dx/dy block is split into N substeps
    """

    verts = np.empty((0,3), dtype=np.float64)
    faces = np.empty((0,3), dtype=np.int64)
    material_ids = np.empty((0), dtype=np.int64)

    for i in range(4*N):
        for j in range(4*N):
            sx, sy, sz = dx / float(N), dy / float(N), srfc_height
            origin = [ i*sx, j*sy, 0]
            verts, faces, material_ids = add_cube(verts, faces, material_ids, origin, sx, sy, sz, materials)
            material_ids[:] = materials.index('GroundSurface')

    Nsurf_verts = len(verts)

    is_high_r = np.zeros((4,4), dtype=np.bool)
    is_high_r[0, 3] = True
    is_high_r[1, 0] = True
    is_high_r[2, 1] = True
    is_high_r[3, 2] = True

    is_high_f = np.zeros((4,4), dtype=np.bool)
    is_high_f[0, 2] = True
    is_high_f[1, 3] = True
    is_high_f[2, 0] = True
    is_high_f[3, 1] = True

    is_high_cr = np.zeros((4,4), dtype=np.bool)
    is_high_cr[0, 1] = True
    is_high_cr[1, 2] = True
    is_high_cr[2, 3] = True
    is_high_cr[3, 0] = True

    for x in range(4):
        for y in range(4):
            if is_high_r[y,x]:
                for j in range(N):
                    for i in range(j+1):
                        sx, sy, sz = dx / float(N), dy / float(N), height-srfc_height
                        origin = [ x*dx + i*sx, y*dy + j*sy, srfc_height ]
                        verts, faces, material_ids = add_cube(verts, faces, material_ids, origin, sx, sy, sz, materials)
            if is_high_f[y,x]:
                for i in range(N):
                    for j in range(N):
                        sx, sy, sz = dx / float(N), dy / float(N), height-srfc_height
                        origin = [ x*dx + i*sx, y*dy + j*sy, srfc_height ]
                        verts, faces, material_ids = add_cube(verts, faces, material_ids, origin, sx, sy, sz, materials)

            if is_high_cr[y,x]:
                for j in range(N):
                    for i in range(j+1,N):
                        sx, sy, sz = dx / float(N), dy / float(N), height-srfc_height
                        origin = [ x*dx + i*sx, y*dy + j*sy, srfc_height ]
                        verts, faces, material_ids = add_cube(verts, faces, material_ids, origin, sx, sy, sz, materials)

    # move all building vertices a wee bit towards the center because mystic does strange things at the boundaries
    verts[Nsurf_verts:,:2] += (verts[Nsurf_verts:,:2] - np.array([2*dx, 2*dy]) ) * -1e-4

    return verts, faces, material_ids


def setup_wall2(height, dx, dy, srfc_height=0.1, materials=None):
    """ Generate a zigzag wall
    E.g.

      *----------M--------------------N---------*
      |          |                    |         |
      |  LOW     |                    |         |
      |  GROUND  |                    |         |
      |          |                    |         |
      A----------O          B---------P         Q
      |                     |                   |
      |           HIGH      |                   |
      |           GROUND    |                   |
      |                     |                   |
      |          E----------F         K---------L
      |          |                    |         |
      |          |       LOW          |         |
      |          |       GROUND       |         |
      |          |                    |         |
      C----------D          I---------J         |
      |                     |                   |
      |                     |        HIGH       |
      |                     |        GROUND     |
      |                     |                   |
      *----------R----------G-------------------H

    """
    Nvg = 9 # ground vertices
    Nverts = Nvg + 18*2 # vertices at the wall edges
    verts = np.empty((Nverts,3))
    verts[:] = np.NaN

    e = 1e-3
    A = [0*dx+e, 3*dy+0, srfc_height]; Ab = Nvg+ 0; At = Nvg+18+ 0
    B = [2*dx+0, 3*dy+0, srfc_height]; Bb = Nvg+ 1; Bt = Nvg+18+ 1
    C = [0*dx+e, 1*dy+0, srfc_height]; Cb = Nvg+ 2; Ct = Nvg+18+ 2
    D = [1*dx+0, 1*dy+0, srfc_height]; Db = Nvg+ 3; Dt = Nvg+18+ 3
    E = [1*dx+0, 2*dy+0, srfc_height]; Eb = Nvg+ 4; Et = Nvg+18+ 4
    F = [2*dx+0, 2*dy+0, srfc_height]; Fb = Nvg+ 5; Ft = Nvg+18+ 5
    G = [2*dx+0, 0*dy+e, srfc_height]; Gb = Nvg+ 6; Gt = Nvg+18+ 6
    H = [4*dx-e, 0*dy+e, srfc_height]; Hb = Nvg+ 7; Ht = Nvg+18+ 7
    I = [2*dx+0, 1*dy+0, srfc_height]; Ib = Nvg+ 8; It = Nvg+18+ 8
    J = [3*dx+0, 1*dy+0, srfc_height]; Jb = Nvg+ 9; Jt = Nvg+18+ 9
    K = [3*dx+0, 2*dy+0, srfc_height]; Kb = Nvg+10; Kt = Nvg+18+10
    L = [4*dx-e, 2*dy+0, srfc_height]; Lb = Nvg+11; Lt = Nvg+18+11

    M = [1*dx+0, 4*dy-e, srfc_height]; Mb = Nvg+12; Mt = Nvg+18+12
    N = [3*dx+0, 4*dy-e, srfc_height]; Nb = Nvg+13; Nt = Nvg+18+13
    O = [1*dx+0, 3*dy+0, srfc_height]; Ob = Nvg+14; Ot = Nvg+18+14
    P = [3*dx+0, 3*dy+0, srfc_height]; Pb = Nvg+15; Pt = Nvg+18+15

    Q = [4*dx-e, 3*dy+0, srfc_height]; Qb = Nvg+16; Qt = Nvg+18+16
    R = [1*dx+0, 0*dy+e, srfc_height]; Rb = Nvg+17; Rt = Nvg+18+17

    verts[0:Nvg,:] = [
            [0   , 0   , srfc_height],
            [2*dx, 0   , srfc_height],
            [4*dx, 0   , srfc_height],
            [0   , 2*dy, srfc_height],
            [2*dx, 2*dy, srfc_height],
            [4*dx, 2*dy, srfc_height],
            [0   , 4*dy, srfc_height],
            [2*dx, 4*dy, srfc_height],
            [4*dx, 4*dy, srfc_height],
            ] # ground vertices
    verts[Nvg:Nvg+18] = [ A, B, C, D, E, F, G, H, I, J, K, L, M, N, O, P, Q, R ]
    verts[Nvg+18:Nvg+2*18] = np.array([ A, B, C, D, E, F, G, H, I, J, K, L, M, N, O, P, Q, R ]) + np.array([0,0,height-srfc_height])

    faces = []
    # ground faces
    faces.append([0,1,4])
    faces.append([0,4,3])
    faces.append([1,2,5])
    faces.append([1,5,4])
    faces.append([3,4,7])
    faces.append([3,7,6])
    faces.append([4,5,8])
    faces.append([4,8,7])
    Nground_faces = len(faces)

    # Side walls upper left block
    faces.append([Ab,Ob,Ot]) # walls between [AO]
    faces.append([Ab,Ot,At])
    faces.append([Ab,Ct,Cb]) # between [AC]
    faces.append([Ab,At,Ct])
    faces.append([Cb,Db,Dt]) # [CD]
    faces.append([Cb,Dt,Ct])
    faces.append([Db,Et,Eb]) # [DE]
    faces.append([Db,Dt,Et])
    faces.append([Eb,Fb,Ft]) # [EF]
    faces.append([Eb,Ft,Et])
    faces.append([Fb,Bt,Bb]) # [FB]
    faces.append([Fb,Ft,Bt])

    # Side walls upper left block
    faces.append([Gb,Hb,Ht]) # [GH]
    faces.append([Gb,Ht,Gt])
    faces.append([Hb,Lb,Lt]) # [HL]
    faces.append([Hb,Lt,Ht])
    faces.append([Lb,Kb,Kt]) # [LK]
    faces.append([Lb,Kt,Lt])
    faces.append([Kb,Jb,Jt]) # [KJ]
    faces.append([Kb,Jt,Kt])
    faces.append([Jb,Ib,It]) # [JI]
    faces.append([Jb,It,Jt])
    faces.append([Ib,Gb,Gt]) # [IG]
    faces.append([Ib,Gt,It])

    faces.append([Ob,Ot,Mt]) # [OM]
    faces.append([Ob,Mt,Mb])
    faces.append([Mb,Nt,Nb]) # [MN]
    faces.append([Mb,Mt,Nt])
    faces.append([Pb,Nb,Nt]) # [PN]
    faces.append([Pb,Nt,Pt])
    faces.append([Bb,Pb,Pt]) # [BP]
    faces.append([Bb,Pt,Bt])
    Nside_faces = len(faces) - Nground_faces

    # Top walls
    faces.append([At,Ct,Dt])
    faces.append([At,Dt,Et])
    faces.append([At,Et,Ft])
    faces.append([At,Ft,Bt])
    faces.append([Mt,Ot,Pt])
    faces.append([Mt,Pt,Nt])

    faces.append([Ht,Lt,Kt])
    faces.append([Ht,Kt,Jt])
    faces.append([Ht,Jt,It])
    faces.append([Ht,It,Gt])
    Ntop_faces = len(faces) - Nground_faces - Nside_faces

    material_ids = [ materials.index('GroundSurface') for _ in range(Nground_faces) ]
    material_ids += [ materials.index('WallSurface') for _ in range(Nside_faces) ]
    material_ids += [ materials.index('RoofSurface') for _ in range(Ntop_faces) ]

    return verts, faces, material_ids


def _main():
    import argparse

    parser = argparse.ArgumentParser(description='Generates libRadtran triangle input with a zigzagging wall')
    parser.add_argument('outfile'     , help="output filename")
    parser.add_argument('-dx'         , type=float, default=25 , help="width in x dimension, domain will have size 4*dx")
    parser.add_argument('-dy'         , type=float, default=25 , help="width in y dimension, domain will have size 4*dx")
    parser.add_argument('-height'     , type=float, default=30 , help="height of wall")
    parser.add_argument('-srfc_height', type=float, default=1  , help="srfc height")
    parser.add_argument('-refine'     , type=int  , default=0  , help="number of simple mesh refinement steps")
    parser.add_argument('-N'          , type=int  , default=1  , help="number of zigzag sub-steps")

    parser.add_argument('-AGroundSurface', type=float, default=0.13, help="broadband surface albedo, e.g. weathered road asphalt")
    parser.add_argument('-ARoofSurface'  , type=float, default=0.3, help="broadband roof albedo, e.g. aged red ceramic roofing tiles .24")
    parser.add_argument('-AWallSurface'  , type=float, default=0.23, help="broadband wall albedo, e.g. weathered concrete .23")

    args = parser.parse_args()

    materials = ['GroundSurface', 'RoofSurface', 'WallSurface']
    material_albedi = [args.AGroundSurface, args.ARoofSurface, args.AWallSurface]

    verts, faces, material_ids = setup_wall(args.height, args.dx, args.dy, N=args.N, srfc_height=args.srfc_height, materials=materials)
    for i in range(args.refine):
        verts, faces, material_ids = mesh_utils.mesh_refinement_div4(verts, faces, material_ids)

    temperature_by_material = [275, 280, 285]
    temperatures = np.array([ temperature_by_material[mid] for mid in material_ids ], dtype=np.float64)


    D = xr.Dataset({
        "vertices": xr.DataArray(verts, dims=("Nvert", "Ndim")),
        "triangles": xr.DataArray(faces, dims=("Ntriangles", "Ncorner")),
        "material_of_triangle": xr.DataArray(material_ids, dims=("Ntriangles")),
        "material_albedo": xr.DataArray(material_albedi, dims=("N_materials")),
        "material_type": xr.DataArray(materials, dims=("N_materials")),
        "temperature_of_triangle": xr.DataArray(temperatures, dims=("Ntriangles")),
        })

    D.to_netcdf(args.outfile)

if __name__ == '__main__':
    _main()
