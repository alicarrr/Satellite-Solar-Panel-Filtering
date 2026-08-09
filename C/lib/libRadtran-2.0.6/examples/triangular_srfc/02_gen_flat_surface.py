#!/usr/bin/env python3

import xarray as xr
import numpy as np

def setup_surface(srfc_height, only_half_area, materials):
    x=1000
    verts = np.empty((9,3), dtype=np.float64)
    verts[:] = [[  0,   0, srfc_height],
                [x/2,   0, srfc_height],
                [  x,   0, srfc_height],
                [  0, x/2, srfc_height],
                [x/2, x/2, srfc_height],
                [  x, x/2, srfc_height],
                [  0,   x, srfc_height],
                [x/2,   x, srfc_height],
                [  x,   x, srfc_height]]

    # move all building vertices a wee bit towards the center because mystic does strange things at the boundaries
    verts[:,:2] += ( verts[:,:2] - np.array([x, x])/2. ) * -1e-6

    # note here that the different winding order of the triangles is wanted as a test case
    if only_half_area:
        faces = np.empty((4,3), dtype=np.int64)
        faces[0] = [0,1,4]
        faces[1] = [0,3,4]
        faces[2] = [4,5,8]
        faces[3] = [4,7,8]
    else:
        faces = np.empty((2,3), dtype=np.int64)
        faces[0] = [0,2,8]
        faces[1] = [0,6,8]

    material_ids = [ materials.index('GroundSurface') for _ in range(len(faces)) ]

    return verts, faces, material_ids

def _main():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('outfile'        , help="output filename")
    parser.add_argument('-srfc_height'   , type=float, default=0.1, help="srfc height")
    parser.add_argument('-half_area'     , action="store_true", help="if we use a full surface or only half of the area")

    parser.add_argument('-A', type=float, default=0.5, help="surface albedo")
    parser.add_argument('-T', type=float, default=300, help="surface temperature")

    args = parser.parse_args()

    materials = ['GroundSurface',]
    material_albedi = [args.A,]

    verts, faces, material_ids = setup_surface(args.srfc_height, args.half_area, materials)

    temperature_by_material = [args.T,]
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
