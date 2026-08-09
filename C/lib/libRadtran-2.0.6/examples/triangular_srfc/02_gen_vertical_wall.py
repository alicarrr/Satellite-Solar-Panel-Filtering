#!/usr/bin/env python3

import xarray as xr
import numpy as np

def setup_surface(wall_height, wall_width, origin, x_axis, y_axis, materials):
    x=wall_width
    o=origin
    verts = np.empty((8,3))
    verts[:] = [[o[0]    , o[1]-x/2, o[2]            ],
                [o[0]    , o[1]-x/2, o[2]+wall_height],
                [o[0]-x/2, o[1]    , o[2]            ],
                [o[0]-x/2, o[1]    , o[2]+wall_height],
                [o[0]    , o[1]+x/2, o[2]            ],
                [o[0]    , o[1]+x/2, o[2]+wall_height],
                [o[0]+x/2, o[1]    , o[2]            ],
                [o[0]+x/2, o[1]    , o[2]+wall_height]]

    faces = []
    if x_axis:
        faces.append([2,6,7])
        faces.append([2,7,3])
    if y_axis:
        faces.append([0,4,5])
        faces.append([0,5,1])

    faces = np.array(faces, dtype=np.int64)

    material_ids = [ materials.index('WallSurface') for _ in range(len(faces)) ]

    return verts, faces, material_ids

def _main():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('outfile'        , help="output filename")
    parser.add_argument('-wall_height'   , type=float, default=50, help="wall height")
    parser.add_argument('-wall_width'    , type=float, default=150, help="wall width, centered around origin")
    parser.add_argument('-x_axis'        , type=bool, default=True, help="create a wall along the x axis")
    parser.add_argument('-y_axis'        , type=bool, default=True, help="create a wall along the x axis")
    parser.add_argument('-origin'        , type=float, nargs=3, default=[500,500,0], help="origin of wall")

    parser.add_argument('-AWallSurface', type=float, default=.05, help="wall albedo")

    args = parser.parse_args()

    materials = ['WallSurface',]
    material_albedi = [args.AWallSurface,]

    verts, faces, material_ids = setup_surface(args.wall_height, args.wall_width, args.origin, args.x_axis, args.y_axis, materials)

    temperature_by_material = [310,]
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
