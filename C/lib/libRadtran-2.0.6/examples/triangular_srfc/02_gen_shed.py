#!/usr/bin/env python3

import xarray as xr
import numpy as np

def setup_building(h, w, l, overhang, angle, srfc_height, translation, materials):
    roof_height = np.tan(np.deg2rad(angle)) * w/2.
    roof_overhang_height = np.tan(np.deg2rad(angle)) * overhang
    x=1000
    verts = np.empty((20,3))
    verts[0:4,:] = [[ 0,  0, srfc_height], [x,  0, srfc_height], [ 0, x, srfc_height], [x, x, srfc_height]] # 4 ground vertices

    # 8 base building corners
    verts[4]  = [0,0,0]
    verts[5]  = [w,0,0]
    verts[6]  = [0,0,h]
    verts[7]  = [w,0,h]
    verts[8]  = [0,l,0]
    verts[9]  = [w,l,0]
    verts[10] = [0,l,h]
    verts[11] = [w,l,h]
    # plus 2 verts at the gable of the roof, inside
    verts[12] = [w/2.,0,h+roof_height]
    verts[13] = [w/2.,l,h+roof_height]
    # 4 corners of roof
    verts[14] = [-overhang,-overhang,h-roof_overhang_height]
    verts[15] = [w+overhang,-overhang,h-roof_overhang_height]
    verts[16] = [-overhang,l+overhang,h-roof_overhang_height]
    verts[17] = [w+overhang,l+overhang,h-roof_overhang_height]

    # plus 2 verts at the gable of the roof, outside
    verts[18] = [w/2.,-overhang ,h+roof_height]
    verts[19] = [w/2.,l+overhang,h+roof_height]

    verts[4:] += np.array(translation)

    faces = np.empty((20,3), dtype=np.int64)
    faces[0] = [0,1,3] # 2 ground faces
    faces[1] = [0,3,2]
    faces[2] = [4,5,7] # 2 front faces
    faces[3] = [4,7,6]
    faces[4] = [8,11,9] # 2 back faces
    faces[5] = [8,10,11]
    faces[6] = [5,9,11] # 4 side faces
    faces[7] = [5,11,7]
    faces[8] = [4,6,10]
    faces[9] = [4,10,8]
    faces[10] = [4,8,9] # 4 bot/top faces
    faces[11] = [4,9,5]
    faces[12] = [6,10,11]
    faces[13] = [6,11,7]
    faces[14] = [6,7,12] # 2 gable front and rear of base wall at roof height
    faces[15] = [10,13,11]
    faces[16] = [14,18,19] # 4 roof faces
    faces[17] = [14,19,16]
    faces[18] = [15,17,19]
    faces[19] = [15,19,18]

    material_ids = [ materials.index('GroundSurface') for _ in range(2) ]
    material_ids += [ materials.index('WallSurface') for _ in range(2,16) ]
    material_ids += [ materials.index('RoofSurface') for _ in range(16,20) ]

    return verts, faces, material_ids

def _main():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('outfile'        , help="output filename")
    parser.add_argument('-height'        , type=float, default=8 , help="height of the core building (without roof)")
    parser.add_argument('-width'         , type=float, default=10, help="width of the building along x-axis")
    parser.add_argument('-length'        , type=float, default=20, help="length of the building along y-axis")
    parser.add_argument('-roof_overshoot', type=float, default=2.5 , help="excess overshoot of the roof")
    parser.add_argument('-roof_angle'    , type=float, default=30 , help="angle of the root")
    parser.add_argument('-srfc_height'   , type=float, default=0.01 , help="srfc height")
    parser.add_argument('-translation'   , type=float, nargs=3, default=[500,500,0] , help="translation of building 3 arguments [x,y,z]")

    parser.add_argument('-AGroundSurface', type=float, default=0.1, help="surface albedo")
    parser.add_argument('-ARoofSurface'  , type=float, default=0.3, help="roof albedo")
    parser.add_argument('-AWallSurface'  , type=float, default=0.7, help="wall albedo")

    args = parser.parse_args()


    materials = ['GroundSurface', 'RoofSurface', 'WallSurface']
    material_albedi = [args.AGroundSurface, args.ARoofSurface, args.AWallSurface]

    verts, faces, material_ids = setup_building(args.height, args.width, args.length, args.roof_overshoot, args.roof_angle, args.srfc_height, args.translation, materials)

    temperature_by_material = [270, 300, 285]
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
