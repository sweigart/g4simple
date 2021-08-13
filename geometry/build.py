# usage: python build.py geometry.json

from functools import wraps
import json
from datetime import date
from materials import *
import pandas as pd
import numpy as np
import math
import sys
import os
from operator import add
import ROOT as r
r.gSystem.Load("libGeom")
r.gSystem.Load("libPhysics.so")

# here are some necessary global lists to keep track of things (and keep them in memory)
# all materials are appended here to keep track of them until they are output
all_root_materials = []
# this is used to keep track of which detectors should have material names in the Auxiliary information
material_dict = {}
targets = []            # all geometries get appended to here
defined_materials = {}  # we can re-use materials with the same name
additional_material_properties = {}

# create geometry manager
manager = r.TGeoManager("manager", "Geometry Manager for PIENUX Detector Creation")

geometry_path = os.environ["GEOMETRY_BUILD_DIR"]

def get_and_register_material(name, debug=False):
    '''
        Because we want to ensure that the materials and mediums are kept in memory, we need to append them to 
        all_root_materials when they are created. This is a little bit of a hack to make sure Pythons memory management
        doesn't make us pay later.
    '''
    if name not in defined_materials:
        this_medium, this_material, props = detector_materials(name, debug)
        # this_material.RegisterYourself()
        all_root_materials.append(this_material)
        all_root_materials.append(this_medium)
        defined_materials[name] = (this_medium, this_material)
        if props is not None:
            additional_material_properties[this_material.GetName()] = props
    else:
        this_medium = defined_materials[name][0]
    return this_medium

def build_world(params):
    this_material = get_and_register_material(params['material'])
    world = r.gGeoManager.MakeBox("world", this_material, params['xWidth'], params['yWidth'], params['zWidth'])
    world.SetVisibility(False)
    world.SetLineColor(1)
    targets.append(world)

    return world

def build_calorimeter(params, world):
    """
    The cool and only slightly broken crystal version of the calorimeter. The
    issue with it right now is overlap in the crystals, particularly at low
    frequencies in pyDome. Adding value to the trim parameter is a bandaid
    solution on this for now
    """

    # gonna use these a lot later on, so make them more reader friendly
    cos = np.cos
    sin = np.sin

    # useful to shorten
    rInn = params['rMin']
    rOut = params['rMax']
    beamAngle = params['thetaMin'] * np.pi / 180

    """
    Do something here to get whatever geometry you want for the inner surface
    """
    """
    Now you have a wrl file with the inner surface geometry you want; you
    need to read the vertices and put the coordinates into a fun format.
    """

    file = open('./pyDomeOutput.wrl', 'r')
    lines = file.readlines()

    # df to hold coordinates of each vertex
    vertex_df = pd.DataFrame(columns=['x', 'y', 'z', 'rho', 'theta', 'phi'])

    coordIndex = []

    index = 0
    # get our xyz points
    for line in lines[13:]:    # start at 13 to skip the header
        this_line = line[:-2].split()    # get rid of ,\n at the ends
        try:
            # grab our cool and good coordinates
            x = float(this_line[0])
            y = float(this_line[1])
            z = float(this_line[2])

            rho = np.sqrt((x**2) + (y**2) + (z**2))
            theta = np.arccos(z/rho)
            try:
                omega = y/x
                if x > 0:  # Q1 and Q4
                    phi = np.arctan(omega)
                elif y >= 0 and x < 0:  # Q2
                    phi = np.pi + np.arctan(omega)
                elif y < 0 and x < 0:   # Q3
                    phi = np.arctan(omega) - np.pi
            except(ZeroDivisionError):
                if y > 0:
                    phi = np.pi/2
                else:
                    phi = -np.pi/2

            vertex_df.loc[index] = {'x': x, 'y': y, 'z': z,
                                    'rho': rho, 'theta': theta, 'phi': phi}
            index += 1
        except(ValueError):   # reached the end of the xyz points
            # we need to handle this line differently, it has a bunch of junk
            this_line = line[line.index('[')+1:].split()
            miniCoordIndex = []
            for item in this_line[:-1]:
                miniCoordIndex.append(int(item[0:-1]))
            coordIndex.append(miniCoordIndex)
            # we are done in this loop, but need the same index for the next,
            # so we still want to iterate it
            index += 1
            break
    # now we are fully in the coordIndex space of the file
    for line in lines[13+index:]:
        this_line = line.split()
        miniCoordIndex = []    # a single line of coordIndex, as a list
        for item in this_line[:-1]:
            try:
                miniCoordIndex.append(int(item[0:-1]))
            except(ValueError):  # the last line is just formatting
                break
        if miniCoordIndex:
            coordIndex.append(miniCoordIndex)
        index += 1

    """
    Ok, we have our vertices, now we want to start making crystals. The overall
    process is, we have vertices that define a surface. We then rotate that
    surface so that its center axis is along the z axis. Root then creates a
    crystal volume for us by creating an outer surface further along z, that is
    a little bigger based on how many crystals we have in our sphere. We then
    rotate the crystal back to the angle of its original center, and continue
    that process until we have built the entire sphere. For a geometry where
    all crystals are the same this is a little inefficient, but it is general
    and should work for sphere-like geometries where that isn't the case.
    """

    # coordIndex tells us which coordinates should connect to each other, so
    # we want to go through every vertex grouping in coordIndex and turn the
    # surface defined by those vertices into the inner surface of the crystal
    for sIndex in range(0, len(coordIndex)):
        # get all our vertices

        # this will be used if one of our vertices is within the area reserved
        # for the beam. If it turns true, we will skip this crystal
        beamFlag = False

        # this is where we keep the points on the inner surface. We don't
        # *need* to put xyz in separate lists like this, but it makes the
        # computation for finding the center vector later a little easier
        surface_x = []
        surface_y = []
        surface_z = []
        surface_vectors = []

        # build those lists
        for pIndex in coordIndex[sIndex]:
            x = vertex_df.loc[pIndex]['x']
            y = vertex_df.loc[pIndex]['y']
            z = vertex_df.loc[pIndex]['z']

            surface_x.append(x)
            surface_y.append(y)
            surface_z.append(z)

            # check if this vertex is in conflict with the beam angle or not
            theta = np.arccos(z/np.sqrt((x**2)+(y**2)+(z**2)))
            # the beam angle is on both sides, so we need to check there
            if (theta <= beamAngle) or (np.pi - theta <= beamAngle):
                if params['cutBeam']:
                    beamFlag = True
                    break

            surface_vectors.append(np.array([x, y, z]))

        # skip this crystal if we are blocking the beam
        if beamFlag:
            continue

        # find the vector normal to this surface at its center
        center = np.array([np.average(surface_x),
                           np.average(surface_y),
                           np.average(surface_z)])

        # define our rotation angles onto the z axis

        # first get the normal vector into radial coordinates
        norm_r = np.sqrt((center[0] ** 2) +
                         (center[1] ** 2) +
                         (center[2] ** 2))
        norm_theta = np.arccos(center[2]/norm_r)
        norm_phi = np.arctan(center[1]/center[0])

        # angle to rotate about Z
        A = -norm_phi

        # angle to rotate about Y
        if center[0] > 0:
            B = -norm_theta
        else:
            B = norm_theta

        # the rotation matrix to put our center vector onto z
        rot = np.array((
            (cos(B)*cos(A), -cos(B)*sin(A), sin(B)),
            (sin(A),         cos(A),        0),
            (-sin(B)*cos(A), sin(B)*sin(A), cos(B))
        ))

        # puts the vector back where it was
        invRot = np.linalg.inv(rot)

        # a place to put rotated coordinates. We don't need z, because
        # everything needs to be flat on z for root's geometry creation
        x_prime_inner = []
        y_prime_inner = []

        # do the transformation
        for vector in surface_vectors:
            transformed_vector = np.matmul(rot, vector)
            x_prime_inner.append(transformed_vector[0])
            y_prime_inner.append(transformed_vector[1])
            if params['debug']:
                print("T Vector: {0}".format(transformed_vector))

        # root wants arrays, not lists
        x_prime_inner = np.array(x_prime_inner)
        y_prime_inner = np.array(y_prime_inner)

        # alpha lets us convert between radius of the sphere and distance to
        # the center of our surfaces. We need this because we are building our
        # crystal along z, where the surface center vector is right now
        alpha = np.arccos(norm_r/rInn)

        if params['wrapping']:
            # we need to know how much to scale the crystal down by
            wrapRatioInnC = 1 - \
                ((params['wrappingThickness'] + 2 *
                 params['trim']) / (rInn*sin(alpha)))
            wrapRatioOutC = 1 - \
                ((params['wrappingThickness'] + 2 *
                 params['trim']) / (rOut*sin(alpha)))
            # outer
            wrapRatioInn = 1-(params['trim']) / (rInn*sin(alpha))
            wrapRatioOut = 1-(params['trim']) / (rOut*sin(alpha))
            # inner
            wrapRatioInn2 = 1 - \
                ((params['wrappingThickness'] + 1 *
                 params['trim']) / (rInn*sin(alpha)))
            wrapRatioOut2 = 1 - \
                ((params['wrappingThickness'] + 1 *
                 params['trim']) / (rOut*sin(alpha)))

            # create the crystal
            this_crystal = r.gGeoManager.MakeXtru("xtal"+str(sIndex),
                                                  get_and_register_material(params['xtalMaterial']), 2)
            this_crystal.GetShape().DefinePolygon(3, x_prime_inner,
                                                  y_prime_inner)

            # the crystal is a little smaller than without wrapping
            this_crystal.GetShape().DefineSection(0, norm_r,
                                                  0, 0, 1*wrapRatioInnC)
            this_crystal.GetShape().DefineSection(1, rOut*cos(alpha),
                                                  0, 0, rOut/rInn*wrapRatioOutC)
            # this_crystal.SetLineColor(2)
            this_crystal.RegisterYourself()

            # we also now need a wrapping volume
            outer = r.gGeoManager.MakeXtru("outer"+str(sIndex),
                                           get_and_register_material(params['wrapMaterial']), 2)
            outer.GetShape().DefinePolygon(3, x_prime_inner, y_prime_inner)
            outer.GetShape().DefineSection(0, norm_r, 0, 0, 1*wrapRatioInn)
            outer.GetShape().DefineSection(1, rOut*cos(alpha),
                                           0, 0, rOut/rInn*wrapRatioOut)
            # outer.RegisterYourself()

            inner = r.gGeoManager.MakeXtru("inner"+str(sIndex),
                                           get_and_register_material(params['wrapMaterial']), 2)
            inner.GetShape().DefinePolygon(3, x_prime_inner, y_prime_inner)
            inner.GetShape().DefineSection(0, norm_r, 0, 0, 1*wrapRatioInn2)
            inner.GetShape().DefineSection(1, rOut*cos(alpha),
                                           0, 0, rOut/rInn*wrapRatioOut2)
            # inner.RegisterYourself()

            addString = "(outer"+str(sIndex) + ") - (" + \
                "inner"+str(sIndex)+")"
            wrapping = r.TGeoCompositeShape("wrap"+str(sIndex), addString)

            wrapVolume = r.TGeoVolume("wrapVolume"+str(sIndex),
                                      wrapping, inner.GetMedium())
            # wrapVolume.SetLineColor(1)
            wrapVolume.RegisterYourself()

            if params['debug']:
                wrapVolume.CheckOverlaps()

        else:
            # create the crystal
            this_crystal = r.gGeoManager.MakeXtru("xtal"+str(sIndex),
                                                  get_and_register_material(params['xtalMaterial']), 2)
            this_crystal.GetShape().DefinePolygon(3, x_prime_inner,
                                                  y_prime_inner)
            this_crystal.GetShape().DefineSection(0, norm_r, 0, 0, 1)
            this_crystal.GetShape().DefineSection(1, rOut*cos(alpha),
                                                  0, 0, (rOut)/rInn)

            this_crystal.RegisterYourself()

        # translate/rotate the crystal as needed to put it back

        # where should the center of the crystal, not just the inner surface,
        # be now that we have built our crystal and want to put it back?
        # Remember root will move the center of the crystal volume, so we need
        # to change our definition of center vector accordingly
        rotCenter = np.matmul(rot, center)
        newCenter = np.matmul(invRot, rotCenter - np.array([0, 0, norm_r]))

        this_translation = r.TGeoTranslation(newCenter[0],
                                             newCenter[1],
                                             newCenter[2])

        # now we need angles to rotate the crystal

        if center[0] >= 0:
            backTheta = norm_theta
        else:
            backTheta = -norm_theta

        backPhi = norm_phi + np.pi/2

        this_rotation = r.TGeoRotation("rot_"+str(sIndex),
                                       backPhi * 180 / np.pi,
                                       backTheta * 180 / np.pi,
                                       -90)

        # put the crystal in the crystal home, and apply the combined movements

        combined = r.TGeoCombiTrans(this_translation, this_rotation)
        world.AddNode(this_crystal, 300000+sIndex, combined)

        if params['wrapping']:
            all_root_materials.append(wrapping)
            all_root_materials.append(wrapVolume)

            #TODO: Clean this up and make it more general
            wrapOptical = r.TGeoOpticalSurface(
                f"wrapFinish_{350000+sIndex}", 
                r.TGeoOpticalSurface.ESurfaceModel.kMunified, 
                r.TGeoOpticalSurface.ESurfaceFinish.kFground)
            # wrapOptical.AddProperty("REFLECTIVITY", 'REF_Wrapping')
            # manager.AddProperty("REF_Wrapping", 0)
            manager.AddOpticalSurface(wrapOptical)
            wrapSkin = r.TGeoSkinSurface(f"wrapSkin_{sIndex + 350000}", wrapOptical.GetName(), wrapOptical, wrapVolume)
            # wrapSkin.Print()
            manager.AddSkinSurface(wrapSkin)

            all_root_materials.append(wrapSkin)
            all_root_materials.append(wrapOptical)

            world.AddNode(wrapVolume, 350000+sIndex, combined)

            world.CheckOverlaps()

def main():
    config_file = sys.argv[1]
    print("Loading geometry configuration from:", config_file)
    with open(config_file, "r") as f:
        params = json.load(f)

    print("Parameters", params)

    print("Building world....")
    world = build_world(params['world'])
    r.gGeoManager.SetTopVolume(world)
    r.gGeoManager.SetVisOption(0)

    if("calorimeter" in params):
        print("Building calorimeter...")
        calorimeter = build_calorimeter(params['calorimeter'], world)

    print("Writing to file...")

    print("Top volume:", r.gGeoManager.GetTopVolume())
    if(params['debug']):
        r.gGeoManager.GetTopVolume().Draw("")
        input("so?")

    outfile = params['filename']
    export_gdml_with_aux(outfile + ".gdml", material_dict, params)

def export_gdml_with_aux(outfile, material_dict, params):
    '''
        ROOT doesn't export any of the auxiliary info by default, so we can do it ourselves.
        material_dict : dict of names of materials to what kind of detector they should be, if not specified then defaults not-sensitive
    '''
    temp = './temp.gdml'
    # r.gGeoManager.Export(temp,"", "vgf")
    r.gGeoManager.Export(temp)
    print("Adding auxiliary information to GDML file...")
    lines = []
    with open(temp, 'r') as f:
        for line in f:
            lines.append(line)
            # add auxilliary information
            spaces = line.split("<")[0]  # +"  "
            # if("<define/>" in line):
            #     lines.pop() # remove the
            if("<define" in line):
                # add material definitions from additional_material_properties
                for name in additional_material_properties:
                    for x in additional_material_properties[name]:
                        print(x)
                        lines.append(x.format())
            if("<material " in line):
                # add link to material definition from additional_material_properties
                name = line.split('name="')[1].split('"')[0]
                if(name.strip() in additional_material_properties):
                    for x in additional_material_properties[name]:
                        lines.append(x.property())
            if('<materialref' in line):
                name = line.split('ref="')[1].split('"')[0]
                print(name)
                try:
                    detector_type = material_dict[name.strip().split("0x", 1)[
                        0]]
                except:
                    detector_type = 'Unspecified'
                auxinfo2 = f'<auxiliary auxtype="SensDet" auxvalue="{detector_type}"/> '+'\n'
                auxinfo = f'<auxiliary auxtype="Material" auxvalue="{name}"/> '+'\n'
                lines.append(spaces+auxinfo)
                lines.append(spaces+auxinfo2)
            if('<solidref' in line):
                name = line.split('ref="')[1].split('"')[0]
                print(name)
                auxinfo = f'<auxiliary auxtype="Solid" auxvalue="{name}"/> '+'\n'
                lines.append(spaces+auxinfo)
            # if('<physvol' in line):
            #     name = line.split('name="')[1].split('"')[0]
            #     copy_number = int(line.split('copynumber="')[1].split('"')[0])
                if('pixel' in name):
                    thislayer = params['target']['nLayers']
                    thisx = params['target']['xPixels']
                    thisy = params['target']['yPixels']
                    auxinfo = f'<auxiliary auxtype="nLayers" auxvalue="{thislayer}"/> '+'\n'
                    auxinfo2 = f'<auxiliary auxtype="nXPixel" auxvalue="{thisx}"/> '+'\n'
                    auxinfo3 = f'<auxiliary auxtype="nYPixel" auxvalue="{thisy}"/> '+'\n'
                    lines.append(spaces+auxinfo)
                    lines.append(spaces+auxinfo2)
                    lines.append(spaces+auxinfo3)

    # write the params in a comment in the GDML code
    if(params is not None):
        lines.append('\n<!--\n')
        lines.append(f'Generated on {date.today()} using: \n')
        for line in json.dumps(params).split("\n"):
            lines.append(line+'\n')
        lines.append('-->\n')

    with open(outfile, 'w') as f:
        for line in lines:
            f.write(line)
    os.remove(temp)

if __name__ == "__main__":
    main()
