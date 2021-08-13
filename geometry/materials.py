import json
import math
import sys
import os
import numpy as np
from operator import add
import ROOT as r
r.gSystem.Load("libGeom")
r.gSystem.Load("libPhysics.so")


class MaterialProperty():
    def __init__(self, typestr, firstUnit, data, secondUnit="", name=""):
        self.type = typestr
        self.name = typestr+'_'+name
        self.ncol = len(data[0]) if (type(data[0]) in [list, tuple]) else 1
        self.firstUnit = firstUnit
        self.secondUnit = secondUnit
        self.units = [str(self.firstUnit), str(self.secondUnit)]
        self.data = data

    def format(self):
        if(self.ncol > 0):
            return self.format_2()
        elif self.ncol == 1:
            return self.format_1()
        else:
            raise NotImplementedError

    def format_1(self):
        return f'''<matrix name="{self.name}"  coldim="1" values="{self.data[0][0]}*{self.firstUnit}" />\n'''

    def format_2(self):
        thisstr = f'''<matrix name="{self.name}" coldim="{self.ncol}" values="'''
        for val in self.data:
            for i in range(self.ncol):
                thisstr += f'{val[i]}'
                thisstr += f'*{self.units[i]} ' if (
                    len(self.units[i]) > 0) else " "
            thisstr += '\n'
        thisstr = thisstr[:-2]
        thisstr += '''"/>\n'''
        # print(thisstr)
        return thisstr

    def property(self):
        return f'''<property name="{self.type}" ref="{self.name}"/>\n'''


def air(name="AIR"):
    # mat_AIR = r.TGeoMixture(name, 4, .00000120479)         # Define a mixture
    mat_AIR = r.TGeoMixture(name, 4, 1e-25)                  # Define a mixture
    # Carbon Dioxide   DefineElement(index, mass (amu), Z, percentage composition)
    mat_AIR.DefineElement(0, 12.01115, 6, 0.000124)
    mat_AIR.DefineElement(1, 14.0067, 7, 0.755268)			 # Nitrogen
    mat_AIR.DefineElement(2, 15.9994, 8, 0.231781)			 # Oxygen
    mat_AIR.DefineElement(3, 39.948, 18, 0.012827)			 # Argon
    med_AIR = r.TGeoMedium(name, 2, mat_AIR)                 # TGeoMedium
    # all_root_materials.append(mat_AIR)
    # all_root_materials.append(med_AIR)

    air_Energy = [0.000070, 7.07, 7000.14]
    air_SCINT = [0.1, 1.0, 0.1]
    air_RIND = [1.00, 1.00, 1.00]
    air_ABSL = [3500., 3500., 3500.]

    props = [
        # MaterialProperty("SCINTILLATIONCOMPONENT1", 'eV', [x for x in zip(air_Energy, air_SCINT)], name=name),
        # MaterialProperty("SCINTILLATIONCOMPONENT2", 'eV', [x for x in zip(air_Energy, air_SCINT)], name=name),
        MaterialProperty("RINDEX", 'eV', [x for x in zip(
            air_Energy, air_RIND)], name=name),
        MaterialProperty("ABSLENGTH", 'eV', [x for x in zip(
            air_Energy, air_ABSL)], name=name, secondUnit='cm'),
        # MaterialProperty("SCINTILLATIONYIELD", '1/MeV', [[12000.]], name=name),
        # MaterialProperty("RESOLUTIONSCALE", '', [[1.0]], name=name),
        # MaterialProperty("SCINTILLATIONTIMECONSTANT1", 'ns', [[20.]], name=name),
        # MaterialProperty("SCINTILLATIONTIMECONSTANT2", 'ns', [[45.]], name=name),
        # MaterialProperty("SCINTILLATIONYIELD1", 'ns', [[1.0]], name=name),
        # MaterialProperty("SCINTILLATIONYIELD2", 'ns', [[0.0]], name=name),
    ]

    return med_AIR, mat_AIR, props


def LXe(name='LXe'):
    mat_LXe = r.TGeoMaterial("mat_"+name,
                             131.29,  # a , g/mole
                             54,     # z
                             3.020,  # g/cm3, rho (at 87 Kelvin)
                             2.87,   # cm, radiation length
                             85.7    # cm, int length
                             )
    med_LXe = r.TGeoMedium(name, 1, mat_LXe)

    # all_root_materials.append(mat_LXe)
    # all_root_materials.append(med_LXe)

    # Taken from LXeDetectorConstruction.cc in geant4 example extended/optical/LXe
    lxe_Energy = [7.0, 7.07, 7.14]
    lxe_SCINT = [0.1, 1.0, 0.1]
    lxe_RIND = [1.59, 1.57, 1.54]
    lxe_ABSL = [35., 35., 35.]

    props = [
        MaterialProperty("FASTCOMPONENT", 'eV', [
                         x for x in zip(lxe_Energy, lxe_SCINT)], name=name),
        MaterialProperty("SLOWCOMPONENT", 'eV', [
                         x for x in zip(lxe_Energy, lxe_SCINT)], name=name),
        MaterialProperty("RINDEX", 'eV', [x for x in zip(
            lxe_Energy, lxe_RIND)], name=name),
        MaterialProperty("ABSLENGTH", 'eV', [x for x in zip(
            lxe_Energy, lxe_ABSL)], name=name, secondUnit='cm'),
        MaterialProperty("SCINTILLATIONYIELD", '1/MeV', [[12000]], name=name),
        MaterialProperty("RESOLUTIONSCALE", '', [[1.0]], name=name),
        MaterialProperty("FASTTIMECONSTANT", 'ns', [[20.]], name=name),
        MaterialProperty("SLOWTIMECONSTANT", 'ns', [[45.]], name=name),
        MaterialProperty("YIELDRATIO", '', [[1.0]], name=name),
    ]

    return med_LXe, mat_LXe, props


def LSO(name):
    # See: https://www.advatech-uk.co.uk/lso_ce.html
    pCe = 0.05
    pLu = (1-pCe)*2/8.
    pSi = (1-pCe)*1/8.
    pO = (1-pCe)*5/8.

    mat_lso = r.TGeoMixture(name, 4, 7.4)                  # Define a mixture
    mat_lso.DefineElement(0, 140.12, 58, pCe)
    mat_lso.DefineElement(1, 15.999, 8, pO)
    mat_lso.DefineElement(2, 28.085, 14, pSi)
    mat_lso.DefineElement(3, 174.97, 71, pLu)
    med_lso = r.TGeoMedium("med_"+name, 1, mat_lso)

    return med_lso, mat_lso, None


def CsI(name='miami'):
    # see https://www.advatech-uk.co.uk/csi_na.html

    density = 4.51

    # Define a mixture
    mat_CsI = r.TGeoMixture(name, 2, density)
    mat_CsI.DefineElement(0, 132.91, 55, 0.5)
    mat_CsI.DefineElement(1, 126.90, 53, 0.5)
    med_CsI = r.TGeoMedium("med_"+name, 1, mat_CsI)

    return med_CsI, mat_CsI, None


def NaI(name):
    # see: https://www.advatech-uk.co.uk/nai_tl.html

    density = 3.67
    # Define a mixture
    mat_NaI = r.TGeoMixture(name, 2, density)
    mat_NaI.DefineElement(0, 22.990, 11, 0.5)
    mat_NaI.DefineElement(1, 126.90, 53, 0.5)
    med_NaI = r.TGeoMedium("med_"+name, 1, mat_NaI)

    return med_NaI, mat_NaI, None


def LYSO(name):
    # (Lu[1-x]Y[x])[2] SiO[5]:0.01Ce
    # See: https://www.sciencedirect.com/science/article/pii/S0030402618300445

    x = 0.05  # ranges from 0 - 0.2 in paper
    pCe = 0.01
    pO = (1-pCe)*5/8.
    pSi = (1-pCe)*1/8.
    pLu = (1-pCe)*2/8.*(1-x)
    pY = (1-pCe)*2/8.*(x)

    mat_lyso = r.TGeoMixture(name, 5, 7.4)                  # Define a mixture
    mat_lyso.DefineElement(0, 140.12, 58, pCe)
    mat_lyso.DefineElement(1, 15.999, 8, pO)
    mat_lyso.DefineElement(2, 28.085, 14, pSi)
    mat_lyso.DefineElement(3, 174.97, 71, pLu)
    mat_lyso.DefineElement(4, 88.906, 39, pY)
    med_lyso = r.TGeoMedium("med_"+name, 1, mat_lyso)

    with open("config_files/lyso_emission.json", 'r') as f:
        emission_spectra = json.load(f)
    lyso_Lambda, lyso_Emission = zip(*emission_spectra['emission'])
    lyso_Energy = 1240./np.array(lyso_Lambda)  # 1240 = hc in eV*nm
    lyso_Emission = np.array(lyso_Emission)/np.amax(lyso_Emission)
    lyso_RIND = np.full_like(lyso_Energy, emission_spectra['refractive_index'])
    lxe_ABSL = np.full_like(lyso_Energy, 20.)

    props = [
        MaterialProperty("FASTCOMPONENT", 'eV', [x for x in zip(
            lyso_Energy, lyso_Emission)], name=name),
        MaterialProperty("SLOWCOMPONENT", 'eV', [x for x in zip(
            lyso_Energy, lyso_Emission)], name=name),
        MaterialProperty("RINDEX", 'eV', [x for x in zip(
            lyso_Energy, lyso_RIND)], name=name),
        MaterialProperty("ABSLENGTH", 'eV', [x for x in zip(
            lyso_Energy, lxe_ABSL)], name=name, secondUnit='cm'),
        MaterialProperty("SCINTILLATIONYIELD", '1/MeV',
                         [[emission_spectra['light_yield']]], name=name),
        MaterialProperty("RESOLUTIONSCALE", '', [[1.0]], name=name),
        MaterialProperty("FASTTIMECONSTANT", 'ns', [
                         [emission_spectra['decay_constant']]], name=name),
        MaterialProperty("SLOWTIMECONSTANT", 'ns', [
                         [emission_spectra['decay_constant']]], name=name),
        MaterialProperty("YIELDRATIO", '', [[1.0]], name=name),
    ]

    return med_lyso, mat_lyso, props


def PbW04(name):
    density = 8.28  # g/cm3
    mat_PbWO4 = r.TGeoMixture(name, 3, density)
    mat_PbWO4.DefineElement(0, 207.2, 82, 1/6.)
    mat_PbWO4.DefineElement(1, 183.84, 74, 1/6.)
    mat_PbWO4.DefineElement(2, 15.999, 8, 4/6.)
    med_PbWO4 = r.TGeoMedium(name, 2, mat_PbWO4)
    return med_PbWO4, mat_PbWO4, None


def carbon():
    density = 1.45  # *g/cm3;
    table = r.gGeoManager.GetElementTable()

    mat_CF = r.TGeoMixture("CarbonFiber", 1, density)
    mat_CF.AddElement(table.GetElement(6), 1.0)

    med_CF = r.TGeoMedium("med_CarbonFiber", 2, mat_CF)
    # all_root_materials.append(mat_CF)
    # all_root_materials.append(med_CF)

    return med_CF, mat_CF, None


def carbonfiber():
    # Modified from: https://gemc.jlab.org/work/doxy/1.8/cpp__materials_8cc_source.html
    table = r.gGeoManager.GetElementTable()

    # mat_CF = r.TGeoMixture("CarbonFiber", 1, density)
    # mat_CF.AddElement( table.GetElement(6), 1.0 )

    Epoxy = r.TGeoMixture("Epoxy", 4,  1.16)  # *g/cm3
    Epoxy.AddElement(table.GetElement(1), 32)  # Hydrogen
    Epoxy.AddElement(table.GetElement(7),  2)  # Nitrogen
    Epoxy.AddElement(table.GetElement(8),  4)  # Oxygen
    Epoxy.AddElement(table.GetElement(6), 15)  # Carbon
    # Epoxy.RegisterYourself()

    # MMats["Epoxy"] = Epoxy;

    percentEpoxy = 0.2550
    mat_CF = r.TGeoMixture("CarbonFiber", 2, 1.750)
    mat_CF.AddElement(Epoxy, percentEpoxy)
    mat_CF.AddElement(table.GetElement(6), 1.0-percentEpoxy)

    med_CF = r.TGeoMedium("med_CarbonFiber", 2, mat_CF)
    # all_root_materials.append(Epoxy)
    # all_root_materials.append(mat_CF)
    # all_root_materials.append(med_CF)

    return med_CF, mat_CF, None


def Beryllium(name):
    density = 1.84  # *g/cm3;
    table = r.gGeoManager.GetElementTable()

    mat_CF = r.TGeoMixture(name, 1, density)
    mat_CF.AddElement(table.GetElement(4), 1.0)

    med_CF = r.TGeoMedium("med_"+name, 2, mat_CF)
    # all_root_materials.append(mat_CF)
    # all_root_materials.append(med_CF)

    return med_CF, mat_CF, None


def Tedlar(name):
    density = 1.84  # *g/cm3;
    table = r.gGeoManager.GetElementTable()

    mat_CF = r.TGeoMixture(name, 3, density)
    mat_CF.AddElement(table.GetElement(6), 0.33333)
    mat_CF.AddElement(table.GetElement(1), 0.5)
    mat_CF.AddElement(table.GetElement(9), 0.16667)

    med_CF = r.TGeoMedium("med_"+name, 2, mat_CF)
    # all_root_materials.append(mat_CF)
    # all_root_materials.append(med_CF)

    return med_CF, mat_CF, None


def kapton(name):
    # From: https://pdg.lbl.gov/2020/AtomicNuclearProperties/HTML/polyimide_film.html
    # name = 'Kapton'
    table = r.gGeoManager.GetElementTable()
    mat_Kapton = r.TGeoMixture(name, 4, 1.420)  # density
    mat_Kapton.AddElement(table.GetElement(6), 0.691133)  # Carbon
    mat_Kapton.AddElement(table.GetElement(7), 0.073270)  # Nitrogen
    mat_Kapton.AddElement(table.GetElement(8), 0.209235)  # Oxygen
    mat_Kapton.AddElement(table.GetElement(1), 0.026362)  # Hydrogen
    med_Kapton = r.TGeoMedium("med_"+name, 4, mat_Kapton)
    # all_root_materials.append(mat_Kapton)
    # all_root_materials.append(med_Kapton)
    return med_Kapton, mat_Kapton, None


def kapton_cable(name):
    # From: https://pdg.lbl.gov/2020/AtomicNuclearProperties/HTML/polyimide_film.html
    # name = 'Kapton'
    # will be kapton with a percentage copper added, to simulate cable
    kapton_percentage = 0.5
    table = r.gGeoManager.GetElementTable()
    kapton_density = 1.420
    copper_density = 8.960
    mat_Kapton = r.TGeoMixture(name, 5,
                               kapton_density*kapton_percentage + copper_density*(1-kapton_percentage))
    mat_Kapton.AddElement(table.GetElement(
        6), kapton_percentage*0.691133)  # Carbon
    mat_Kapton.AddElement(table.GetElement(
        7), kapton_percentage*0.073270)  # Nitrogen
    mat_Kapton.AddElement(table.GetElement(
        8), kapton_percentage*0.209235)  # Oxygen
    mat_Kapton.AddElement(table.GetElement(
        1), kapton_percentage*0.026362)  # Hydrogen
    mat_Kapton.AddElement(table.GetElement(29), 1.0 -
                          kapton_percentage)  # Hydrogen
    med_Kapton = r.TGeoMedium("med_"+name, 4, mat_Kapton)
    # all_root_materials.append(mat_Kapton)
    # all_root_materials.append(med_Kapton)
    return med_Kapton, mat_Kapton, None


def SuperDense(name):
    density = 1000.84  # *g/cm3;
    table = r.gGeoManager.GetElementTable()

    mat_CF = r.TGeoMixture(name, 1, density)
    mat_CF.AddElement(table.GetElement(82), 1.0)

    med_CF = r.TGeoMedium("med_"+name, 2, mat_CF)
    # all_root_materials.append(mat_CF)
    # all_root_materials.append(med_CF)

    return med_CF, mat_CF, None


def detector_materials(name, debug=False):
    r.gSystem.Load("libGeom")
    if(debug):
        print("Scintillator name:", name)

    valid_types = [
        'Bicron408',
        'Silicon',
        'Silicon2',
        'SiliconDead',
        'IBMS',
        'LXe',
        'CaloAbsorber',
        'Lead',
        'Aluminum',
        'CaloGhostPlane',
        'CarbonFiber',
        'Kapton',
        'KaptonTarget',
        'KaptonCable',
        'TargetGhostPlane',
        'CarbonFiber',
        'Be',
        'Beryllium',
        'SuperDense',
        'Copper',
        'air',
        'Tedlar',
        'LYSO',
        'LSO',
        'NaI',
        'CsI',
        'PbWO4',
    ]

    name = name.strip()
    if(name.strip() not in valid_types):
        raise ValueError(
            "ERROR:", name, 'is not an implemented type of detector material. Please choose from', valid_types)

    table = r.gGeoManager.GetElementTable()

    if name == 'Bicron408':
        mat_BC408 = r.TGeoMixture(name, 2, 1.032)
        mat_BC408.DefineElement(0, 1.00797, 1, 0.08475)  # ;  // H
        mat_BC408.DefineElement(1, 12.01115, 6, 0.91525)  # ; // C
        med_BC408 = r.TGeoMedium("med_"+name, 18, mat_BC408)
        # all_root_materials.append(mat_BC408)
        # all_root_materials.append(med_BC408)
        return med_BC408, mat_BC408, None
    elif(name == 'Silicon' or name == 'Silicon2' or name == 'SiliconDead'):
        mat_Si = r.TGeoMixture(name, 1, 2.33)
        mat_Si.AddElement(table.GetElement(14), 1.0)

        med_Si = r.TGeoMedium("med_"+name, 2, mat_Si)
        # all_root_materials.append(mat_Si)
        # all_root_materials.append(med_Si)
        return med_Si, mat_Si, None
    elif(name == 'IBMS'):
        # From: http://kuraraypsf.jp/psf/
        mat_IBMS = r.TGeoMixture(name, 2, 1.05)  # density should be 1.05
        mat_IBMS.AddElement(table.GetElement(6), 0.5)  # Carbon
        mat_IBMS.AddElement(table.GetElement(1), 0.5)  # Hydrogen
        med_IBMS = r.TGeoMedium("med_"+name, 2, mat_IBMS)
        # all_root_materials.append(mat_IBMS)
        # all_root_materials.append(med_IBMS)
        return med_IBMS, mat_IBMS, None
    elif name == 'LXe':
        return LXe()
    elif name == "CaloGhostPlane":
        return air(name)
    elif name == "TargetGhostPlane":
        return air(name)
    elif name == 'CarbonFiber':
        return carbonfiber()
    elif name == 'KaptonCable':
        return kapton_cable(name)
    elif name == 'Kapton' or 'Kapton' in name:
        return kapton(name)
    elif name == 'Tedlar':
        return Tedlar(name)
    elif name == 'Lead':
        mat_Pb = r.TGeoMixture(name, 1, 11.35)
        mat_Pb.AddElement(table.GetElement(82), 1.0)

        med_Pb = r.TGeoMedium("med_"+name, 2, mat_Pb)
        # all_root_materials.append(mat_Pb)
        # all_root_materials.append(med_Pb)
        return med_Pb, mat_Pb, None
    elif name == 'Aluminum':
        mat_Al = r.TGeoMixture(name, 1, 2.7)
        mat_Al.AddElement(table.GetElement(13), 1.0)

        med_Al = r.TGeoMedium("med_"+name, 2, mat_Al)
        # all_root_materials.append(mat_Al)
        # all_root_materials.append(med_Al)
        return med_Al, mat_Al, None
    elif name == 'CaloAbsorber':
        mat_Pb = r.TGeoMixture(name, 1, 0.00001)
        mat_Pb.AddElement(table.GetElement(1), 1.0)

        med_Pb = r.TGeoMedium("med_"+name, 2, mat_Pb)
        # all_root_materials.append(mat_Pb)
        # all_root_materials.append(med_Pb)
        return med_Pb, mat_Pb, None
    elif name == 'Be' or name == 'Beryllium':
        return Beryllium(name)
    elif name == 'SuperDense':
        return SuperDense(name)
    elif name == 'air':
        return air(name)
    elif name == 'Copper':
        mat_Pb = r.TGeoMixture(name, 1, 8.96)
        mat_Pb.AddElement(table.GetElement(29), 1.0)
        med_Pb = r.TGeoMedium("med_"+name, 2, mat_Pb)
        # all_root_materials.append(mat_Pb)
        # all_root_materials.append(med_Pb)
        return med_Pb, mat_Pb, None
    elif name == 'LYSO':
        return LYSO(name)
    elif name == 'LSO':
        return LSO(name)
    elif name == 'NaI':
        return NaI(name)
    elif name == 'CsI':
        return CsI(name)
    elif name == 'PbWO4':
        return PbW04(name)
    else:
        raise NotImplementedError()
