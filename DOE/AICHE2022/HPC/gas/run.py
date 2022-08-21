
# General
import os
import sys
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import time
from copy import deepcopy

# For building things
from ase import Atom, Atoms
from ase.io import read, write
from ase.io.trajectory import Trajectory
from ase.build import molecule, sort, add_adsorbate, surface
from ase.visualize import view
from ase.neighborlist import NeighborList, natural_cutoffs

# Unit Conversions and Fixing Atoms
from ase.units import Bohr,Rydberg,kJ,kB,fs,Hartree,mol,kcal
from ase.constraints import FixedPlane, FixedLine, FixAtoms

# ASE Calculators
from plumed import Plumed
from ase.calculators.cp2k import CP2K
from ase.calculators.lj import LennardJones
from ase.calculators.plumed import Plumed
from ase.calculators.idealgas import IdealGas

# Geometry Optimizations and Normal Mode Analysis
from ase.optimize import LBFGS, FIRE
from ase.vibrations import Vibrations
from ase.thermochemistry import IdealGasThermo

# EOS fitting for Unit Cells
from ase.eos import EquationOfState, calculate_eos

# Molecular Dynamics
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md.verlet import VelocityVerlet
from ase.md.langevin import Langevin
from ase.md.npt import NPT

# Import my calculators
from mycalculators import *

def main():

    nametag = sys.argv[1]
    ligninfile=sys.argv[2]
    zeolitefile=sys.argv[3]

    acidity = sys.argv[4]

    bo4 = read(ligninfile)
    traj = Trajectory(zeolitefile)
    unitsurface = traj[-1]

    calc = CP2KCalculator(ecut=400, functional="rVV10", scf=20, orbital_transform=True, voronoi=True)

    # Atom Indices of the Si atoms we want to substitute for Al, as well as the O atoms for the BA Site
    T1 = 209
    T4 = 185
    T1_O = 201
    T4_O = 191

    if acidity == "pureSi":
        supercell = unitsurface
    elif acidity == "Al1":
        supercell = create_bronsted_acid_site(unitsurface, T1, T1_O, rOH=1.0)
    elif acidity == "Al4":
        supercell = create_bronsted_acid_site(unitsurface, T4, T4_O, rOH=1.0)
    elif acidity == "Al4-Al1":
        supercell = create_bronsted_acid_site(unitsurface, T4, T4_O, rOH=1.0)

    bo4_supercell = create_probe_on_zeo(supercell, bo4, nrepeats=(2,2,1))

    bo4_supercell.calc = calc
    set_tags(bo4_supercell, len(supercell.repeat((2,2,1))), len(bo4))

    # Freeze two atoms in the middle of the zeolite
    idx_1 = (len(supercell.repeat((2,2,1))) - len(bo4)) // 2
    idx_2 = (len(supercell.repeat((2,2,1))) - len(bo4)) // 2 + 1

    c = FixAtoms(indices=[idx_1, idx_2])
    bo4_supercell.set_constraint(c)

    runMD(bo4_supercell, name=nametag, temp=450.0, runlength=2000)


def set_tags(atoms, slablength, probelength):

    for i in range(0, slablength):
        atoms[i].tag = 0

    for i in range(slablength, probelength+slablength):
        atoms[i].tag = 1

def runMD(atoms, name, temp=298.0, runlength=40000):

    MaxwellBoltzmannDistribution(atoms, temperature_K=300)

    md = Langevin(atoms, 0.5*fs, temperature_K = temp, friction=0.01,
              trajectory=name+"_langevin.traj", logfile=name+'_langevin.log')

    md.run(runlength)


def create_probe_on_zeo(zeolite, probe, height=1.5, nrepeats=(1,1,1)):

    system = zeolite.repeat(nrepeats)
    system = sort(system, system.positions[:,2])

    box = system.get_cell()
    x = 0.5*(box[0][0] + box[1][0] + box[2][0])
    y = 0.5*(box[0][1] + box[1][1] + box[2][1])

    # Get the center of mass of the probe
    r_com = probe.get_center_of_mass()

    # Find out how far the COM is from index 0
    delta = probe.positions[0] - r_com

    add_adsorbate(system, probe, height, position=(x+delta[0], y+delta[1]), offset=None, mol_index=0)

    return system


def create_bronsted_acid_site(atoms, Si_index, O_index, rOH=1.0):

    # Make the Symbol for the Si_index Aluminum

    tmp = deepcopy(atoms)
    original_elements = tmp.get_chemical_symbols()

    original_elements[Si_index] = "Al"
    tmp.set_chemical_symbols(original_elements)

    # Get the coordinates of the Al atom and the Si atom next to the oxygen
    #r_Al = atoms.positions[Si_index]
    cuts = natural_cutoffs(tmp)
    nl = NeighborList(cuts, self_interaction=False, bothways=True)
    nl.update(tmp)

    # Get neighbors for atom i
    indices, offsets = nl.get_neighbors(O_index)

    # Get corrdinates of neighboring atoms
    neighbor_O_xyzs = []
    for i, offset in zip(indices, offsets):
        neighbor_O_xyzs.append(tmp.positions[i] + np.dot(offset, tmp.get_cell()))

    midpoint = 0.5*(neighbor_O_xyzs[0]+neighbor_O_xyzs[1])
    vector_midpoint_to_oxygen = tmp.positions[O_index] - midpoint
    unit_vector = vector_midpoint_to_oxygen / np.linalg.norm(vector_midpoint_to_oxygen)
    H_position = midpoint + vector_midpoint_to_oxygen + rOH*unit_vector

    Hatom = Atom('H', H_position)
    return tmp + Hatom



main()
