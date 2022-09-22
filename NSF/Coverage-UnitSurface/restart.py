#!/usr/bin/env python
# coding: utf-8


# General
import os
import sys
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import time

# For building things
from ase import Atom, Atoms
from ase.io import read, write
from ase.io.trajectory import Trajectory
from ase.build import molecule, surface, add_adsorbate, add_vacuum, sort
from ase.visualize import view
from ase.db import connect
from ase.geometry import get_layers

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

cwd = os.getcwd()

from mycalculators import *


def main():

    rootdir = "/work/hpc/users/wnw36/THESIS/MSELProjects/"
    database = rootdir + "NSF/SurfacesAndBulk-HPC/bulk.db"


    #metals = ["Alpha-Mo2C_mp-1552","Beta-Mo2C_mp-1221498"]
    metals = ["Beta-Mo2C_mp-1221498"]
    functional = "rVV10"
    unitcelldb = database

    calc = CP2KCalculator(
            ecut=400,
            functional=functional,
            scf=30,
            orbital_transform=True,
            voronoi=False
            )
    for metal in metals:
        if metal == "Alpha-Mo2C_mp-1552":
            miller=(1,2,0)
            repeats=(2,3,1)
            repeats=(1,1,1)
            layers = 3

        elif metal == "Beta-Mo2C_mp-1221498":
            miller=(1,0,1)
            repeats=(3,3,1)
            #repeats=(1,1,1)
            layers = 3

        db = connect(unitcelldb)
        row = db.get(xcfunctional=functional, metal=metal)
        uc = row.toatoms()
        slab = surface(lattice=uc, indices=miller,vacuum=10.0, layers=layers, periodic=True)

        slab = slab.repeat(repeats)
        slab = sort(slab, slab.positions[:,2])

        maxatom = slab[-1]

        system = slab + Atoms("H", [maxatom.position + np.array([0., 0., 1.5])])
        c = FixAtoms(indices=[atom.index for atom in system if atom.position[2] < 0.5*system.get_cell()[2][2]])
        system.set_constraint(c)

        # Use the Duterium mass for Hydrogen
        masses = system.get_masses()
        for i in range(0, len(masses)):
            if masses[i] == 1.008:
                masses[i] *= 2
        system.set_masses(masses)

        system.calc = calc
        md = Langevin(
            system,
            timestep = 1.0 * fs,
            temperature_K=450.0,
            friction=0.01,
            logfile = metal+"_langevin.log",
            trajectory = metal+"_langevin.traj"
        )

        md.run(1000)

        view(system)

def set_tags(atoms, slablength, probelength):

    for i in range(0, slablength):
        atoms[i].tag = 0

    for i in range(slablength, probelength+slablength):
        atoms[i].tag = 1


def create_probe_on_surface(zeolite, probe, height=1.5, nrepeats=(1,1,1)):

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

main()


