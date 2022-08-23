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
from copy import deepcopy

# For building things
from ase import Atom, Atoms
from ase.io import read, write
from ase.io.trajectory import Trajectory
from ase.build import molecule, sort, add_adsorbate, surface
from ase.neighborlist import NeighborList, natural_cutoffs

# Unit Conversions and Fixing Atoms
from ase.units import Bohr,Rydberg,kJ,kB,fs,Hartree,mol,kcal
from ase.constraints import FixedPlane, FixedLine, FixAtoms

# ASE Calculators
#from plumed import Plumed
from ase.calculators.cp2k import CP2K
from ase.calculators.lj import LennardJones
#from ase.calculators.plumed import Plumed
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


zeo = read("MWW.cif")


calc = CP2KCalculator(ecut=400, functional="rVV10", scf=20, orbital_transform=True, voronoi=False)


# # EOS for Zeolite

zeo.calc = calc

v0 = zeo.get_volume()
cell0 = zeo.get_cell()

start = time.time()
eos = calculate_eos(zeo,npoints=5, eps=0.04,trajectory="zeo-eos.traj")
v, e, B = eos.fit()  # find minimum
end = time.time()
elapsed = end - start

# Do one more calculation at the minimum and write to database:
s = (v/v0)**(1./3.)
zeo.set_cell(cell0*s, scale_atoms=True)
zeo.get_potential_energy()
write("zeo-optimized.xyz", zeo)
B_GPa = B / kJ * 1.0e24

eos.plot(filename="zeo.png")


# ## Creating the Nanosheet

rOH = 1.0
rSiO = 1.61

angleSiOH = np.deg2rad(109.5)
angleOSiO = np.deg2rad(109.5)

supercell = surface(zeo, indices=(0,0,1), layers=1, vacuum=10.0, tol=1e-10, periodic=True)

# Find out how much of the bottom half of the zeolite we are taking out
cutoffz = 0.5*np.sum(supercell.get_cell()[:,2])
minz = np.min(supercell.get_positions()[:,2])

remove_vacuum = cutoffz - minz
a, b, c = supercell.get_cell()
del supercell[[atom.index for atom in supercell if atom.position[2] < cutoffz]]
supercell.set_cell([a, b, c - np.array([0., 0., remove_vacuum])])
supercell.center()

cuts = natural_cutoffs(supercell)
nl = NeighborList(cuts, self_interaction=False, bothways=True)
nl.update(supercell)

# Loop over all atoms
natoms = len(supercell)
for idx in range(0, natoms):

        # Get neighbors for atom i
        indices, offsets = nl.get_neighbors(idx)

        # If atom i is and Oxygen with only 1 neighbor, it needs H to satisfy its valency
        if supercell[idx].symbol == 'O' and len(indices) < 2:

                # Get corrdinates of Oxygen Atom
                terminal_O_xyz = supercell[idx].position

                # Get corrdinates of neighboring atoms
                for i, offset in zip(indices, offsets):
                    neighboring_Si_xyz = supercell.positions[i] + np.dot(offset, supercell.get_cell())


                # Generate vector
                position_vector_Si_to_O = terminal_O_xyz - neighboring_Si_xyz


                # Use spherical coordinates to make coordinates for H atom
                if neighboring_Si_xyz[2] >= 0:
                    theta = 2.0*np.random.rand()*np.pi
                    terminal_H_x = rOH * np.cos(theta) * np.sin(angleSiOH)
                    terminal_H_y = rOH * np.sin(theta) * np.sin(angleSiOH)
                    terminal_H_z = rOH * np.cos(angleSiOH)
                    terminal_H_x = neighboring_Si_xyz[0] + position_vector_Si_to_O[0] + terminal_H_x
                    terminal_H_y = neighboring_Si_xyz[1] + position_vector_Si_to_O[1] + terminal_H_y
                    terminal_H_z = neighboring_Si_xyz[2] + position_vector_Si_to_O[2] + terminal_H_z
                else:
                    theta = 2.0*np.random.rand()*np.pi
                    terminal_H_x = rOH * np.cos(theta) * np.sin(angleSiOH)
                    terminal_H_y = rOH * np.sin(theta) * np.sin(angleSiOH)
                    terminal_H_z = rOH * np.cos(angleSiOH)
                    terminal_H_x = neighboring_Si_xyz[0] + position_vector_Si_to_O[0] - terminal_H_x
                    terminal_H_y = neighboring_Si_xyz[1] + position_vector_Si_to_O[1] - terminal_H_y
                    terminal_H_z = neighboring_Si_xyz[2] + position_vector_Si_to_O[2] - terminal_H_z

                Hatom = Atom('H', [terminal_H_x, terminal_H_y, terminal_H_z])

                supercell = supercell + Hatom

        # If atom i is Si with 3 neighbors, it needs OH to satisfy its valency
        if supercell[idx].symbol == 'Si' and len(indices) < 4:
            # Get corrdinates of Si Atom
            terminal_Si_xyz = supercell[idx].position

            # Get corrdinates of neighboring atoms
            neighbor_O_xyzs = []
            for i, offset in zip(indices, offsets):
                neighbor_O_xyzs.append(supercell.positions[i] + np.dot(offset, supercell.get_cell()))

            # Find the center of the circle made by the three oxygen atoms
            A = neighbor_O_xyzs[0]
            B = neighbor_O_xyzs[1]
            C = neighbor_O_xyzs[2]
            u1 = B - A
            w1 = np.cross(C-A, u1)

            u = u1 / np.linalg.norm(u1)
            w = w1 / np.linalg.norm(w1)
            v = np.cross(w,u)

            b = [np.dot(B-A,u), 0.0]
            c = [np.dot(C-A,u),np.dot(C-A,v)]

            h = ((c[0] - 0.5*b[0])**2 + c[1]**2 - (0.5*b[0])**2) / (2.0 * c[1])
            center = A + 0.5*b[0]*u + h*v

            # Generate vector
            vector_center_to_Si = terminal_Si_xyz - center
            unit_vector = vector_center_to_Si / np.linalg.norm(vector_center_to_Si)
            position_O = terminal_Si_xyz + rSiO*unit_vector

            Oatom =  Atoms('O', [position_O])

            vector_Si_to_O = position_O - terminal_Si_xyz

            # Use spherical coordinates to make coordinates for H atom
            theta = 2.0*np.random.rand()*np.pi
            terminal_H_x = rOH * np.cos(theta) * np.sin(angleSiOH)
            terminal_H_y = rOH * np.sin(theta) * np.sin(angleSiOH)
            terminal_H_z = rOH * np.cos(angleSiOH)

            if terminal_Si_xyz[2] >= 0:
                terminal_H_x = terminal_Si_xyz[0] + vector_Si_to_O[0] - terminal_H_x
                terminal_H_y = terminal_Si_xyz[1] + vector_Si_to_O[1] - terminal_H_y
                terminal_H_z = terminal_Si_xyz[2] + vector_Si_to_O[2] - terminal_H_z

            else:
                terminal_H_x = terminal_Si_xyz[0] + vector_Si_to_O[0] + terminal_H_x
                terminal_H_y = terminal_Si_xyz[1] + vector_Si_to_O[1] + terminal_H_y
                terminal_H_z = terminal_Si_xyz[2] + vector_Si_to_O[2] + terminal_H_z

            Hatom = Atom('H', [terminal_H_x, terminal_H_y, terminal_H_z])

            OH = Oatom + Hatom
            supercell = supercell + OH
            print()


supercell= sort(supercell, supercell.positions[:,2])

supercell.calc = calc

# Freeze every atom that is not a hydrogen, then optimize
c = FixAtoms(indices=[atom.index for atom in supercell if atom.symbol != "H"])
supercell.set_constraint(c)

minimizer = FIRE(supercell, trajectory="fire.traj", logfile="fire.log")
minimizer.run(fmax=0.03)

# After optimizing the Hydrogens, freeze the bottom half and do one more optimization
supercell.set_constraint()
c = FixAtoms(indices=[atom.index for atom in supercell if atom.position[2] < 0.5*supercell.get_cell()[2][2]])
supercell.set_constraint(c)
minimizer.run(fmax=0.03)


