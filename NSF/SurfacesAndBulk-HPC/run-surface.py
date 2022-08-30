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
import nglview as nv

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


from mycalculators import *
cwd = os.getcwd()

        
metal = sys.argv[1]
functional = sys.argv[2]

miller = sys.argv[3]

miller_a = int(miller[0])
miller_b = int(miller[1])
miller_c = int(miller[2])
#miller_indices = [(0,0,1), (1,0,1), (1,2,1), (1,2,0), (0,0,2)]
miller_indices = [(miller_a, miller_b, miller_c)]

db1 = connect('bulk.db')
db2 = connect('surface.db')

# Read the metal we want to test
atoms = read("../../Resources/cif/tmc/"+metal+"_conventional_standard.cif")
calc = CP2KCalculator( 400, functional, kpoints=(4,4,1), dipole_correction=False, orbital_transform=False, smearing=True, added_mos=30)

for idx in miller_indices:
    midx_str = str(idx[0]) + str(idx[1]) + str(idx[2])
    fh = metal + "_" + functional +"_" + midx_str

    start = time.time()
    
    row = db1.get(xcfunctional=functional, metal=metal)
    uc = row.toatoms()
    uc_formula = uc.get_chemical_formula()
    
    # Get the Bulk Energy and Formula from the database
    e_bulk = uc.get_potential_energy()
    uc_symbols = uc.get_chemical_symbols()
    # Create the Surface
    unit_surface = surface(uc, indices=idx, layers=2, vacuum=5)
    unit_surface.set_pbc(True)
    unit_surface = sort(unit_surface, unit_surface.positions[:,2])
    unit_surface.calc = calc
    
    c = FixAtoms(indices=[atom.index for atom in unit_surface if atom.position[2] < 0.5*np.sum(unit_surface.get_cell()[:,2])])
    unit_surface.set_constraint(c)
                
    #Figure out how many unitcells make up this surface
    unit_surface_symbols = unit_surface.get_chemical_symbols()
    n_uc = len(unit_surface_symbols) // len(uc_symbols)
                
    #Get the Surface Area 
    surface_area = np.linalg.norm(np.cross(unit_surface.get_cell()[0], unit_surface.get_cell()[1]))
                
    #Do one SCF cycle to get the unrelaxed structure energy
    e_surface_unrelaxed = unit_surface.get_potential_energy()
                
    #Minimize the surface to get the relaxed structure
    minimizer = FIRE(unit_surface, trajectory=fh+"_fire.traj", logfile=fh+"_fire.log")
    minimizer.run(fmax=0.03)
    e_surface_relaxed = unit_surface.get_potential_energy()
                
    #calculate surface enrgy
    e_cleavage = 0.5*(e_surface_unrelaxed - n_uc*e_bulk)
    e_relaxation = e_surface_relaxed - e_surface_unrelaxed
    surface_energy = (e_cleavage + e_relaxation)/surface_area
    end = time.time()
    elapsed = end - start

    db2.write(unit_surface, miller=str(idx), metal=metal, xcfunctional=functional, totaltime=elapsed, ce=e_cleavage, re=e_relaxation, se=surface_energy)
