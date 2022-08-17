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


def MyEOSFitting(atoms, db, fh, metal, functional):
    
    v0 = atoms.get_volume()
    cell0 = atoms.get_cell()
    
    start = time.time()
    eos = calculate_eos(atoms,npoints=5, eps=0.04,trajectory=fh+".traj")
    v, e, B = eos.fit()  # find minimum
    end = time.time()
    elapsed = end - start
    
    # Do one more calculation at the minimum and write to database:
    s = (v/v0)**(1./3.)
    atoms.set_cell(cell0*s, scale_atoms=True)
    atoms.get_potential_energy()
    
    B_GPa = B / kJ * 1.0e24
    
    db.write(atoms, bm=B_GPa, scaler=s, metal=metal, xcfunctional=functional, totaltime=elapsed)
    eos.plot(filename=fh+".png")
    
    
    

metal = sys.argv[1]
functional = sys.argv[2]
atoms.calc = CP2KCalculator(
        400, 
        functional,
        kpoints=(4,4,4),
        dipole_correction=False,
        orbital_transform=False,
        smearing=True)

db = connect('bulk.db')

# Read the metal we want to test
atoms = read("../../Resources/cif/tmc/"+metal+"_conventional_standard.cif")
MyEOSFitting(atoms,db,metal+"_"+functional,metal, functional)
