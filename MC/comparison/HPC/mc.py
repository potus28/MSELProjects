
from ase import Atoms
from ase.calculators.cp2k import CP2K
from ase.mc import GCMC
from ase.units import kJ, mol, kB
from ase.io.trajectory import Trajectory
import numpy as np


def myCP2K():
    inp = """
&FORCE_EVAL
&DFT

&SCF
&OUTER_SCF .TRUE.
MAX_SCF 20
&END OUTER_SCF
&OT .TRUE.
MINIMIZER DIIS
PRECONDITIONER FULL_SINGLE_INVERSE
&END OT
&END SCF

&XC

&XC_GRID
XC_DERIV NN10_SMOOTH
XC_SMOOTH_RHO NN10
&END XC_GRID

&XC_FUNCTIONAL
&GGA_X_RPW86
&END GGA_X_RPW86
&GGA_C_PBE
&END GGA_C_PBE
&END XC_FUNCTIONAL
&vdW_POTENTIAL
DISPERSION_FUNCTIONAL NON_LOCAL
&NON_LOCAL
TYPE RVV10
VERBOSE_OUTPUT
KERNEL_FILE_NAME rVV10_kernel_table.dat
&END NON_LOCAL
&END vdW_POTENTIAL

&END
&END
&END
    """
    return CP2K(xc=None, inp=inp)


species0 = Atoms("Ar", [[0, 0, 0]])

vol = 7.5
species0.set_cell((vol, vol, vol))
species0.center()

species_list = [species0]

#atoms = species0.repeat((4, 4, 4))
#atoms.set_pbc(True)

traj = Trajectory("mc-cp2k.traj")
atoms = traj[-1]
atoms.wrap()
ncycles = len(traj)



tag = -1
tags = []
for atom in atoms:
    if atom.index % 1 == 0:
        tag += 1
        tags.append(tag)
    atom.tag = tag

species_tags = []
species_tags.append(tags)

atoms.calc = myCP2K()


mc = GCMC(
    atoms=atoms,
    species=species_list,
    species_tags=species_tags,
    temperature_K=100,
    chemical_potential=np.array([-9.5]) * kJ / mol,
    calc_function=myCP2K,
    thermal_interval=1000,
    logfile="mc-cp2k.log",
    trajectory="mc-cp2k.traj",
    loginterval=64,
    append_trajectory = True
    )

mc.run((1000 - ncycles) * 64)

