from ase.calculators.cp2k import CP2K
from ase.units import Rydberg

def CP2KCalculator(ecut, functional="LDA", kpoints=None, dipole_correction=False,
                   orbital_transform=False,smearing=False,
                   voronoi=False, cube=False, bqb=False,v_hartree=False, added_mos=None,scf=50):
    """Creates a CP2K calculator object with Settings for Production Runs"""

    # By Default, assume we want to have the walltime as just shy of 48 hours
    inp = '''
&GLOBAL
WALLTIME 47:59:00
&END GLOBAL
&FORCE_EVAL
&DFT
'''

    ### DFT SECTION
    if dipole_correction:
        inp += '''
SURFACE_DIPOLE_CORRECTION .TRUE.
SURF_DIP_DIR Z
'''
    if kpoints is not None:

        s = "SCHEME MONKHORST-PACK " + str(kpoints[0]) + " " + str(kpoints[1]) + " " + str(kpoints[2]) + "\n"
        inp += '''
&KPOINTS
'''
        inp += s
        inp += '''
&END KPOINTS
'''

    ### SCF SECTION
    inp += '''
&SCF
&OUTER_SCF .TRUE.
MAX_SCF 50
&END OUTER_SCF
'''

    if orbital_transform:
        inp +='''
&OT .TRUE.
MINIMIZER DIIS
ALGORITHM IRAC
PRECONDITIONER FULL_SINGLE_INVERSE
&END OT
'''

    if smearing:
        if added_mos is not None:
            mos = added_mos
        else:
            mos = 10
        inp += "ADDED_MOS " + str(mos) + "\n"
        inp +='''
&SMEAR ON
METHOD FERMI_DIRAC
ELECTRONIC_TEMPERATURE [K] 300
&END SMEAR
&MIXING .TRUE.
METHOD BROYDEN_MIXING
&END MIXING
'''


    ###CLOSE SCF SECTION
    inp += '''
&END SCF
    '''

    ### XC Section
    inp += '''
&XC
&XC_GRID
XC_DERIV NN10_SMOOTH
XC_SMOOTH_RHO NN10
&END
'''

    if functional == "PBE+D3":
        functional="PBE"
        inp += '''
''''''
&VDW_POTENTIAL
POTENTIAL_TYPE PAIR_POTENTIAL
&PAIR_POTENTIAL
R_CUTOFF 15.0
TYPE DFTD3
VERBOSE_OUTPUT
CALCULATE_C9_TERM .FALSE.
REFERENCE_FUNCTIONAL PBE
PARAMETER_FILE_NAME dftd3.dat
&END PAIR_POTENTIAL
&END VDW_POTENTIAL
'''

    if functional == "optB88-vdw":
        functional = None
        inp += '''
&XC_FUNCTIONAL

&GGA_X_OPTB88_VDW
&END GGA_X_OPTB88_VDW
&PW92
&END PW92
&END XC_FUNCTIONAL
&vdW_POTENTIAL
DISPERSION_FUNCTIONAL NON_LOCAL
&NON_LOCAL
TYPE DRSLL
VERBOSE_OUTPUT
KERNEL_FILE_NAME vdW_kernel_table.dat
&END NON_LOCAL
&END vdW_POTENTIAL
'''

    if functional == "optB86B-vdw":
        functional = None
        inp += '''
&XC_FUNCTIONAL
&GGA_X_OPTB86B_VDW
&END
&PW92
&END PW92
&END XC_FUNCTIONAL
&vdW_POTENTIAL
DISPERSION_FUNCTIONAL NON_LOCAL
&NON_LOCAL
TYPE DRSLL
VERBOSE_OUTPUT
KERNEL_FILE_NAME vdW_kernel_table.dat
&END NON_LOCAL
&END vdW_POTENTIAL
'''
    if functional == "optPBE-vdw":
        functional = None
        inp += '''
&XC_FUNCTIONAL
&GGA_X_OPTPBE_VDW
&END
&PW92
&END PW92
&END XC_FUNCTIONAL
&vdW_POTENTIAL
DISPERSION_FUNCTIONAL NON_LOCAL
&NON_LOCAL
TYPE DRSLL
VERBOSE_OUTPUT
KERNEL_FILE_NAME vdW_kernel_table.dat
&END NON_LOCAL
&END vdW_POTENTIAL
'''

    if functional == "rVV10":
        functional = None
        inp += '''
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
'''


    ### CLOSE OFF XC
    inp += '''
&END XC
'''
        ### DFT Print SECTION
    if voronoi or cube or bqb or v_hartree:
        inp += '''
&PRINT
'''
        if cube:
            inp += '''
&E_DENSITY_CUBE
&END E_DENSITY_CUBE
'''
        if bqb:
            inp += '''
&E_DENSITY_BQB
&END E_DENSITY_BQB
'''
        if v_hartree:
            inp += '''
&V_HARTREE_CUBE
&END V_HARTREE_CUBE
'''
        if voronoi:
            inp += '''
&VORONOI ON
OUTPUT_TEXT .TRUE.
&END VORONOI
'''
        inp += '''
&END PRINT
'''

    #### CLOSE EVERYTHING
    inp +='''
&END DFT
&END FORCE_EVAL
'''

    calc = CP2K(
        auto_write=False,
        basis_set="DZVP-MOLOPT-SR-GTH",
        basis_set_file="BASIS_MOLOPT",
        charge=0,
        uks = True,
        cutoff = ecut*Rydberg,
        debug = False,
        force_eval_method = "Quickstep",
        xc = functional,
        inp = inp,
        max_scf = scf,
        poisson_solver ="auto",
        potential_file = "POTENTIAL",
        pseudo_potential = "GTH-PBE",
        stress_tensor = True,
        print_level = "LOW",

    )

    return calc
