###############################################################################
# import the dipole hyperfine coupling package
import pyDipoleHF
# import further packages
import numpy as np
import time


###############################################################################
###############################################################################
# input parameters, do not modify code outside this section if you dont #######
# understand what you are doing ###############################################
###############################################################################
###############################################################################
# instantiate the classes and define the relevant attributes therein. See the 
# following file for details: dipole_hyperfine_lattice_sum/ipynb_examples/
# Mn226/Mn226_2021-08-30.ipynb

input_pars = pyDipoleHF.input_parameters()
# path to cif file
input_pars.cif_file_name = 'Mn2P2S6_Ouvrard_1985_Structuraldeterminationsome_EntryWithCollCode61391.cif'

# create a supercell if one is required to define the magnetic structure
# the supercell_dimensions parameter should be a list of form [Na, Nb, Nc],
# where Na, Nb, and Nc are integers that represent the number of unit cells
# in the crystallographic a, b, and c directions, respectivly, contained in
# the supercell. For example, [1, 1, 2] will create a supercell that is 
# 1 x 1 x 2 unit cells in size. Note: creating a supercell modifies the
# lattice as well, therefore one must specify the magnetic_site_indices,
# probe indices, etc. in the supercell. see the following for deteails
# dipole_hyperfine_lattice_sum/ipynb_examples/LiCuVO4/
# LiCuVO4_dipole_hyperfine_calc_2021-08-30.ipynb
#supercell_dimensions = [1, 1, 1]

# define the magnetic sites by index in the cif file (list of ints)
# Note: that site indices in pymatgen start at zero.
# Further Note: if creating a supercell, one must include all magetic 
# site indices in the supercell. therefore the user must first know the 
# site indices.
input_pars.magnetic_site_indices = [0, 1, 2, 3]

# static magnetic order input parameters
# magnetic_structure_input_method has two possible (string) values:
# 'set_moments_by_site_index' or 'set_moments_by_function'
# If the user does not care about static magnetic order or internal static
# dipolar hyperfine fields, then leave it as the default value None
# The first input method is direct assignment of magnetic moment vectors
# by magnetic site index. This method is easier, but does not allow the 
# user to define complicated and incommensurate magnetic orders. The 
# second method is more general and can be used for any case, but requires 
# the user to first calculate the magnetic propagation vector.
input_pars.magnetic_structure_input_method = 'set_moments_by_site_index'

# first_cell_magnetic_moment_vectors should be set to be a numpy array of
# shape = (len(magnetic_site_indices), 3), where the 0th axis shape is the 
# number of magnetic sites, and the 1st axis defines the magnetic moment vectors in 
# cartesian or lattice basis. the order of the magnetic moment vectors should
# be the same as the magnetic_site_indices.
# Example:
# define antiferromagnetic order on the  magnetic sites with indices 0 and 1,
# so above magnetic_site_indices = [0, 1], and there is a spin on site 0 with a 1.2 Bohr
# magneton moment pointing in the z direction and a -1.2 Bohr magneton in the
# z direction, then first_cell_magnetic_moment_vectors = np.array([[0, 0, 1.2], [0, 0, -1.2]])
input_pars.first_cell_magnetic_moment_vectors = np.array([[0, 0,  4.1], 
                                                             [0, 0, -4.1],
                                                             [0, 0,  4.1], 
                                                             [0, 0, -4.1]])

# bool to identify the basis in which the input magnetic moments are defined. if 
# set to true, then the magnetic moment vectors, given in the lattice basis 
# will be transformed to cartesian vectors for the calculations, but often neutron 
# scattering reports the moments in the crystal lattice basis, and therefore the 
# option to supply input moments in the crystal lattice basis will save time.
#moment_vectors_in_lattice_basis = False

# defines the magnetic propagation (wave) vector(s) (numpy array with the same shape 
# as first_cell_magnetic_moment_vectors, fractional coordinates)
# example: if we have two magnetic sites, with the first being ferromagnetic 
# in c and afm in a and b, and the second being afm in all directions then 
# k_frac = np.array([[0.5, 0.5, 0], [0.5, 0.5, 0.5]])
# Note: users needs to figure this out themselves, should be in fractional coords
# of the chosen supercell
#k_frac = None
# magnetic moment vector in a certain direction (in cartesian or lattice basis
# depending in the same way as first_cell_magnetic_moment_vectors on the bool
# moment_vectors_in_lattice_basis. this allows for general description of 
# spin/spiral incommensurate order, numpy array with the same shape as 
# first_cell_magnetic_moment_vectors
#I = None

# choose the radius around which we will include spins in angstroms with 
# respect to the probe position (float)
input_pars.sphere_radius = 150

# define the probe site number (nuclear spin site) (int)
input_pars.probe_site_index = 7

# for performing large calculations, we want to attempt to keep enough 
# memory. To do this we'll trash some arrays that we don't need
# running on 32bit python would cause memory allocation errors, so to 
# make the program more portable, we can allow for that memory to be 
# freed by the garbage collector by setting (no longer needed arrays 
# to None)
input_pars.delete_plotting_arrays = False


# a class for returning results of the calculation
# all results will be in cartesian coordinates!!!
results = pyDipoleHF.results()
#contains the following attributes:
# results used for visualization:
# the structure object created from the cif file
#structure
# a matrix representing the the lattice unit cell as row vectors in cartesian coordinates, angstroms
#lattice
# cartesian positions of the electron magnetic moments, angstroms
#positions 
# the electron magnertic moment vectors in the cartesian basis in mu_B
#magnetic_moment_vectors
# main results
# hyperfine coupling tensor in T/mu_B in cartesian basis
#A
# dipolar hyperfine field in Tesla cartesian basis
#B
# magnetization in mu_B/atom
#M


# plot or not
run_visualization = True
# plot inclusion minimum/maximum in angtroms, all linked, but one can modify the values individually below
minmax=10


###############################################################################
###############################################################################
# end input section; do not modify below here unless you understand what ######
# you are doing! ##############################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################

#the output of this code should be:

##calc took 4.618 s to run for 1.362560e+05 sites
##structure
##Full Formula (Mn4 P4 S12)
##Reduced Formula: MnPS3
##abc   :   6.077000  10.524000   6.796000
##angles:  90.000000 107.350000  90.000000
##Sites (20)
##  #  SP         a        b       c
##---  ----  ------  -------  ------
##  0  Mn2+  0       0.33258  0
##  1  Mn2+  0       0.66742  0
##  2  Mn2+  0.5     0.83258  0
##  3  Mn2+  0.5     0.16742  0
##  4  P4+   0.9444  0        0.8314
##  5  P4+   0.0556  0        0.1686
##  6  P4+   0.4444  0.5      0.8314
##  7  P4+   0.5556  0.5      0.1686
##  8  S2-   0.2407  0        0.7503
##  9  S2-   0.7593  0        0.2497
## 10  S2-   0.7407  0.5      0.7503
## 11  S2-   0.2593  0.5      0.2497
## 12  S2-   0.7562  0.1612   0.7484
## 13  S2-   0.7562  0.8388   0.7484
## 14  S2-   0.2438  0.8388   0.2516
## 15  S2-   0.2438  0.1612   0.2516
## 16  S2-   0.2562  0.6612   0.7484
## 17  S2-   0.2562  0.3388   0.7484
## 18  S2-   0.7438  0.3388   0.2516
## 19  S2-   0.7438  0.6612   0.2516
##probe position (Angstroms)
##[[3.03469286]
## [5.262     ]
## [1.09367251]]
##M (mu_B/atom) =
##[[ 0.00000000e+00]
## [ 0.00000000e+00]
## [-9.90804954e-19]]
##A (T/mu_B) =
##[[ 3.21911763e-02 -4.57533317e-16  1.24904516e-03]
## [-4.57533317e-16  3.48326092e-02  5.08273978e-16]
## [ 1.24904516e-03  5.08273978e-16 -6.70237855e-02]]
##B (T) =
##[[ 4.92661467e-16]
## [-8.09533488e-03]
## [-1.38777878e-17]]



###############################################################################
# time and run the calculation
start_time = time.time()

pyDipoleHF.calc_dip_A(input_pars, results)

end_time = time.time()

if not input_pars.delete_plotting_arrays:
    N_sites = results.positions.shape[0]
    print('calc took', 
          round(end_time - start_time, 3), 
          's to run for', 
          '{:e}'.format(N_sites), 
          'sites')
else:
    print('calc took', 
          round(end_time - start_time, 3), 
          's to run')



###############################################################################
# print the results (all in Cartesian coordinates)

print('structure')
print(results.structure)
print('probe position (Angstroms)')
print(results.probe_position) # note, cartesian origin is at unit cell origin
# spin-lattice magnetization
print('M (mu_B/atom) =')
print(results.M)
# dipole hyperfine coupling tensor
print('A (T/mu_B) =')
print(results.A)
# internal hyperfine field at the probe site
print('B (T) =')
print(results.B)



###############################################################################
# display a visualization (note, should only do this for less than approx
# 50 atoms or it is not helpful

# define the point-inclusion limits as a tuple
x_min_rel = -minmax
x_max_rel = minmax
y_min_rel = -minmax
y_max_rel = minmax
z_min_rel = -minmax
z_max_rel = minmax
pos_limits = (x_min_rel, x_max_rel, y_min_rel, y_max_rel, z_min_rel, z_max_rel)

# set the axes limits as a tuple
axes_min = -minmax
axes_max = minmax
ax_range = (axes_min, axes_max, axes_min, axes_max, axes_min, axes_max)

# run visualization if desired
if run_visualization and not input_pars.delete_plotting_arrays:
    pyDipoleHF.visualization(input_parameters=input_pars,
                             results=results,
                             relative_position_limits=pos_limits,
                             axes_range=None,
                             elev=0.,
                             azim=0.,
                             plot_sphere=False,
                             plot_lattice_vectors=True)
