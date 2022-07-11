###############################################################################
# import the EFG point charge lattice sum package
import pyEFGPointCharge
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
# following file for details: EFG_point_charge_lattice_sum/ipynb_examples/
# aAlO2/aAl2O3_EFG_point_charge_calc_2021-10-18.ipynb

input_pars = pyEFGPointCharge.input_parameters()
# path to cif file
input_pars.cif_file_name = 'Ishizawa_1980_aAl2O3_300K_EntryWithCollCode10425.cif'

# create a supercell if one is required to define the magnetic structure
# the supercell_dimensions parameter should be a list of form [Na, Nb, Nc],
# where Na, Nb, and Nc are integers that represent the number of unit cells
# in the crystallographic a, b, and c directions, respectivly, contained in
# the supercell. For example, [1, 1, 2] will create a supercell that is 
# 1 x 1 x 2 unit cells in size. Note: creating a supercell modifies the
# lattice as well, therefore one must specify the magnetic_site_indices,
# probe indices, etc. in the supercell. see the following for details
# dipole_hyperfine_lattice_sum/ipynb_examples/LiCuVO4/
# LiCuVO4_dipole_hyperfine_calc_2021-08-30.ipynb
#input_pars.supercell_dimensions = [1, 1, 1]

# choose the radius around which we will include spins in angstroms with 
# respect to the probe position (float)
input_pars.sphere_radius = 3

# define the probe site number (nuclear spin site) (int)
input_pars.probe_site_index = 12

# for performing large calculations, we want to attempt to keep enough 
# memory. To do this we'll trash some arrays that we don't need
# running on 32bit python would cause memory allocation errors, so to 
# make the program more portable, we can allow for that memory to be 
# freed by the garbage collector by setting (no longer needed arrays 
# to None)
input_pars.delete_plotting_arrays = False

# probe nucleus identifying dictionary key (from isotope_data_dict) (string)
# e.g. '11B'
#input_pars.probe_nucleus = "17O"

# or if you want to set a custom nuclear quadrupole moment, do that here and
# the default for the probe nucleus will be overridden (barns)
input_pars.manual_nuclear_quad_moment = -0.0265

# in that case one also needs to define the nuclear spin I
input_pars.manual_nuclear_spin = 5/2

# define the charge sites by index in the cif file (list of ints)
# Note: that site indices in pymatgen start at zero.
# Further Note: if creating a supercell, one must include all magetic 
# site indices in the supercell. therefore the user must first know the 
# site indices.
input_pars.charge_site_indices = [i for i in range(30)]

# charges should be set to be a numpy array of shape = (len(charge_site_indices),), 
# the order of the charges should be the same as the charge_site_indices. units 
# are in multiples of the elementary charge e.
# Example:
# define antiferromagnetic order on the  magnetic sites with indices 0 and 1,
# so above charge_site_indices = [0, 1], and there is a charge of -2 on site 0
# and a +1 charge on site 1, then charges = np.array([-2.0, 1.0])
input_pars.first_cell_charges = [3, 3, 3, 3, 3,
                                 3, 3, 3, 3, 3,
                                 3, 3,-2,-2,-2,
                                -2,-2,-2,-2,-2,
                                -2,-2,-2,-2,-2,
                                -2,-2,-2,-2,-2]

# Sternheimer antishielding factor (\gamma_\inf)
#input_pars.gamma_sternheimer = -2.2

results = pyEFGPointCharge.results()
# a class for returning results of the calculation
# all results are in cartesian coordinates!!!
#class results:
    ## the structure object created from the cif file
    #structure = None
    ## a matrix representing the the lattice unit cell as row vectors in 
    ## cartesian coordinates (angstroms)
    #lattice = None
    ## Charge positions (angstroms)
    #charge_positions = None
    ## charges, in units of the elementary charge e
    #charges = None
    ## raw EFG tensor, ie not diagonalized (V m^-2)
    #V_ab_raw = None
    ## raw eigenvalues of the EFG tensor, ordered such that 
    ## |V_xx| <= |V_yy| <= |V_zz|; where x, y, and z no longer correspond to the 
    ## cartesian axes, but indicate the directions of the principal axes of the
    ## EFG tensor (V m^-2)
    #V_aa = None
    ## eigenvectors of the EFG tensor  (principal axes directions)
    ## (normalized, cartesian coordinates, a matrix of column vectors, 
    ## order is the same as the V_aa)
    #principal_axes = None
    ## eta = (V_xx - V_yy)/V_zz
    #eta = None
    ## nu_z = 3eQV_zz/(2I(2I-1)h)
    #nu_z_MHz = None
    ## nu_Q_MHz = nu_z_MHz(1 - eta**2/3)**0.5
    #nu_Q_MHz = None

###############################################################################


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
##calc took 0.439 s to run for 1.600000e+01 sites
##structure =
##Full Formula (Al12 O18)
##Reduced Formula: Al2O3
##abc   :   4.754000   4.754000  12.990000
##angles:  90.000000  90.000000 120.000000
##Sites (30)
##  #  SP           a         b         c    charge
##---  ----  --------  --------  --------  --------
##  0  Al3+  0         0         0.14772          3
##  1  Al3+  0         0         0.64772          3
##  2  Al3+  0         0         0.85228          3
##  3  Al3+  0         0         0.35228          3
##  4  Al3+  0.666667  0.333333  0.481053         3
##  5  Al3+  0.666667  0.333333  0.981053         3
##  6  Al3+  0.666667  0.333333  0.185613         3
##  7  Al3+  0.666667  0.333333  0.685613         3
##  8  Al3+  0.333333  0.666667  0.814387         3
##  9  Al3+  0.333333  0.666667  0.314387         3
## 10  Al3+  0.333333  0.666667  0.518947         3
## 11  Al3+  0.333333  0.666667  0.018947         3
## 12  O2-   0.3064    0         0.25            -2
## 13  O2-   0.6936    0.6936    0.25            -2
## 14  O2-   0         0.3064    0.25            -2
## 15  O2-   0.3064    0.3064    0.75            -2
## 16  O2-   0         0.6936    0.75            -2
## 17  O2-   0.6936    0         0.75            -2
## 18  O2-   0.973067  0.333333  0.583333        -2
## 19  O2-   0.360267  0.026933  0.583333        -2
## 20  O2-   0.666667  0.639733  0.583333        -2
## 21  O2-   0.973067  0.639733  0.083333        -2
## 22  O2-   0.666667  0.026933  0.083333        -2
## 23  O2-   0.360267  0.333333  0.083333        -2
## 24  O2-   0.639733  0.666667  0.916667        -2
## 25  O2-   0.026933  0.360267  0.916667        -2
## 26  O2-   0.333333  0.973067  0.916667        -2
## 27  O2-   0.639733  0.973067  0.416667        -2
## 28  O2-   0.333333  0.360267  0.416667        -2
## 29  O2-   0.026933  0.666667  0.416667        -2
##
##lattice (angstroms) =
##[[ 4.75400000e+00  0.00000000e+00  0.00000000e+00]
## [-2.37700000e+00  4.11708477e+00 -4.93038066e-32]
## [ 7.95408096e-16  1.37768724e-15  1.29900000e+01]]
##
##probe_position (angstroms) =
##[1.45662560e+00 3.44421809e-16 3.24750000e+00]
##
##V_ab_raw (V m^-2) =
##[[ 2.88260527e+20  1.43769600e+06 -1.11411200e+06]
## [ 1.43769600e+06 -1.71363640e+20 -1.29119619e+21]
## [-1.11411200e+06 -1.29119619e+21 -1.16896887e+20]]
##
##V_aa (V m^-2) =
##[ 2.88260527e+20  1.14735309e+21 -1.43561361e+21]
##
##principal_axes (unitless) =
##[[ 1.00000000e+00  2.09743158e-15 -1.43758398e-16]
## [-1.36466899e-15  6.99611715e-01  7.14523232e-01]
## [ 1.59923865e-15 -7.14523232e-01  6.99611715e-01]]
##
##eta (unitless) = (V_xx - V_yy)/V_zz =
##0.5984148884763557
##
##Q (m^2) =
##-2.65e-30
##
##nu_z (MHz) = 3eQV_zz/(2I(2I - 1)h)1e-6 =
##0.13798410650381887
##
##nu_Q (MHz) = nu_z(1 - eta^2/3)^0.5 =
##0.1294871268820299


###############################################################################
# time and run the calculation
start_time = time.time()

pyEFGPointCharge.calc_EFG_point_charge(input_pars, results)

end_time = time.time()

if not input_pars.delete_plotting_arrays:
    N_charges = results.charge_positions.shape[0]
    print('calc took', 
          round(end_time - start_time, 3), 
          's to run for', 
          '{:e}'.format(N_charges), 
          'sites')
else:
    print('calc took', 
          round(end_time - start_time, 3), 
          's to run')



###############################################################################
# print the results (all in Cartesian coordinates)

# print the results (all in Cartesian coordinates)
print('structure =')
print(results.structure)
print('')
print('lattice (angstroms) =')
print(results.lattice)
print('')
print('probe_position (angstroms) =')
print(results.probe_position)
print('')
print('V_ab_raw (V m^-2) =')
print(results.V_ab_raw)
print('')
print('V_aa (V m^-2) =')
print(results.V_aa)
print('')
print('principal_axes (unitless) =')
print(results.principal_axes)
print('')
print('eta (unitless) = (V_xx - V_yy)/V_zz =')
print(results.eta)
print('')
print('Q (m^2) =')
print(results.Q)
print('')
print('nu_z (MHz) = 3eQV_zz/(2I(2I - 1)h)1e-6 =')
print(results.nu_z_MHz)
print('')
print('nu_Q (MHz) = nu_z(1 - eta^2/3)^0.5 =')
print(results.nu_Q_MHz)
print('')



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
    pyEFGPointCharge.visualization(input_parameters=input_pars,
                             results=results,
                             relative_position_limits=pos_limits,
                             axes_range=None,
                             elev=0,
                             azim=0,
                             plot_sphere=True,
                             plot_lattice_vectors=False)
