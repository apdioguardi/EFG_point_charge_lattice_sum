
###############################################################################
# imports #####################################################################
###############################################################################

import numpy as np
from pymatgen.core import Structure
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D # this may not be required in python 3.9 and up
from isotopeDict import isotope_data_dict

###############################################################################
###############################################################################



###############################################################################
# data structures #############################################################
###############################################################################

# class to hold the input parameters and pass them to the function that 
# performs the calculation
class input_parameters:
    ##########################
    # general input parameters
    ##########################
    # path to cif file
    cif_file_name = None

    # create a supercell if one is required to define the charge structure
    # the supercell_dimensions parameter should be a list of form [Na, Nb, Nc],
    # where Na, Nb, and Nc are integers that represent the number of unit cells
    # in the crystallographic a, b, and c directions, respectively, contained in
    # the supercell. For example, [1, 1, 2] will create a supercell that is
    # 1 x 1 x 2 unit cells in size. Note: creating a supercell modifies the
    # lattice as well, therefore one must specify the charge_site_indices,
    # probe indices, etc. in the supercell.
    supercell_dimensions = [1, 1, 1]

    # choose the radius around which we will include spins in angstroms with 
    # respect to the probe position (float)
    sphere_radius = None

    # define the probe site number (nuclear spin site) (int)
    probe_site_index = None

    # for performing large calculations, we want to attempt to keep enough 
    # memory. To do this we'll trash some arrays that we don't need
    # running on 32bit python would cause memory allocation errors, so to 
    # make the program more portable, we can allow for that memory to be 
    # freed by the garbage collector by setting (no longer needed arrays 
    # to None)
    delete_plotting_arrays = True

    # probe nucleus identifying dictionary key (from isotope_data_dict) (string)
    # e.g. '11B'
    probe_nucleus = None

    # or if you want to set a custom nuclear quadrupole moment, do that here and
    # the default for the probe nucleus will be overridden (barns)
    manual_nuclear_quad_moment = None

    # in that case one also needs to define the nuclear spin I
    manual_nuclear_spin = None

    # define the charge sites by index in the cif file (list of ints)
    # Note: that site indices in pymatgen start at zero.
    # Further Note: if creating a supercell, one must include all charge
    # site indices in the supercell. therefore the user must first know the
    # site indices.
    charge_site_indices = None

    # charges should be set to be a numpy array of shape = (len(charge_site_indices),), 
    # the order of the charges should be the same as the charge_site_indices. units 
    # are in multiples of the elementary charge e.
    # Example:
    # define antiferromagnetic order on the  magnetic sites with indices 0 and 1,
    # so above charge_site_indices = [0, 1], and there is a charge of -2 on site 0
    # and a +1 charge on site 1, then charges = np.array([-2.0, 1.0])
    first_cell_charges = None

    # Sternheimer antishielding factor (\gamma_\inf)
    gamma_sternheimer = 0

    # number of Monte Carlo configurations to sample when disorder is present
    n_monte_carlo_samples = 1000

    # random seed for reproducibility of MC sampling (None for random)
    random_seed = None

    # correlated disorder groups. each entry is a dict with keys:
    #   'site_indices': list of site indices involved in this group
    #   'configurations': list of dicts, each with:
    #       'species': list of species strings or None (vacancy), aligned with site_indices
    #       'probability': float (optional; uniform assumed if omitted for all configs)
    # sites in a group are sampled together as a unit. a site may belong to at most
    # one group. sites in charge_site_indices but not in any group are sampled
    # independently based on their CIF occupancy. defaults to [] (no correlated disorder).
    # example:
    #   disorder_groups = [
    #       {
    #           'site_indices': [D, A, B, C],
    #           'configurations': [
    #               {'species': ['X', 'Al', None, None], 'probability': 0.5},
    #               {'species': ['Y', None, 'Si', 'Si'], 'probability': 0.5},
    #           ]
    #       }
    #   ]
    disorder_groups = []


# a class for returning results of the calculation
# all results are in cartesian coordinates!!!
class results:
    # the structure object created from the cif file
    structure = None
    # a matrix representing the the lattice unit cell as row vectors in 
    # cartesian coordinates (angstroms)
    lattice = None
    # probe position (angstroms)
    probe_position = None
    # Charge positions (angstroms)
    charge_positions = None
    # charges, in units of the elementary charge e
    charges = None
    # raw EFG tensor, ie not diagonalized (V m^-2)
    V_ab_raw = None
    # eigenvalues of the EFG tensor, ordered so |V_xx| <= |V_yy| <= |V_zz|; 
    # where x, y, and z no longer correspond to the 
    # cartesian axes, but indicate the directions of the principal axes of the
    # EFG tensor (V m^-2)
    V_aa = None
    # eigenvectors of the EFG tensor  (principal axes directions)
    # (normalized, cartesian coordinates, a matrix of column vectors, 
    # order is the same as the V_aa)
    principal_axes = None
    # eta = (V_xx - V_yy)/V_zz
    eta = None
    # nuclear quadrupole moment
    Q = None
    # nu_z = 3eQV_zz/(2I(2I-1)h)
    nu_z_MHz = None
    # nu_Q_MHz = nu_z_MHz(1 - eta**2/3)**0.5
    nu_Q_MHz = None
    # --- Monte Carlo disorder results (None for ordered structures) ---
    # eigenvalues for each MC sample, shape (n_mc, 3), ordered |V_xx|<=|V_yy|<=|V_zz|
    V_aa_samples = None
    # asymmetry parameter for each MC sample, shape (n_mc,)
    eta_samples = None
    # quadrupole frequency for each MC sample, shape (n_mc,); None if no Q given
    nu_Q_MHz_samples = None

###############################################################################
###############################################################################



###############################################################################
# functions ###################################################################
###############################################################################

def R_y(phi):
    """
    rotation matrix about the y-axis by angle phi. phi is expected to be in degrees.
    """
    return np.array([[ np.cos(phi*np.pi/180), 0, np.sin(phi*np.pi/180)],
                     [ 0,                     1, 0                    ],
                     [-np.sin(phi*np.pi/180), 0, np.cos(phi*np.pi/180)]])


def R_x(theta):
    """
    rotation matrix about the x-axis by angle theta. theta is expected to be in degrees.
    """
    return np.array([[1, 0,                        0                      ],
                     [0, np.cos(theta*np.pi/180), -np.sin(theta*np.pi/180)],
                     [0, np.sin(theta*np.pi/180),  np.cos(theta*np.pi/180)]])


def R_z(phi):
    """
    rotation matrix about the z-axis by angle phi. phi is expected to be in degrees.
    """
    return np.array([[np.cos(phi*np.pi/180), -np.sin(phi*np.pi/180), 0],
                     [np.sin(phi*np.pi/180),  np.cos(phi*np.pi/180), 0],
                     [0,                      0,                     1]])


def rotate_lattice(lattice_matrix):
    """
    function to rotate the lattice matrix from the default orientation, which 
    is that z and c are aligned, to a more reasonable default orientation, especially
    for dealing with monoclinic structures, but that is also consistent with triclinic
    structures. rotate such that $a||x$, $b$ in $x-y$ plane, and $c$ arbitrary
    we want to perform two rotations. first about the y-axis to align 
    the a direction with x. we'll get the angle via arctan for consistency
    with the second rotation. we could have just used 90 - beta, but 
    the second rotation about the x-axis will require some trig anyway,
    so lets just be consistent
    """
    # get original lattice vectors
    a_vec = lattice_matrix[0,:]
    b_vec = lattice_matrix[1,:]
    c_vec = lattice_matrix[2,:]
    # define rotation matrix about y that will bring a parallel to x
    R_y_a_par_x = R_y(np.arctan(a_vec[2]/a_vec[0])*180/np.pi)
    # perform first rotation about y
    a_vec = R_y_a_par_x @ a_vec
    b_vec = R_y_a_par_x @ b_vec
    c_vec = R_y_a_par_x @ c_vec
    # define rotation about the x-axis that will bring the b lattice vector
    # to lie in the xy plane
    R_x_b_in_xy = R_x(-np.arctan(b_vec[2]/b_vec[1])*180/np.pi)
    # perform second rotation about x
    a_vec = R_x_b_in_xy @ a_vec
    b_vec = R_x_b_in_xy @ b_vec
    c_vec = R_x_b_in_xy @ c_vec
    
    # return the rotated lattice matrix
    return np.vstack((a_vec, b_vec, c_vec))


def calc_EFG_point_charge(input_parameters, results):
    """
    Calculate the electric field gradient (EFG) tensor at a probe nucleus
    via a point-charge lattice-sum method.
    ############################################################################
    ahoy! see input parameter and results classes for details on the inputs to
    this function
    """
    # unpack the necessary attributes from the input_parameters class:
    cif_file_name = input_parameters.cif_file_name
    supercell_dimensions = input_parameters.supercell_dimensions
    sphere_radius = input_parameters.sphere_radius
    probe_site_index = input_parameters.probe_site_index
    delete_plotting_arrays = input_parameters.delete_plotting_arrays
    probe_nucleus = input_parameters.probe_nucleus
    manual_nuclear_quad_moment = input_parameters.manual_nuclear_quad_moment
    manual_nuclear_spin = input_parameters.manual_nuclear_spin
    charge_site_indices = input_parameters.charge_site_indices
    first_cell_charges = input_parameters.first_cell_charges
    gamma_sternheimer = input_parameters.gamma_sternheimer
    n_monte_carlo_samples = input_parameters.n_monte_carlo_samples
    random_seed = input_parameters.random_seed
    disorder_groups = input_parameters.disorder_groups

    # load cif file and create a pymatgen.core.structure.Structure object with
    # the info from said file.
    structure = Structure.from_file(cif_file_name)
    structure.make_supercell(supercell_dimensions)

    # modify the lattices of the structure and magnetic structure to the default
    # orientation such that $a \parallel x$, $b$ in $x-y$ plane, and $c$ arbitrary
    # we want to perform two rotations. first about the y-axis to align
    # the a direction with x. we'll get the angle via arctan for consistency
    # with the second rotation. we could have just used 90 - beta, but
    # the second rotation about the x-axis will require some trig anyway,
    # so lets just be consistent
    rotated_lattice = rotate_lattice(structure.lattice.matrix)
    # set the new lattice matrix
    structure.lattice = rotated_lattice
    # set the new structure and lattice into the results
    results.structure = structure
    results.lattice = rotated_lattice

    # define constants to get the units correct
    epsilon_0 = 8.8541878128e-12 # permittivity of free space. units: F m^-1 = kg^-1 m^-3 s^4 A^2
    e = 1.602176634e-19 # elementary charge. units: C = A s
    h = 6.62607015e-34 # Planck constant. units: J s

    # resolve nuclear quadrupole moment and spin
    Q = None
    I = None
    if probe_nucleus is not None:
        Q = isotope_data_dict[probe_nucleus]["Q"]*1e-28 #1 barn = 10^-28 m^2
        I = isotope_data_dict[probe_nucleus]["I0"]
    elif manual_nuclear_quad_moment is not None:
        Q = manual_nuclear_quad_moment*1e-28 #1 barn = 10^-28 m^2
        I = manual_nuclear_spin

    ###########################################################################
    # disorder setup: tag each charge site with its sampling info
    ###########################################################################

    # build a set of site indices that belong to a disorder group so we can
    # distinguish them from independently-disordered sites
    grouped_site_indices = set()
    for group in disorder_groups:
        for idx in group['site_indices']:
            grouped_site_indices.add(idx)

    # validate disorder groups
    for g, group in enumerate(disorder_groups):
        g_indices = group['site_indices']
        g_configs = group['configurations']
        for idx in g_indices:
            if idx not in charge_site_indices:
                raise ValueError(
                    f"disorder_groups[{g}]: site index {idx} is not in charge_site_indices")
        seen = set()
        for idx in g_indices:
            if idx in seen:
                raise ValueError(
                    f"disorder_groups[{g}]: site index {idx} appears more than once")
            seen.add(idx)
        n_sites_in_group = len(g_indices)
        for c, config in enumerate(g_configs):
            if len(config['species']) != n_sites_in_group:
                raise ValueError(
                    f"disorder_groups[{g}] config {c}: 'species' length "
                    f"({len(config['species'])}) != number of site_indices ({n_sites_in_group})")
            for k, sp in enumerate(config['species']):
                if sp is not None:
                    site_idx = g_indices[k]
                    charge_entry = first_cell_charges[site_idx]
                    if not isinstance(charge_entry, dict) or sp not in charge_entry:
                        raise ValueError(
                            f"disorder_groups[{g}] config {c}: species '{sp}' at site "
                            f"{site_idx} not found in first_cell_charges[{site_idx}]")
        probs = [cfg.get('probability', None) for cfg in g_configs]
        if any(p is not None for p in probs):
            total = sum(p for p in probs if p is not None)
            if abs(total - 1.0) > 1e-6:
                raise ValueError(
                    f"disorder_groups[{g}]: probabilities sum to {total}, expected 1.0")

    # also validate that no site index appears in more than one group
    all_grouped = []
    for g, group in enumerate(disorder_groups):
        for idx in group['site_indices']:
            if idx in all_grouped:
                raise ValueError(
                    f"Site index {idx} appears in more than one disorder group")
            all_grouped.append(idx)

    # determine whether any disorder is present (from groups or independent sites)
    has_group_disorder = len(disorder_groups) > 0
    has_independent_disorder = any(
        idx not in grouped_site_indices and not structure[idx].is_ordered
        for idx in charge_site_indices
    )
    has_disorder = has_group_disorder or has_independent_disorder

    # tag each charge site with sampling metadata stored in site.properties so
    # the info propagates through get_sites_in_sphere
    for index in charge_site_indices:
        site = structure[index]
        charge_entry = first_cell_charges[index]

        if index in grouped_site_indices:
            # group-sampled sites: find which group and position within it,
            # store a reference so the MC loop can look up the config
            for g, group in enumerate(disorder_groups):
                if index in group['site_indices']:
                    pos_in_group = group['site_indices'].index(index)
                    site.properties['_disorder_type'] = 'group'
                    site.properties['_group_index'] = g
                    site.properties['_pos_in_group'] = pos_in_group
                    # placeholder charge (will be overwritten per MC step)
                    site.properties['charge'] = 0.0
                    break

        elif isinstance(charge_entry, dict):
            # independently-disordered site: build species->charge map and
            # occupancy weights from the CIF
            species_charges = charge_entry  # e.g. {'Al': +3.0, 'Si': +4.0}
            occ_weights = {}
            for sp, occ in site.species.items():
                sp_str = sp.symbol
                if sp_str in species_charges:
                    occ_weights[sp_str] = float(occ)
            total_occ = sum(occ_weights.values())
            if total_occ < 1.0 - 1e-6:
                # implicit vacancy
                occ_weights['_vacancy'] = 1.0 - total_occ
            charge_map = dict(species_charges)
            charge_map['_vacancy'] = 0.0
            site.properties['_disorder_type'] = 'independent'
            site.properties['_species_charges'] = charge_map
            site.properties['_occ_weights'] = occ_weights
            # placeholder charge
            site.properties['charge'] = 0.0

        else:
            # fully ordered site: just store the scalar charge
            site.properties['_disorder_type'] = 'ordered'
            site.properties['charge'] = float(charge_entry)

    ###########################################################################
    # geometry: computed once, independent of charges
    ###########################################################################

    # get the position in cartesian coordinates of the probe site
    probe_position = structure[probe_site_index].coords

    #500 angstrom radius: executed in 13.0s, finished 15:54:31 2021-08-03
    sites_in_sphere = structure.get_sites_in_sphere(probe_position, sphere_radius)

    # loop over the sites and extract positions of sites in cartesian
    # coordinates. Couldn't find a way to make this faster than by using a list
    # comprehension. it is possible that one could use some numpy vectorization magic
    # here but for now it is fast enough for our purposes
    positions = np.array([site.coords for site in sites_in_sphere])
    # calc relative positions to the probe position
    results.probe_position = probe_position
    relative_positions = positions - probe_position

    # calculate the relative distances
    relative_distances = np.linalg.norm(relative_positions, axis=1)
    # get the indices of the nonzero relative distances, as we don't want to include
    # the probe-site charge (although numpy properly deals with this with just warnings and nan values)
    # this is actually memory inefficient, because the index array is the same size as the relative distance
    # array so if we run into memory issues, could try a different method here... this was just the first
    # idea that i had
    nonzero_indices = np.nonzero(relative_distances)
    # retain just the relative distances and positions with nonzero values
    relative_distances = relative_distances[nonzero_indices]
    positions = positions[nonzero_indices]

    # set the orignal positions array to None to free memory
    if delete_plotting_arrays:
        positions = None
    else:
        results.charge_positions = positions
        positions = None

    # extract the x, y, and z (cartesian coordinates) components of the r_vecs
    relative_positions_x = relative_positions[:, 0].flatten()
    relative_positions_y = relative_positions[:, 1].flatten()
    relative_positions_z = relative_positions[:, 2].flatten()
    # retain just the relative positions except for the probe site
    relative_positions_x = relative_positions_x[nonzero_indices]
    relative_positions_y = relative_positions_y[nonzero_indices]
    relative_positions_z = relative_positions_z[nonzero_indices]

    # delete this array that is no longer needed
    relative_positions = None
    # note: nonzero_indices is kept alive until after qi / sphere_props are filtered below

    # calculate the independent EFG tensor geometric factors (charge-independent)
    Vi_xx = (3*relative_positions_x**2 - relative_distances**2)/relative_distances**5
    Vi_xy = 3*relative_positions_x*relative_positions_y/relative_distances**5
    Vi_xz = 3*relative_positions_x*relative_positions_z/relative_distances**5
    Vi_yy = (3*relative_positions_y**2 - relative_distances**2)/relative_distances**5
    Vi_yz = 3*relative_positions_y*relative_positions_z/relative_distances**5
    Vi_zz = (3*relative_positions_z**2 - relative_distances**2)/relative_distances**5

    # delete arrays that are no longer needed
    relative_distances = None
    relative_positions_x = None
    relative_positions_y = None
    relative_positions_z = None

    # stack geometric factors into shape (N_sites, 3, 3) with SI prefactor baked in.
    # units: V m^-2 per unit charge e, with positions in Angstroms converted to meters.
    # scale = e / (4 pi epsilon_0) * (1e10)^3  [Angstrom^-3 -> m^-3]
    scale = 1e10**3 * e / (4*np.pi*epsilon_0)
    Vi_geom = np.column_stack((Vi_xx, Vi_xy, Vi_xz,
                               Vi_xy, Vi_yy, Vi_yz,
                               Vi_xz, Vi_yz, Vi_zz)).reshape(-1, 3, 3) * scale

    Vi_xx = None; Vi_xy = None; Vi_xz = None
    Vi_yy = None; Vi_yz = None; Vi_zz = None

    ###########################################################################
    # helper: compute EFG tensor from a charges array
    ###########################################################################

    def _compute_V_ab(qi):
        Vi = Vi_geom * qi[:, np.newaxis, np.newaxis]
        return np.sum(np.nan_to_num(Vi), axis=0) * (1 - gamma_sternheimer)

    def _diagonalize(V_ab_raw):
        evals, evecs = np.linalg.eigh(V_ab_raw)
        idx = np.argsort(np.abs(evals))
        V_aa = evals[idx]
        principal_axes = evecs[:, idx]
        eta = (V_aa[0] - V_aa[1]) / V_aa[2]
        return V_aa, principal_axes, eta

    ###########################################################################
    # ordered-only path (no disorder): single calculation, fully backwards
    # compatible with the original function behaviour
    ###########################################################################

    if not has_disorder:
        # extract charges from tagged sites (all ordered), filtering out the probe site
        qi = np.array([site.properties['charge'] for site in sites_in_sphere])
        qi = qi[nonzero_indices]
        nonzero_indices = None
        sites_in_sphere = None

        V_ab_raw = _compute_V_ab(qi)
        results.V_ab_raw = V_ab_raw

        if delete_plotting_arrays:
            qi = None
        else:
            results.charges = qi
            qi = None

        V_aa, principal_axes, eta = _diagonalize(V_ab_raw)
        results.V_aa = V_aa
        results.principal_axes = principal_axes
        results.eta = eta

        if Q is not None:
            nu_z_MHz = 3*e*Q*V_aa[2]/(2*I*(2*I - 1)*h)*1e-6
            nu_Q_MHz = nu_z_MHz*(1 - eta**2/3)**0.5
            results.Q = Q
            results.nu_z_MHz = nu_z_MHz
            results.nu_Q_MHz = nu_Q_MHz

        return results

    ###########################################################################
    # disordered path: Monte Carlo sampling
    ###########################################################################

    rng = np.random.default_rng(random_seed)

    # build a list of site property dicts, filtered to non-probe sites only
    # (same ordering as Vi_geom rows, which used nonzero_indices)
    sphere_props = [sites_in_sphere[i].properties for i in nonzero_indices[0]]
    nonzero_indices = None
    sites_in_sphere = None

    # pre-build per-group probability arrays (normalised) for fast sampling
    group_prob_arrays = []
    for group in disorder_groups:
        probs = [cfg.get('probability', None) for cfg in group['configurations']]
        if all(p is None for p in probs):
            n = len(group['configurations'])
            probs = np.ones(n) / n
        else:
            probs = np.array([p if p is not None else 0.0 for p in probs], dtype=float)
            probs /= probs.sum()
        group_prob_arrays.append(probs)

    # number of sites in the sphere (excluding probe)
    n_sphere_sites = Vi_geom.shape[0]

    V_aa_list = []
    eta_list = []
    nu_Q_list = []
    # store last sample's qi and V_ab_raw for results.charges / results.V_ab_raw
    qi_last = None
    V_ab_raw_last = None

    for _ in range(n_monte_carlo_samples):
        # for each MC step, first sample a configuration index for every group
        group_config_indices = [
            int(rng.choice(len(disorder_groups[g]['configurations']),
                           p=group_prob_arrays[g]))
            for g in range(len(disorder_groups))
        ]

        # build qi for this sample
        qi_sample = np.zeros(n_sphere_sites)
        for i, props in enumerate(sphere_props):
            dtype = props.get('_disorder_type', None)
            if dtype == 'ordered' or dtype is None:
                qi_sample[i] = props['charge']
            elif dtype == 'group':
                g = props['_group_index']
                pos = props['_pos_in_group']
                config_idx = group_config_indices[g]
                sp = disorder_groups[g]['configurations'][config_idx]['species'][pos]
                if sp is None:
                    qi_sample[i] = 0.0
                else:
                    site_index = disorder_groups[g]['site_indices'][pos]
                    qi_sample[i] = first_cell_charges[site_index][sp]
            elif dtype == 'independent':
                charge_map = props['_species_charges']
                occ_weights = props['_occ_weights']
                symbols = list(occ_weights.keys())
                weights = np.array([occ_weights[s] for s in symbols])
                chosen = symbols[int(rng.choice(len(symbols), p=weights/weights.sum()))]
                qi_sample[i] = charge_map.get(chosen, 0.0)

        V_ab_raw = _compute_V_ab(qi_sample)
        V_aa, principal_axes, eta = _diagonalize(V_ab_raw)

        V_aa_list.append(V_aa)
        eta_list.append(eta)
        if Q is not None:
            nu_z = 3*e*Q*V_aa[2]/(2*I*(2*I - 1)*h)*1e-6
            nu_Q_list.append(nu_z*(1 - eta**2/3)**0.5)

        qi_last = qi_sample
        V_ab_raw_last = V_ab_raw

    # store sample distributions
    results.V_aa_samples = np.array(V_aa_list)
    results.eta_samples = np.array(eta_list)

    # scalar results are means over the MC ensemble
    results.V_aa = np.mean(results.V_aa_samples, axis=0)
    results.eta = float(np.mean(results.eta_samples))
    # principal axes from the last sample (representative orientation)
    results.principal_axes = principal_axes
    results.V_ab_raw = V_ab_raw_last

    if Q is not None:
        results.nu_Q_MHz_samples = np.array(nu_Q_list)
        results.nu_Q_MHz = float(np.mean(results.nu_Q_MHz_samples))
        # nu_z_MHz from the mean V_zz
        nu_z_mean = 3*e*Q*results.V_aa[2]/(2*I*(2*I - 1)*h)*1e-6
        results.nu_z_MHz = nu_z_mean
        results.Q = Q

    if not delete_plotting_arrays:
        results.charges = qi_last

    return results


###############################################################################
###############################################################################



###############################################################################
# visualization ###############################################################
###############################################################################

def visualization(input_parameters,
                  results, 
                  relative_position_limits=None,
                  axes_range=None,
                  elev=0., 
                  azim=0., 
                  plot_sphere=False,
                  plot_lattice_vectors=False,
                  charge_scale_factor=50,
                  show_labels=True):
    """
    plot the lattice, and EFG tensor elements, with bells and whistles
    input_parameters and results
        input and output class data structures from the EFG tensor calculation
    relative_position_limits
        a tuple of the form: (x_min_rel, x_max_rel, y_min_rel, y_max_rel, 
        z_min_rel, z_max_rel) where the values are relative to to the probe
        position
    axes_range
        tuple of (x_min, x_max, y_min, y_max, z_min, z_max)
    elev and azim
        angles describing the initial viewpoint 
        of the visualization
    plot_sphere
        boolean variable that toggles the plotting of the sphere over which
        the calculation was performed
    plot_lattice_vectors
        boolean variable that toggles plotting the lattice vectors (of only
        the first lattice)
    charge_scale_factor
        used to scale the value of charge to the size on the plot for
        help of visualization. plotted size is proportional to the absolute
        value of the charge on that site
    show_labels
        show the charge values on the plot
    
    Note that the visualization is only useful when plotting less than aprox
    50 sites. If more sites are desired, one should use mayavi, but because
    this package is difficult to install due to dependencies and problems 
    with vtk istallation, we decided to use this matplotlib setup for ease
    of use.
    """
    # unpack the results for the visualization
    probe_position = results.probe_position
    positions = results.charge_positions
    charges = results.charges
    lattice = results.lattice
    V_aa = results.V_aa
    principal_axes = results.principal_axes
    sphere_radius = input_parameters.sphere_radius

    probe_position_x = probe_position[0]
    probe_position_y = probe_position[1]
    probe_position_z = probe_position[2]

    # make 1D array copies from all the input arrays (so we don't change the input data)
    positions_x = np.copy(positions[:,0])
    positions_y = np.copy(positions[:,1])
    positions_z = np.copy(positions[:,2])

    charges = np.copy(charges)

    # if we have limits apply them
    if relative_position_limits is not None:
        # unpack the limits
        x_min_rel, x_max_rel, y_min_rel, y_max_rel, z_min_rel, z_max_rel = relative_position_limits
        x_min = x_min_rel + probe_position_x
        x_max = x_max_rel + probe_position_x
        y_min = y_min_rel + probe_position_y
        y_max = y_max_rel + probe_position_y
        z_min = z_min_rel + probe_position_z
        z_max = z_max_rel + probe_position_z
        # generate the boolean filter arrays
        vecs_filter_x = np.logical_and((positions_x > x_min), (positions_x < x_max))
        vecs_filter_y = np.logical_and((positions_y > y_min), (positions_y < y_max))
        vecs_filter_z = np.logical_and((positions_z > z_min), (positions_z < z_max))
        vecs_filter = np.logical_and(np.logical_and(vecs_filter_x, vecs_filter_y), vecs_filter_z)
        # apply the filters to the position vectors and charges
        positions_x = positions_x[vecs_filter]
        positions_y = positions_y[vecs_filter]
        positions_z = positions_z[vecs_filter]
        charges = charges[vecs_filter]


    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')

    if plot_lattice_vectors:
        # make unit cell vectors and plot them:
        a = lattice[0,:]
        b = lattice[1,:]
        c = lattice[2,:]
        ax.quiver(0, 0, 0, a[0], a[1], a[2], color='red')
        ax.quiver(0, 0, 0, b[0], b[1], b[2], color='green')
        ax.quiver(0, 0, 0, c[0], c[1], c[2], color='blue')
        
    # draw the probe position
    ax.scatter(probe_position_x, 
               probe_position_y, 
               probe_position_z)
    
    # draw the EFG tensor at the probe position
    rgb_list = ['red', 'green', 'blue']
    for i in range(3):
        ax.quiver(probe_position_x, 
                  probe_position_y, 
                  probe_position_z,
                  V_aa[i]*principal_axes[0, i],
                  V_aa[i]*principal_axes[1, i],
                  V_aa[i]*principal_axes[2, i],
                  color=rgb_list[i],
                  length=3,
                  normalize=True)
    
    # draw the lattice and charges
    ax.scatter(positions_x,
               positions_y,
               positions_z,
               c=charges,
               cmap='bwr',
               s=charge_scale_factor*np.abs(charges))
    # add labels if desired
    if show_labels:
        it = np.nditer(charges, flags=['multi_index'])
        for charge in it:
            label = '$~~$%+d' % charge
            ax.text(positions_x[it.multi_index[0]],
                    positions_y[it.multi_index[0]],
                    positions_z[it.multi_index[0]],
                    label)

    # Create a sphere and plot for double check
    if plot_sphere:
        r = sphere_radius
        phi, theta = np.mgrid[0.0:np.pi:25j, 0.0:2.0*np.pi:50j] # note that fewer points will make it faster to rotate
        x = r*np.sin(phi)*np.cos(theta) + probe_position_x
        y = r*np.sin(phi)*np.sin(theta) + probe_position_y
        z = r*np.cos(phi) + probe_position_z
        ax.plot_surface(x, y, z,  rstride=1, cstride=1, color='c', alpha=0.15, linewidth=0)

    # label the axes
    ax.set_xlabel('x (angstroms)')
    ax.set_ylabel('y (angstroms)')
    ax.set_zlabel('z (angstroms)')
    
    # axes limits if given
    if axes_range is not None:
        # unpack the axes limits limits
        x_axis_min, x_axis_max, y_axis_min, y_axis_max, z_axis_min, z_axis_max = axes_range
        ax.set_xlim3d(x_axis_min, x_axis_max)
        ax.set_ylim3d(y_axis_min, y_axis_max)
        ax.set_zlim3d(z_axis_min, z_axis_max)
        
    # set the initial viewing angle
    ax.view_init(elev=elev, azim=azim)

    plt.show()

###############################################################################
###############################################################################

