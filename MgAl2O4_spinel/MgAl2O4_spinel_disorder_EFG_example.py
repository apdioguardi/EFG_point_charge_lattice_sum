###############################################################################
# MgAl2O4 spinel -- disordered EFG point-charge calculation example
###############################################################################
#
# Material: MgAl2O4 (magnesia spinel)
# Structure: Redfern et al. (1999) Am. Min. 84, 299-310
#            In-situ neutron diffraction at T = 673 K (heating cycle)
#            COD entry 9002084
#
# Disorder: cation inversion -- a fraction i of Al3+ occupies the normally
#           Mg2+ tetrahedral (8a) site, and the same fraction i of Mg2+
#           occupies the normally Al3+ octahedral (16d) site.
#           At 673 K: i = 0.155 (from Al1 occupancy at the 8a site)
#
#           8a  (tetrahedral): Mg (occ=0.845) + Al (occ=0.155)
#           16d (octahedral):  Al (occ=0.9225) + Mg (occ=0.0775)
#           32e (oxygen):      O  (occ=1.000)  -- ordered
#
# Probe nucleus: 27Al (I = 5/2, Q = 0.1466 barns)
#
# Experimental 27Al NMR reference values (from literature):
#   Octahedral Al [16d] in near-normal spinel:
#       Cq ~ 3.3-3.7 MHz, eta ~ 0  (D3d local symmetry, ordered limit)
#   Tetrahedral Al [8a] (inverted Al):
#       Cq ~ 1.6 MHz,    eta ~ 0  (S4 local symmetry, ordered limit)
#   Both sites show a DISTRIBUTION of Cq when disorder is present,
#   extending up to ~5-7 MHz for octahedral sites near inverted neighbors.
#   References:
#       Millard et al. (1992) Am. Min. 77, 44-52
#       Kashii et al. (1999) J. Am. Ceram. Soc. 82, 1229-1235
#
###############################################################################

import sys
import os
import numpy as np
import matplotlib.pyplot as plt

# add the parent directory to the path so pyEFGPointCharge is importable
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
import pyEFGPointCharge

###############################################################################
# Step 0: explore the structure to identify site indices
###############################################################################
# Run this block once to learn the site layout, then comment it out.

from pymatgen.core import Structure

cif_path = os.path.join(os.path.dirname(__file__),
                        'Redfern_1999_MgAl2O4_673K_COD9002084.cif')

structure_explore = Structure.from_file(cif_path)
print(f"\nStructure loaded: {structure_explore.formula}")
print(f"Number of sites: {len(structure_explore)}\n")
print(f"{'Index':>5}  {'Species':>30}  {'is_ordered':>10}  Frac coords")
for i, site in enumerate(structure_explore):
    print(f"{i:>5}  {str(site.species):>30}  {str(site.is_ordered):>10}  "
          f"{site.frac_coords}")

# Expected output (pymatgen Fd-3m, conventional cell, 56 sites):
#   Sites 0-7:   8a tetrahedral -- disordered {Mg0+: 0.845, Al0+: 0.155}
#   Sites 8-23:  16d octahedral -- disordered {Al0+: 0.9225, Mg0+: 0.0775}
#   Sites 24-55: 32e oxygen     -- ordered    {O0+: 1.0}
# NOTE: the exact indexing depends on pymatgen version; verify by running.

###############################################################################
# Step 1: define site index ranges (adjust if needed based on Step 0 output)
###############################################################################

n_tet = 8    # number of 8a tetrahedral sites
n_oct = 16   # number of 16d octahedral sites
n_O   = 32   # number of 32e oxygen sites

tet_indices = list(range(0, n_tet))               # 8a  sites: 0-7
oct_indices = list(range(n_tet, n_tet + n_oct))   # 16d sites: 8-23
O_indices   = list(range(n_tet + n_oct,
                         n_tet + n_oct + n_O))    # 32e sites: 24-55

charge_site_indices = tet_indices + oct_indices + O_indices

###############################################################################
# Step 2: assign charges
###############################################################################
# first_cell_charges is a list indexed by site index.
# For ordered (O) sites: scalar charge.
# For disordered (Mg/Al) sites: dict keyed by element symbol.

first_cell_charges = [None] * (n_tet + n_oct + n_O)

# tetrahedral (8a) sites: normally Mg2+, minority Al3+ from inversion
for i in tet_indices:
    first_cell_charges[i] = {'Mg': +2.0, 'Al': +3.0}

# octahedral (16d) sites: normally Al3+, minority Mg2+ from inversion
for i in oct_indices:
    first_cell_charges[i] = {'Al': +3.0, 'Mg': +2.0}

# oxygen (32e) sites: fully ordered O2-
for i in O_indices:
    first_cell_charges[i] = -2.0

###############################################################################
# Example A: probe = Al on an octahedral (16d) site
###############################################################################
# Pick site index 8 (first 16d site). This site has 92.25% Al, 7.75% Mg.
# The probe IS the Al nucleus on this site; the surrounding lattice is disordered.

probe_site_index_oct = n_tet   # = 8 (first 16d site)

###############################################################################
# Example B: probe = Al on a tetrahedral (8a) site (inverted Al)
###############################################################################
# Pick site index 0 (first 8a site). Only 15.5% of these sites have Al,
# but when they do, the EFG is what we want to capture.

probe_site_index_tet = 0

###############################################################################
# Part 1: independent disorder (sites sampled independently from CIF occupancies)
###############################################################################
# This is the simplest model: each disordered site independently draws a
# species at random based on its CIF occupancy fractions.
# No disorder_groups needed; pyEFGPointCharge detects disorder automatically
# from structure[i].is_ordered.

print("\n" + "="*70)
print("PART 1: Independent disorder -- octahedral Al probe [16d site]")
print("="*70)

input_pars_indep = pyEFGPointCharge.input_parameters()
input_pars_indep.cif_file_name        = cif_path
input_pars_indep.sphere_radius        = 15.0      # Angstroms
input_pars_indep.probe_site_index     = probe_site_index_oct
input_pars_indep.charge_site_indices  = charge_site_indices
input_pars_indep.first_cell_charges   = first_cell_charges
input_pars_indep.probe_nucleus        = '27Al'
input_pars_indep.n_monte_carlo_samples = 5000
input_pars_indep.random_seed          = 42        # for reproducibility
input_pars_indep.delete_plotting_arrays = True

results_indep = pyEFGPointCharge.results()
pyEFGPointCharge.calc_EFG_point_charge(input_pars_indep, results_indep)

print(f"\nMean V_zz:    {results_indep.V_aa[2]:.4e} V/m^2")
print(f"Mean eta:     {results_indep.eta:.4f}")
print(f"Mean nu_Q:    {results_indep.nu_Q_MHz:.4f} MHz")
print(f"\nDistribution (n={input_pars_indep.n_monte_carlo_samples} samples):")
print(f"  |nu_Q| mean  = {np.mean(np.abs(results_indep.nu_Q_MHz_samples)):.4f} MHz")
print(f"  |nu_Q| std   = {np.std(np.abs(results_indep.nu_Q_MHz_samples)):.4f} MHz")
print(f"  |nu_Q| min   = {np.min(np.abs(results_indep.nu_Q_MHz_samples)):.4f} MHz")
print(f"  |nu_Q| max   = {np.max(np.abs(results_indep.nu_Q_MHz_samples)):.4f} MHz")
print(f"  eta   mean   = {np.mean(np.abs(results_indep.eta_samples)):.4f}")
print(f"  eta   std    = {np.std(results_indep.eta_samples):.4f}")

print("\nLiterature (Millard 1992, Kashii 1999):")
print("  Octahedral Al: Cq ~ 3.3-3.7 MHz (ordered), up to ~7 MHz (disordered)")
print("  eta ~ 0 (ordered), 0-0.4 (disordered)")

###############################################################################
# Part 2: correlated disorder -- pairs of (8a Al, 16d Mg) are coupled
###############################################################################
# Physical motivation: when Al inverts from 16d to 8a, it swaps with Mg.
# The pair (one 8a site, one neighbouring 16d site) forms a correlated unit:
#   Config 0 (prob = 1-i): 8a = Mg,  16d = Al   (normal spinel)
#   Config 1 (prob = i):   8a = Al,  16d = Mg    (inverted pair)
#
# Here i = 0.155 (from CIF at 673 K). We demonstrate with one correlated pair
# of adjacent sites (site 0 paired with site 8). In a real calculation you
# would define all 8 such pairs across the unit cell.
#
# NOTE: in the correlated model the marginal occupancies of each site
# individually reproduce the CIF values exactly:
#   P(8a = Al)  = i = 0.155
#   P(16d = Mg) = i = 0.155   (but note: 16d has 16 sites and 8a has 8,
#    so stoichiometry requires each 8a-Al is correlated with 2 16d-Mg sites
#    if the pairing is 1:2. For simplicity this example uses 1:1 pairs.)

inversion_parameter = 0.155   # = Al1 occupancy at 8a from CIF

print("\n" + "="*70)
print("PART 2: Correlated disorder -- octahedral Al probe [16d site]")
print("="*70)

input_pars_corr = pyEFGPointCharge.input_parameters()
input_pars_corr.cif_file_name        = cif_path
input_pars_corr.sphere_radius        = 15.0
input_pars_corr.probe_site_index     = probe_site_index_oct
input_pars_corr.charge_site_indices  = charge_site_indices
input_pars_corr.first_cell_charges   = first_cell_charges
input_pars_corr.probe_nucleus        = '27Al'
input_pars_corr.n_monte_carlo_samples = 5000
input_pars_corr.random_seed          = 42

# Define correlated pairs: each 8a site (index 0-7) is paired with one 16d
# site (index 8-15) such that when the 8a has Al, the paired 16d has Mg.
# Stoichiometry is handled approximately here (1:1 pairing).
i = inversion_parameter

disorder_groups = []
for tet_idx, oct_idx in zip(tet_indices[:8], oct_indices[:8]):
    disorder_groups.append({
        'site_indices': [tet_idx, oct_idx],
        'configurations': [
            {'species': ['Mg', 'Al'], 'probability': 1.0 - i},  # normal
            {'species': ['Al', 'Mg'], 'probability': i},         # inverted
        ]
    })

# The remaining 8 octahedral sites (indices 16-23) are not paired and sample
# independently. To keep them at the right average occupancy, we leave them
# as independently-disordered (which uses the CIF occupancy automatically).

input_pars_corr.disorder_groups = disorder_groups

results_corr = pyEFGPointCharge.results()
pyEFGPointCharge.calc_EFG_point_charge(input_pars_corr, results_corr)

print(f"\nMean V_zz:    {results_corr.V_aa[2]:.4e} V/m^2")
print(f"Mean eta:     {results_corr.eta:.4f}")
print(f"Mean nu_Q:    {results_corr.nu_Q_MHz:.4f} MHz")
print(f"\nDistribution (n={input_pars_corr.n_monte_carlo_samples} samples):")
print(f"  |nu_Q| mean  = {np.mean(np.abs(results_corr.nu_Q_MHz_samples)):.4f} MHz")
print(f"  |nu_Q| std   = {np.std(np.abs(results_corr.nu_Q_MHz_samples)):.4f} MHz")

###############################################################################
# Part 3: plot the nu_Q distributions
###############################################################################

fig, axes = plt.subplots(1, 2, figsize=(12, 5))

for ax, res, label in zip(axes,
                           [results_indep, results_corr],
                           ['Independent disorder', 'Correlated disorder (paired groups)']):
    nu_Q = np.abs(res.nu_Q_MHz_samples)
    ax.hist(nu_Q, bins=80, density=True, color='steelblue', edgecolor='none', alpha=0.8)
    ax.axvline(np.mean(nu_Q), color='red', lw=1.5, linestyle='--', label=f'Mean = {np.mean(nu_Q):.2f} MHz')
    # literature range for octahedral Al
    ax.axvspan(3.3, 3.7, alpha=0.15, color='green', label='Lit. ordered Cq (3.3-3.7 MHz)')
    ax.set_xlabel(r'$|\nu_Q|$ (MHz)', fontsize=12)
    ax.set_ylabel('Probability density', fontsize=12)
    ax.set_title(f'27Al EFG distribution [{label}]\nOctahedral (16d) probe site, MgAl$_2$O$_4$ 673 K', fontsize=10)
    ax.legend(fontsize=9)

plt.tight_layout()
plt.savefig(os.path.join(os.path.dirname(__file__), 'MgAl2O4_nuQ_distribution.png'),
            dpi=150, bbox_inches='tight')
plt.show()

print("\nDone. Distribution plot saved to MgAl2O4_nuQ_distribution.png")
