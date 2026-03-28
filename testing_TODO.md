# Testing TODO -- pyEFGPointCharge disorder support

## Test material: MgAl2O4 spinel
**CIF:** `MgAl2O4_spinel/Redfern_1999_MgAl2O4_673K_COD9002084.cif`
**Source:** Redfern et al. (1999) *Am. Min.* 84, 299-310, COD entry 9002084
**Disorder:** Cation inversion -- Mg2+/Al3+ swap between 8a (tetrahedral) and 16d (octahedral) sites
**Inversion parameter** at 673 K: i = 0.155 (Al occupancy at 8a)
**Probe nucleus:** 27Al (I = 5/2, Q = 0.1466 barns, key `'27Al'` in `isotopeDict`)
**Example script:** `MgAl2O4_spinel/MgAl2O4_spinel_disorder_EFG_example.py`

---

## Literature reference values for comparison

| Site | Probe | Cq (MHz) | eta | Reference |
|------|-------|-----------|-----|-----------|
| Tetrahedral Al (8a, inverted) | 27Al | ~1.6 | ~0 | Millard et al. 1992 *Am. Min.* 77, 44 |
| Octahedral Al (16d, normal) | 27Al | ~3.3-3.7 | ~0 | Millard et al. 1992; Kashii et al. 1999 |
| Octahedral Al near inverted neighbor | 27Al | ~5-7 (broad) | 0-0.4 | Kashii et al. 1999 |

Note: Cq = e*Q*V_zz / h is the quadrupole coupling constant.
The code computes nu_z = 3*e*Q*V_zz / (2I(2I-1)*h) = 3*Cq / (2I(2I-1)),
and nu_Q = nu_z * (1 - eta^2/3)^(1/2).
For I = 5/2: nu_z = 3*Cq/20, so Cq = 20*nu_z/3.

---

## Test checklist

### 1. Backwards compatibility -- ordered structure
- [ ] Run `aAl2O3_EFG_point_charge_calc_example.py` and verify output is **unchanged** from before the disorder refactor
- [ ] Check that `results.V_aa_samples` is `None`, `results.eta_samples` is `None`
- [ ] Check `results.V_aa`, `results.eta`, `results.nu_Q_MHz` are populated and match prior values

### 2. Structure loading and site identification
- [ ] Load `Redfern_1999_MgAl2O4_673K_COD9002084.cif` with pymatgen and print site list
- [ ] Confirm total site count = 56 (8 tet + 16 oct + 32 O in the conventional cell)
- [ ] Confirm sites 0-7 are disordered (`site.is_ordered == False`) with species containing Mg and Al
- [ ] Confirm sites 8-23 are disordered with species containing Al and Mg
- [ ] Confirm sites 24-55 are ordered O2-
- [ ] If site ordering differs, update index ranges in the example script accordingly

### 3. Independent disorder -- basic function
- [ ] Run `MgAl2O4_spinel_disorder_EFG_example.py` Part 1
- [ ] Confirm `has_disorder == True` path is triggered (no silent fall-through to ordered path)
- [ ] Confirm `results.V_aa_samples.shape == (5000, 3)`
- [ ] Confirm `results.eta_samples.shape == (5000,)`
- [ ] Confirm `results.nu_Q_MHz_samples.shape == (5000,)`
- [ ] Confirm scalar `results.V_aa` == `np.mean(results.V_aa_samples, axis=0)` (to floating-point tolerance)
- [ ] Confirm scalar `results.eta` == `np.mean(results.eta_samples)`
- [ ] Confirm `results.nu_Q_MHz` == `np.mean(results.nu_Q_MHz_samples)`

### 4. Independent disorder -- physical sanity
- [ ] Mean `|nu_Q|` for octahedral Al probe should be in the range **1-10 MHz**
  (point-charge overestimates EFG; exact match to DFT not expected)
- [ ] Distribution of `|nu_Q|` should be **broader** than the ordered-spinel single value
  (disorder broadening is the key physics we're testing)
- [ ] Distribution of eta should span 0-1 (breaking of cubic symmetry by random neighbors)
- [ ] Mean |nu_Q| for **tetrahedral** Al probe (Part 1, probe at 8a site) should be
  noticeably **smaller** than for octahedral probe (~1.6 MHz expected vs. ~3.5 MHz)

### 5. Reproducibility
- [ ] Run the calculation twice with `random_seed=42` -> `nu_Q_MHz_samples` arrays must be **identical**
- [ ] Run with `random_seed=None` twice -> arrays must **differ**
- [ ] Run with `random_seed=123` -> arrays must differ from `random_seed=42` arrays

### 6. Convergence with n_monte_carlo_samples
- [ ] Run with n = [100, 500, 1000, 5000, 10000] and the same `random_seed`
- [ ] Plot mean |nu_Q| vs. n -- should stabilize (converge) by ~1000-2000 samples
- [ ] Plot std |nu_Q| vs. n -- should also stabilize

### 7. Correlated disorder groups
- [ ] Run Part 2 of the example script
- [ ] Confirm no ValueError is raised for valid `disorder_groups` input
- [ ] Confirm that paired sites (8a-Al / 16d-Mg or 8a-Mg / 16d-Al) are sampled together:
  add a debug print inside the MC loop to verify group sampling is applied
- [ ] Compare the `nu_Q` distribution from Part 2 (correlated) to Part 1 (independent):
  - The means should be **similar** (same average occupancy)
  - The shapes may differ -- correlated disorder enforces charge balance within each pair,
    which can narrow or shift the distribution relative to fully independent sampling

### 8. Input validation
Test that `ValueError` is raised for each of the following:

- [ ] A site index in `disorder_groups` that is NOT in `charge_site_indices`
- [ ] The same site index appearing in **two different** disorder groups
- [ ] `species` list length != `site_indices` length in a configuration
- [ ] A species string in a configuration that is NOT a key in `first_cell_charges[site_index]`
- [ ] Group probabilities summing to a value other than 1.0 (e.g. 0.7 + 0.7 = 1.4)

### 9. Edge cases
- [ ] `disorder_groups = []` with independently disordered sites -> works (no groups path)
- [ ] All sites in `charge_site_indices` are in groups -> works (no independent path)
- [ ] A group with only one configuration and probability 1.0 -> behaves like a fixed charge
- [ ] A site with total CIF occupancy exactly 1.0 but `isinstance(charge_entry, dict)` ->
  no vacancy should be added (check `total_occ < 1.0 - 1e-6` threshold)
- [ ] A site with total occupancy < 1.0 (e.g. 0.75 Mg + 0.0 Al) -> vacancy probability
  should be 0.25 and charge = 0.0 for vacancy draws
- [ ] `n_monte_carlo_samples = 1` -> still returns arrays of shape (1, 3), (1,), (1,)

### 10. Comparison with DFT/literature (longer-term)
- [ ] Compare mean point-charge nu_Q with published DFT-computed values:
  - For the **ordered** end-member spinel (i = 0): run with `disorder_groups=[]` and
    purely scalar charges, compare to DFT (e.g. Vosegaard et al. 2001 or Mauri GIPAW papers)
  - Note: point-charge will likely overestimate Cq by 2-5x compared to DFT/experiment;
    a Sternheimer anti-shielding factor gamma_inf ~ -11 (for Al) can partly correct this.
    Try `gamma_sternheimer = -11.0` and compare to literature.
- [ ] For the **disordered** case at 673 K (i = 0.155): confirm the distribution width
  increases relative to the ordered case, consistent with Kashii et al. (1999) observations.

---

## Files

| File | Purpose |
|------|---------|
| `MgAl2O4_spinel/Redfern_1999_MgAl2O4_673K_COD9002084.cif` | Test CIF with cation disorder (COD 9002084) |
| `MgAl2O4_spinel/MgAl2O4_spinel_disorder_EFG_example.py` | Full worked example (independent + correlated) |
| `pyEFGPointCharge.py` | Main module (all disorder logic lives here) |
| `aAl2O3_EFG_point_charge_calc_example.py` | Ordered reference example (backwards compat. check) |
