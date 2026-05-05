# mGluR8_Barr1
Code and Data for mGluR8_Barr1 paper
.dcd and .psf files can be accessed in Zenodo at the following urls:
Simulations 0-35: https://zenodo.org/records/19898710
Simulations 36-71: https://zenodo.org/records/20029951


# MD Simulation Analysis Scripts

Analysis scripts for molecular dynamics simulations of an **arrestin-1 (Barr1) / mGluR dimer** complex. Includes VMD/TCL scripts for trajectory analysis and Python scripts for post-processing and visualization.

---

## System

| Segment | Molecule | Key regions |
|---|---|---|
| `PROA` | mGluR chain A | LBD (1–499), CRD (500–579), TMD (580–850) |
| `PROD` | mGluR chain B | LBD (1–499), CRD (500–579), TMD (580–850) |
| `PROF` | Arrestin-1 (Barr1) | N-lobe (1–175), C-lobe (176–360), finger loop (62–74) |

Topology: `now_0000.psf` | Trajectories: `trajectory_now_<k>.dcd` (k = 0–71)

---

## Requirements

**VMD scripts**
- [VMD](https://www.ks.uiuc.edu/Research/vmd/) 1.9.x
- `Orient` and `tempoUserVMD` packages
- `hbonds` plugin (`loadall_hbonds.tcl` only)

**Python scripts**
- Python 3.8+
- `MDAnalysis`, `numpy`, `pandas`, `matplotlib`, `openpyxl`

```bash
pip install MDAnalysis numpy pandas matplotlib openpyxl
```

---

## VMD/TCL Scripts

All batch scripts load the PSF and trajectories internally. Run headlessly:
```bash
vmd -dispdev none -e <script>.tcl
```

| Script | Description | Output |
|---|---|---|
| `loadall_Barractangle.tcl` | Arrestin inter-lobe bending angle (angle between first principal axes of N-lobe and C-lobe) | `Barr1_prinax1_angle_rep<k>.dat` |
| `loadall_xyz.tcl` | Arrestin orientation relative to mGluR: XY-plane angle and Z-tilt angle (replicas 20–71) | `Barr1TMD_newxy_angle_rep<k>.dat`, `Barr1chainA_z_angle_rep<k>.dat` |
| `loadall_COMdist.tcl` | Center-of-mass distances: CRD inter-subunit, TMD inter-subunit, finger loop to TMD-A | `CRD_com_rep<k>.txt`, `TMD_com_rep<k>.txt`, `FL_com_rep<k>.txt` |
| `loadall_CRDdist.tcl` | Inter-subunit CA–CA distances at 9 CRD residue positions (chain A vs B) | `mGluR8_CRD_<resid>_dist_rep<k>.dat` |
| `loadall_ECDdist.tcl` | ECD/LBD distances: lobe closure (res 155–284), inter-subunit distances at res 229 and 553, lower-lobe COM distance | `LBDa/b_closure_dist_rep<k>.txt`, `LBD_229_dist_rep<k>.txt`, `CRD_553_dist_rep<k>.txt`, `LBDcom_dist_rep<k>.txt` |
| `loadall_TM6dist.tcl` | Inter-subunit CA–CA distances at TM6 residues 792 and 793 (cytoplasmic kink region) | `mGluR8_TM6_792_dist_rep<k>.dat`, `mGluR8_TM6_793_dist_rep<k>.dat` |
| `loadall_dists.tcl` | Inter-subunit CA–CA distances across TM6 residues 784–804 (21 values per frame, sequential) | `TM6_dist_rep<k>.txt` |
| `loadall_hbonds.tcl` | Hydrogen bonds between arrestin residues 48–50 and N-lobe (res 1–175) and six lipid species (POPS, SAPI24, SAPI25, POPC, POPE, SSM) | `<residue>_<lipid>_hbonds_rep<k>.dat` |
| `loadall_rmsd.tcl` | Backbone RMSD (heavy atoms) for four alignment/measurement combinations across both mGluR subunits | `allchainA/B_RMSD_rep<k>.dat`, `struc_chainA/B_RMSD_rep<k>.dat` |
| `loadall_sasa.tcl` | Buried SASA of arrestin N-lobe loops (res 238–248, 305–317) upon contact with mGluR TMD chain A (currently replica 39, every other frame) | `Nloops_chainA_rep<k>.txt` |

All output files are plain text, one value per line (one per frame), readable with `numpy.loadtxt`.

---

## Python Scripts

### `rmsd_to_prism.py`

Computes Cα RMSD over time for each simulation replica and outputs a wide-format Excel table ready for import into GraphPad Prism.

```
Time (ns) | sim_0 | sim_1 | sim_2 | ...
```

```bash
# Basic usage
python rmsd_to_prism.py --input /path/to/folder/

# With chain and residue selection
python rmsd_to_prism.py --input /path/to/folder/ --chains PROA --resid 50:200

# Multiple chains, multiple residue ranges
python rmsd_to_prism.py --input /path/to/folder/ \
    --chains PROA PROD --resid 50:100,150:300 --output rmsd_TM.xlsx
```

**Options**

| Flag | Default | Description |
|---|---|---|
| `--input` | *(required)* | Folder containing PSF + DCD files |
| `--psf` | `now_0000.psf` | PSF filename |
| `--dcd-pattern` | `trajectory_now_{:d}.dcd` | DCD filename pattern |
| `--n-sims` | `72` | Number of simulations |
| `--chains` | all protein | One or more segIDs (e.g. `PROA PROD`) |
| `--resid` | all | Residue range(s), e.g. `50:200` or `50:100,150:200` |
| `--output` | `rmsd_prism.xlsx` | Output Excel filename |

---

### `membrane_area_convergence.py`

Tracks membrane XY box area (Å²) over time across all replicas and assesses convergence using coefficient of variation (CV%) in the production window. Outputs a formatted Excel workbook and an overlay plot.

```bash
# Basic usage
python membrane_area_convergence.py --input /path/to/folder/

# Custom thresholds
python membrane_area_convergence.py --input /path/to/folder/ \
    --area-threshold 3.0 --eq-fraction 0.33 --output membrane.xlsx
```

**Options**

| Flag | Default | Description |
|---|---|---|
| `--input` | *(required)* | Folder containing PSF + DCD files |
| `--psf` | `now_0000.psf` | PSF filename |
| `--dcd-pattern` | `trajectory_now_{:d}.dcd` | DCD filename pattern |
| `--n-sims` | `72` | Number of simulations |
| `--area-threshold` | `5.0` | CV (%) above which a simulation is flagged FAIL |
| `--eq-fraction` | `0.5` | Fraction of trajectory treated as production window |
| `--output` | `membrane_area_convergence.xlsx` | Output Excel filename |

**Output**
- Excel workbook: summary convergence matrix (PASS/WARN/FAIL per replica) + per-simulation time series sheets
- PNG overlay plot: all replica membrane area traces colored by convergence status

