# mGluR8_Barr1
Code and Data for mGluR8_Barr1 paper
.dcd and .psf files can be accessed in Zenodo at the following urls:
Simulations 0-35: https://zenodo.org/records/19898710
Simulations 36-71: https://zenodo.org/records/20029951


# VMD MD Simulation Analysis Scripts

TCL scripts for analyzing molecular dynamics (MD) simulations of an arrestin-1 (Barr1) / mGluR complex in [VMD](https://www.ks.uiuc.edu/Research/vmd/). Each script loads one or more DCD trajectory replicas and extracts structural metrics frame-by-frame, writing results to plain-text `.dat` / `.txt` files for downstream analysis.

---

## System Description

These scripts assume a simulation system containing:

| Segment ID | Molecule | Residues used |
|---|---|---|
| `PROF` | Arrestin-1 (Barr1) | 1–175 (N-lobe), 176–360 (C-lobe) |
| `PROA` | mGluR chain A | 580–900 |

Topology is provided by a PSF file (`now_0000.psf`). Trajectories are DCD files named `trajectory_now_<k>.dcd`.

---

## Scripts

### `loadall_Barractangle.tcl`
Computes the **arrestin inter-lobe bending angle** across all trajectory replicas (k = 0–71).

The angle is defined as the angle between the first principal axes of the arrestin N-lobe and C-lobe. Angles greater than 100° are reflected to their supplement to handle principal axis sign ambiguity and keep values in a consistent 0–90° range.

**Output (per replica):**
- `Barr1_prinax1_angle_rep<k>.dat` — inter-lobe bending angle in degrees, one value per frame

---

### `loadall_xyz.tcl`
Computes two angular metrics describing the **orientation of arrestin relative to the mGluR receptor** across all trajectory replicas (k = 20–71).

1. **XY-plane angle** (`adj_ang`): angle between the arrestin inter-lobe vector (projected onto XY) and the second principal axis of mGluR chain A (projected onto XY). Adjusted so that 0° = aligned with the receptor axis.
2. **Z-tilt angle** (`ang_deg_z`): angle between the arrestin inter-lobe vector and the membrane-normal (Z-axis). Adjusted so that 0° = membrane-parallel.

**Output (per replica):**
- `Barr1TMD_newxy_angle_rep<k>.dat` — XY-plane angle (degrees)
- `Barr1chainA_z_angle_rep<k>.dat` — Z-tilt angle (degrees)

---

### `xy_pos.tcl`
Tracks the **XY position of a single residue's alpha-carbon (CA)** across all frames of a trajectory already loaded in VMD. The trajectory is aligned to a reference frame (frame 0 of mGluR chain A residues 580–900) before coordinate extraction.

Edit the `residue` variable at the top of the script to change the target residue (default: `resid 777`).

**Output:**
- `long_chainA_777_x.txt` — CA X coordinates, one value per frame
- `long_chainA_777_y.txt` — CA Y coordinates, one value per frame

---

## Requirements

- [VMD](https://www.ks.uiuc.edu/Research/vmd/) (tested with VMD 1.9.x)
- VMD packages:
  - `Orient` — for principal axis calculations
  - `tempoUserVMD` — for PBC correction via `dopbc`

---

## Usage

### Batch scripts (loadall_*.tcl)
These scripts load the PSF and trajectories themselves. Run headlessly from the command line:

```bash
vmd -dispdev none -e loadall_Barractangle.tcl
vmd -dispdev none -e loadall_xyz.tcl
```

Or source interactively within a VMD session:
```tcl
source loadall_Barractangle.tcl
```

### Single-trajectory script (xy_pos.tcl)
This script operates on a molecule already loaded in VMD (`mol top`). Load your trajectory first, then:
```tcl
source xy_pos.tcl
```

---

## Output Format

All output files are plain text with one numeric value per line, corresponding to one trajectory frame. They can be read directly into Python (e.g. `numpy.loadtxt`), R, or MATLAB for plotting and further analysis.

---

## Notes

- `loadall_Barractangle.tcl` contains commented-out blocks for computing angles along principal axes 2 and 3. Uncomment the relevant sections to enable them.
- PBC correction (`dopbc`) is applied in both batch scripts, wrapping fragments relative to `segid PROA`.
- `animate delete all` is called between replicas to free memory.
