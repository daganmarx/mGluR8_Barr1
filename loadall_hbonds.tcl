# loadall_hbonds.tcl
#
# Description:
#   Batch analysis script for VMD. Iterates over all trajectory replicas and
#   counts hydrogen bonds between specific arrestin residues (or the full
#   N-lobe) and individual lipid species in the membrane.
#
#   For each replica, the VMD `hbonds` plugin is called with standard geometric
#   criteria (default cutoffs: distance < 3.5 Å, angle < 30°). Results are
#   written directly to output files by the plugin (-writefile yes).
#
#   Residues analyzed (all from segid PROF, arrestin-1):
#     - Residue 48  (individual)
#     - Residue 49  (individual)
#     - Residue 50  (individual)
#     - N-lobe      (residues 1-175)
#
#   Lipid species analyzed:
#     - POPS  (phosphatidylserine)
#     - SAPI24 (phosphatidylinositol, 24:0 chain)
#     - SAPI25 (phosphatidylinositol, 25:0 chain)
#     - POPC  (phosphatidylcholine)
#     - POPE  (phosphatidylethanolamine)
#     - SSM   (sphingomyelin)
#
#   Output (per replica, per residue, per lipid):
#       <residue>_<lipid>_hbonds_rep<k>.dat
#   e.g.: 48_PS_hbonds_rep0.dat, Nlobe_SAPI24_hbonds_rep3.dat
#
# Requirements:
#   - VMD with Orient, tempoUserVMD, and hbonds plugins installed
#   - PSF file: now_0000.psf
#   - DCD trajectory files: trajectory_now_<k>.dcd (k = 0 to 71)
#
# Segment IDs used:
#   PROF : Arrestin-1 (Barr1) — residues 48, 49, 50, and N-lobe (1-175)
#
# Note:
#   PBC wrapping uses the full protein as the centering reference. If
#   membrane-centered wrapping is preferred, change -ref to "lipid".
#
# Usage:
#   vmd -dispdev none -e loadall_hbonds.tcl

# Load the structural topology (PSF) file
mol load psf now_0000.psf

# Load required VMD packages
package require Orient        ;# Principal axis utilities
package require tempoUserVMD  ;# Provides dopbc for PBC correction
namespace import Orient::orient

# --- Main loop over trajectory replicas ---
for { set k 0 } { $k <= 71 } { incr k 1 } {

    # Load trajectory with PBC correction, centering on the full protein
    dopbc -file trajectory_now_$k.dcd -frames 0:1:1000 \
          -ref "protein" -wrapby "fragment"

    # --- Hydrogen bond analysis: Arrestin residue 48 vs each lipid species ---
    # hbonds writes one value per frame to the output file automatically
    hbonds -sel1 [atomselect top "segid PROF and resid 48"] \
           -sel2 [atomselect top "resname POPS"]   \
           -plot no -writefile yes -outfile 48_PS_hbonds_rep$k.dat

    hbonds -sel1 [atomselect top "segid PROF and resid 48"] \
           -sel2 [atomselect top "resname SAPI24"] \
           -plot no -writefile yes -outfile 48_SAPI24_hbonds_rep$k.dat

    hbonds -sel1 [atomselect top "segid PROF and resid 48"] \
           -sel2 [atomselect top "resname SAPI25"] \
           -plot no -writefile yes -outfile 48_SAPI25_hbonds_rep$k.dat

    hbonds -sel1 [atomselect top "segid PROF and resid 48"] \
           -sel2 [atomselect top "resname POPC"]   \
           -plot no -writefile yes -outfile 48_PC_hbonds_rep$k.dat

    hbonds -sel1 [atomselect top "segid PROF and resid 48"] \
           -sel2 [atomselect top "resname POPE"]   \
           -plot no -writefile yes -outfile 48_PE_hbonds_rep$k.dat

    hbonds -sel1 [atomselect top "segid PROF and resid 48"] \
           -sel2 [atomselect top "resname SSM"]    \
           -plot no -writefile yes -outfile 48_SSM_hbonds_rep$k.dat

    # --- Hydrogen bond analysis: Arrestin residue 49 vs each lipid species ---
    hbonds -sel1 [atomselect top "segid PROF and resid 49"] \
           -sel2 [atomselect top "resname POPS"]   \
           -plot no -writefile yes -outfile 49_PS_hbonds_rep$k.dat

    hbonds -sel1 [atomselect top "segid PROF and resid 49"] \
           -sel2 [atomselect top "resname SAPI24"] \
           -plot no -writefile yes -outfile 49_SAPI24_hbonds_rep$k.dat

    hbonds -sel1 [atomselect top "segid PROF and resid 49"] \
           -sel2 [atomselect top "resname SAPI25"] \
           -plot no -writefile yes -outfile 49_SAPI25_hbonds_rep$k.dat

    hbonds -sel1 [atomselect top "segid PROF and resid 49"] \
           -sel2 [atomselect top "resname POPC"]   \
           -plot no -writefile yes -outfile 49_PC_hbonds_rep$k.dat

    hbonds -sel1 [atomselect top "segid PROF and resid 49"] \
           -sel2 [atomselect top "resname POPE"]   \
           -plot no -writefile yes -outfile 49_PE_hbonds_rep$k.dat

    hbonds -sel1 [atomselect top "segid PROF and resid 49"] \
           -sel2 [atomselect top "resname SSM"]    \
           -plot no -writefile yes -outfile 49_SSM_hbonds_rep$k.dat

    # --- Hydrogen bond analysis: Arrestin residue 50 vs each lipid species ---
    hbonds -sel1 [atomselect top "segid PROF and resid 50"] \
           -sel2 [atomselect top "resname POPS"]   \
           -plot no -writefile yes -outfile 50_PS_hbonds_rep$k.dat

    hbonds -sel1 [atomselect top "segid PROF and resid 50"] \
           -sel2 [atomselect top "resname SAPI24"] \
           -plot no -writefile yes -outfile 50_SAPI24_hbonds_rep$k.dat

    hbonds -sel1 [atomselect top "segid PROF and resid 50"] \
           -sel2 [atomselect top "resname SAPI25"] \
           -plot no -writefile yes -outfile 50_SAPI25_hbonds_rep$k.dat

    hbonds -sel1 [atomselect top "segid PROF and resid 50"] \
           -sel2 [atomselect top "resname POPC"]   \
           -plot no -writefile yes -outfile 50_PC_hbonds_rep$k.dat

    hbonds -sel1 [atomselect top "segid PROF and resid 50"] \
           -sel2 [atomselect top "resname POPE"]   \
           -plot no -writefile yes -outfile 50_PE_hbonds_rep$k.dat

    hbonds -sel1 [atomselect top "segid PROF and resid 50"] \
           -sel2 [atomselect top "resname SSM"]    \
           -plot no -writefile yes -outfile 50_SSM_hbonds_rep$k.dat

    # --- Hydrogen bond analysis: Arrestin N-lobe (res 1-175) vs each lipid species ---
    hbonds -sel1 [atomselect top "segid PROF and resid 1 to 175"] \
           -sel2 [atomselect top "resname POPS"]   \
           -plot no -writefile yes -outfile Nlobe_PS_hbonds_rep$k.dat

    hbonds -sel1 [atomselect top "segid PROF and resid 1 to 175"] \
           -sel2 [atomselect top "resname SAPI24"] \
           -plot no -writefile yes -outfile Nlobe_SAPI24_hbonds_rep$k.dat

    hbonds -sel1 [atomselect top "segid PROF and resid 1 to 175"] \
           -sel2 [atomselect top "resname SAPI25"] \
           -plot no -writefile yes -outfile Nlobe_SAPI25_hbonds_rep$k.dat

    hbonds -sel1 [atomselect top "segid PROF and resid 1 to 175"] \
           -sel2 [atomselect top "resname POPC"]   \
           -plot no -writefile yes -outfile Nlobe_PC_hbonds_rep$k.dat

    hbonds -sel1 [atomselect top "segid PROF and resid 1 to 175"] \
           -sel2 [atomselect top "resname POPE"]   \
           -plot no -writefile yes -outfile Nlobe_PE_hbonds_rep$k.dat

    hbonds -sel1 [atomselect top "segid PROF and resid 1 to 175"] \
           -sel2 [atomselect top "resname SSM"]    \
           -plot no -writefile yes -outfile Nlobe_SSM_hbonds_rep$k.dat

    # Clear trajectory frames before loading the next replica
    animate delete all
}

exit
