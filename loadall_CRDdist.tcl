# loadall_CRDdist.tcl
#
# Description:
#   Batch analysis script for VMD. Iterates over all trajectory replicas and
#   computes Ca-Ca distances between symmetry-equivalent residues in the
#   cysteine-rich domain (CRD) of mGluR chain A (PROA) and chain B (PROD).
#
#   For each of 9 CRD residue positions, the Euclidean distance between the
#   alpha-carbon of that residue in chain A and the corresponding residue in
#   chain B is computed per frame. This serves as a measure of inter-subunit
#   separation or asymmetry at the CRD interface.
#
#   Output (per replica, per residue pair):
#       mGluR8_CRD_<resid>_dist_rep<k>.dat — inter-subunit CA distance (Å)
#
# Requirements:
#   - VMD with Orient and tempoUserVMD packages installed
#   - PSF file: now_0000.psf
#   - DCD trajectory files: trajectory_now_<k>.dcd (k = 0 to 71)
#
# Segment IDs used:
#   PROA : mGluR chain A — CRD residues measured: 516, 520, 534, 535, 583, 541, 553, 556, 569
#   PROD : mGluR chain B — same residue positions as PROA
#
# Note:
#   atomsel5 selects resid 583 (not 538 as suggested by the output filename).
#   Verify residue numbering if adding new positions.
#   Additional residue pairs (803, 804, 792, 793) are implemented but
#   commented out — uncomment to enable.
#
# Usage:
#   vmd -dispdev none -e loadall_CRDdist.tcl

# Load the structural topology (PSF) file
mol load psf now_0000.psf

# Load required VMD packages
package require Orient        ;# Principal axis utilities
package require tempoUserVMD  ;# Provides dopbc for PBC correction
namespace import Orient::orient

# --- Main loop over trajectory replicas ---
for { set k 0 } { $k <= 71 } { incr k 1 } {

    # Load trajectory with PBC correction, wrapping fragments relative to chain A
    dopbc -file trajectory_now_$k.dcd -frames 0:1:1000 \
          -ref "segid PROA" -wrapby "fragment"

    set nf [molinfo top get numframes]  ;# Total number of loaded frames

    # Open one output file per CRD residue pair
    set outfile1 [open mGluR8_CRD_516_dist_rep$k.dat w]
    set outfile2 [open mGluR8_CRD_520_dist_rep$k.dat w]
    set outfile3 [open mGluR8_CRD_534_dist_rep$k.dat w]
    set outfile4 [open mGluR8_CRD_535_dist_rep$k.dat w]
    set outfile5 [open mGluR8_CRD_538_dist_rep$k.dat w]  ;# Note: selects resid 583
    set outfile6 [open mGluR8_CRD_541_dist_rep$k.dat w]
    set outfile7 [open mGluR8_CRD_553_dist_rep$k.dat w]
    set outfile8 [open mGluR8_CRD_556_dist_rep$k.dat w]
    set outfile9 [open mGluR8_CRD_569_dist_rep$k.dat w]
    # Uncomment to add additional residue pairs:
    #set outfile10 [open mGluR8_CRD_803_dist_rep$k.dat w]
    #set outfile11 [open mGluR8_CRD_804_dist_rep$k.dat w]

    # Alpha-carbon selections for chain A CRD residues
    set atomsel1  [atomselect top "segid PROA and resid 516 and name CA"]
    set atomsel2  [atomselect top "segid PROA and resid 520 and name CA"]
    set atomsel3  [atomselect top "segid PROA and resid 534 and name CA"]
    set atomsel4  [atomselect top "segid PROA and resid 535 and name CA"]
    set atomsel5  [atomselect top "segid PROA and resid 583 and name CA"]  ;# Note: resid 583
    set atomsel6  [atomselect top "segid PROA and resid 541 and name CA"]
    set atomsel7  [atomselect top "segid PROA and resid 553 and name CA"]
    set atomsel8  [atomselect top "segid PROA and resid 556 and name CA"]
    set atomsel9  [atomselect top "segid PROA and resid 569 and name CA"]
    #set atomsel10 [atomselect top "segid PROA and resid 803 and name CA"]
    #set atomsel11 [atomselect top "segid PROA and resid 804 and name CA"]

    # Alpha-carbon selections for the corresponding chain B CRD residues
    set atomsel14 [atomselect top "segid PROD and resid 516 and name CA"]
    set atomsel15 [atomselect top "segid PROD and resid 520 and name CA"]
    set atomsel16 [atomselect top "segid PROD and resid 534 and name CA"]
    set atomsel17 [atomselect top "segid PROD and resid 535 and name CA"]
    set atomsel18 [atomselect top "segid PROD and resid 583 and name CA"]
    set atomsel19 [atomselect top "segid PROD and resid 541 and name CA"]
    set atomsel20 [atomselect top "segid PROD and resid 553 and name CA"]
    set atomsel21 [atomselect top "segid PROD and resid 556 and name CA"]
    set atomsel22 [atomselect top "segid PROD and resid 569 and name CA"]
    #set atomsel23 [atomselect top "segid PROD and resid 803 and name CA"]
    #set atomsel24 [atomselect top "segid PROD and resid 804 and name CA"]

    # --- Per-frame analysis ---
    for {set i 0} {$i < $nf} {incr i} {

        # Update all CA selections to the current frame
        $atomsel1  frame $i
        $atomsel2  frame $i
        $atomsel3  frame $i
        $atomsel4  frame $i
        $atomsel5  frame $i
        $atomsel6  frame $i
        $atomsel7  frame $i
        $atomsel8  frame $i
        $atomsel9  frame $i
        $atomsel14 frame $i
        $atomsel15 frame $i
        $atomsel16 frame $i
        $atomsel17 frame $i
        $atomsel18 frame $i
        $atomsel19 frame $i
        $atomsel20 frame $i
        $atomsel21 frame $i
        $atomsel22 frame $i

        # Get XYZ coordinates for each CA atom
        set coord1  [lindex [$atomsel1  get {x y z}] 0]
        set coord2  [lindex [$atomsel2  get {x y z}] 0]
        set coord3  [lindex [$atomsel3  get {x y z}] 0]
        set coord4  [lindex [$atomsel4  get {x y z}] 0]
        set coord5  [lindex [$atomsel5  get {x y z}] 0]
        set coord6  [lindex [$atomsel6  get {x y z}] 0]
        set coord7  [lindex [$atomsel7  get {x y z}] 0]
        set coord8  [lindex [$atomsel8  get {x y z}] 0]
        set coord9  [lindex [$atomsel9  get {x y z}] 0]
        set coord14 [lindex [$atomsel14 get {x y z}] 0]
        set coord15 [lindex [$atomsel15 get {x y z}] 0]
        set coord16 [lindex [$atomsel16 get {x y z}] 0]
        set coord17 [lindex [$atomsel17 get {x y z}] 0]
        set coord18 [lindex [$atomsel18 get {x y z}] 0]
        set coord19 [lindex [$atomsel19 get {x y z}] 0]
        set coord20 [lindex [$atomsel20 get {x y z}] 0]
        set coord21 [lindex [$atomsel21 get {x y z}] 0]
        set coord22 [lindex [$atomsel22 get {x y z}] 0]

        # Compute inter-subunit difference vectors (chain A - chain B)
        set v1 [vecsub $coord1  $coord14]
        set v2 [vecsub $coord2  $coord15]
        set v3 [vecsub $coord3  $coord16]
        set v4 [vecsub $coord4  $coord17]
        set v5 [vecsub $coord5  $coord18]
        set v6 [vecsub $coord6  $coord19]
        set v7 [vecsub $coord7  $coord20]
        set v8 [vecsub $coord8  $coord21]
        set v9 [vecsub $coord9  $coord22]

        # Compute scalar inter-subunit CA-CA distances
        set length_vector1 [veclength $v1]
        set length_vector2 [veclength $v2]
        set length_vector3 [veclength $v3]
        set length_vector4 [veclength $v4]
        set length_vector5 [veclength $v5]
        set length_vector6 [veclength $v6]
        set length_vector7 [veclength $v7]
        set length_vector8 [veclength $v8]
        set length_vector9 [veclength $v9]

        # Write distances to output files
        puts $outfile1 "$length_vector1"
        puts $outfile2 "$length_vector2"
        puts $outfile3 "$length_vector3"
        puts $outfile4 "$length_vector4"
        puts $outfile5 "$length_vector5"
        puts $outfile6 "$length_vector6"
        puts $outfile7 "$length_vector7"
        puts $outfile8 "$length_vector8"
        puts $outfile9 "$length_vector9"
    }

    # Close output files for this replica
    close $outfile1
    close $outfile2
    close $outfile3
    close $outfile4
    close $outfile5
    close $outfile6
    close $outfile7
    close $outfile8
    close $outfile9
    #close $outfile10
    #close $outfile11

    # Clear trajectory frames before loading the next replica
    animate delete all
}

exit
