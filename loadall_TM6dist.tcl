# loadall_TM6dist.tcl
#
# Description:
#   Batch analysis script for VMD. Iterates over all trajectory replicas and
#   computes inter-subunit Ca-Ca distances for TM6 residues 792 and 793 of
#   the mGluR dimer (chain A vs chain B). These residues are located near the
#   cytoplasmic end of TM6 and report on the opening of the TM6 kink, which
#   is associated with receptor activation and G protein/arrestin coupling.
#
#   The script defines CA selections for residues 781-793 on both chains,
#   computes inter-subunit distances for all of them, but currently writes
#   output only for residues 792 and 793. Distances for residues 781-791
#   are computed but commented out.
#
#   Output (per replica):
#       mGluR8_TM6_792_dist_rep<k>.dat — inter-subunit CA-CA distance at res 792 (Å)
#       mGluR8_TM6_793_dist_rep<k>.dat — inter-subunit CA-CA distance at res 793 (Å)
#
# Requirements:
#   - VMD with Orient and tempoUserVMD packages installed
#   - PSF file: now_0000.psf
#   - DCD trajectory files: trajectory_now_<k>.dcd (k = 0 to 71)
#
# Segment IDs used:
#   PROA : mGluR chain A — TM6 residues 781-793
#   PROD : mGluR chain B — TM6 residues 781-793
#
# Note:
#   To enable output for residues 781-791, uncomment the corresponding
#   outfile declarations, puts statements, and close statements.
#
# Usage:
#   vmd -dispdev none -e loadall_TM6dist.tcl

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

    # Open output files for residues 792 and 793 only
    # Uncomment additional outfile lines to enable other residues (781-791)
    #set outfile1  [open mGluR8_TM6_781_dist_rep$k.dat w]
    #set outfile2  [open mGluR8_TM6_782_dist_rep$k.dat w]
    #set outfile3  [open mGluR8_TM6_783_dist_rep$k.dat w]
    #set outfile4  [open mGluR8_TM6_784_dist_rep$k.dat w]
    #set outfile5  [open mGluR8_TM6_785_dist_rep$k.dat w]
    #set outfile6  [open mGluR8_TM6_786_dist_rep$k.dat w]
    #set outfile7  [open mGluR8_TM6_787_dist_rep$k.dat w]
    #set outfile8  [open mGluR8_TM6_788_dist_rep$k.dat w]
    #set outfile9  [open mGluR8_TM6_789_dist_rep$k.dat w]
    #set outfile10 [open mGluR8_TM6_790_dist_rep$k.dat w]
    #set outfile11 [open mGluR8_TM6_791_dist_rep$k.dat w]
    set outfile12 [open mGluR8_TM6_792_dist_rep$k.dat w]
    set outfile13 [open mGluR8_TM6_793_dist_rep$k.dat w]

    # CA selections for TM6 residues 781-793 on chain A (PROA)
    set atomsel1  [atomselect top "segid PROA and resid 781 and name CA"]
    set atomsel2  [atomselect top "segid PROA and resid 782 and name CA"]
    set atomsel3  [atomselect top "segid PROA and resid 783 and name CA"]
    set atomsel4  [atomselect top "segid PROA and resid 784 and name CA"]
    set atomsel5  [atomselect top "segid PROA and resid 785 and name CA"]
    set atomsel6  [atomselect top "segid PROA and resid 786 and name CA"]
    set atomsel7  [atomselect top "segid PROA and resid 787 and name CA"]
    set atomsel8  [atomselect top "segid PROA and resid 788 and name CA"]
    set atomsel9  [atomselect top "segid PROA and resid 789 and name CA"]
    set atomsel10 [atomselect top "segid PROA and resid 790 and name CA"]
    set atomsel11 [atomselect top "segid PROA and resid 791 and name CA"]
    set atomsel12 [atomselect top "segid PROA and resid 792 and name CA"]
    set atomsel13 [atomselect top "segid PROA and resid 793 and name CA"]

    # CA selections for the corresponding residues on chain B (PROD)
    set atomsel14 [atomselect top "segid PROD and resid 781 and name CA"]
    set atomsel15 [atomselect top "segid PROD and resid 782 and name CA"]
    set atomsel16 [atomselect top "segid PROD and resid 783 and name CA"]
    set atomsel17 [atomselect top "segid PROD and resid 784 and name CA"]
    set atomsel18 [atomselect top "segid PROD and resid 785 and name CA"]
    set atomsel19 [atomselect top "segid PROD and resid 786 and name CA"]
    set atomsel20 [atomselect top "segid PROD and resid 787 and name CA"]
    set atomsel21 [atomselect top "segid PROD and resid 788 and name CA"]
    set atomsel22 [atomselect top "segid PROD and resid 789 and name CA"]
    set atomsel23 [atomselect top "segid PROD and resid 790 and name CA"]
    set atomsel24 [atomselect top "segid PROD and resid 791 and name CA"]
    set atomsel25 [atomselect top "segid PROD and resid 792 and name CA"]
    set atomsel26 [atomselect top "segid PROD and resid 793 and name CA"]

    # --- Per-frame analysis ---
    for {set i 0} {$i < $nf} {incr i} {

        # Update all CA selections to the current frame
        $atomsel1  frame $i ; $atomsel2  frame $i ; $atomsel3  frame $i
        $atomsel4  frame $i ; $atomsel5  frame $i ; $atomsel6  frame $i
        $atomsel7  frame $i ; $atomsel8  frame $i ; $atomsel9  frame $i
        $atomsel10 frame $i ; $atomsel11 frame $i ; $atomsel12 frame $i
        $atomsel13 frame $i ; $atomsel14 frame $i ; $atomsel15 frame $i
        $atomsel16 frame $i ; $atomsel17 frame $i ; $atomsel18 frame $i
        $atomsel19 frame $i ; $atomsel20 frame $i ; $atomsel21 frame $i
        $atomsel22 frame $i ; $atomsel23 frame $i ; $atomsel24 frame $i
        $atomsel25 frame $i ; $atomsel26 frame $i

        # Get XYZ coordinates for all CA atoms
        set coord1  [lindex [$atomsel1  get {x y z}] 0]
        set coord2  [lindex [$atomsel2  get {x y z}] 0]
        set coord3  [lindex [$atomsel3  get {x y z}] 0]
        set coord4  [lindex [$atomsel4  get {x y z}] 0]
        set coord5  [lindex [$atomsel5  get {x y z}] 0]
        set coord6  [lindex [$atomsel6  get {x y z}] 0]
        set coord7  [lindex [$atomsel7  get {x y z}] 0]
        set coord8  [lindex [$atomsel8  get {x y z}] 0]
        set coord9  [lindex [$atomsel9  get {x y z}] 0]
        set coord10 [lindex [$atomsel10 get {x y z}] 0]
        set coord11 [lindex [$atomsel11 get {x y z}] 0]
        set coord12 [lindex [$atomsel12 get {x y z}] 0]
        set coord13 [lindex [$atomsel13 get {x y z}] 0]
        set coord14 [lindex [$atomsel14 get {x y z}] 0]
        set coord15 [lindex [$atomsel15 get {x y z}] 0]
        set coord16 [lindex [$atomsel16 get {x y z}] 0]
        set coord17 [lindex [$atomsel17 get {x y z}] 0]
        set coord18 [lindex [$atomsel18 get {x y z}] 0]
        set coord19 [lindex [$atomsel19 get {x y z}] 0]
        set coord20 [lindex [$atomsel20 get {x y z}] 0]
        set coord21 [lindex [$atomsel21 get {x y z}] 0]
        set coord22 [lindex [$atomsel22 get {x y z}] 0]
        set coord23 [lindex [$atomsel23 get {x y z}] 0]
        set coord24 [lindex [$atomsel24 get {x y z}] 0]
        set coord25 [lindex [$atomsel25 get {x y z}] 0]
        set coord26 [lindex [$atomsel26 get {x y z}] 0]

        # Compute inter-subunit difference vectors (chain A - chain B)
        set v1  [vecsub $coord1  $coord14]
        set v2  [vecsub $coord2  $coord15]
        set v3  [vecsub $coord3  $coord16]
        set v4  [vecsub $coord4  $coord17]
        set v5  [vecsub $coord5  $coord18]
        set v6  [vecsub $coord6  $coord19]
        set v7  [vecsub $coord7  $coord20]
        set v8  [vecsub $coord8  $coord21]
        set v9  [vecsub $coord9  $coord22]
        set v10 [vecsub $coord10 $coord23]
        set v11 [vecsub $coord11 $coord24]
        set v12 [vecsub $coord12 $coord25]
        set v13 [vecsub $coord13 $coord26]

        # Compute scalar inter-subunit CA-CA distances
        set length_vector1  [veclength $v1]
        set length_vector2  [veclength $v2]
        set length_vector3  [veclength $v3]
        set length_vector4  [veclength $v4]
        set length_vector5  [veclength $v5]
        set length_vector6  [veclength $v6]
        set length_vector7  [veclength $v7]
        set length_vector8  [veclength $v8]
        set length_vector9  [veclength $v9]
        set length_vector10 [veclength $v10]
        set length_vector11 [veclength $v11]
        set length_vector12 [veclength $v12]
        set length_vector13 [veclength $v13]

        # Write output for active residues 792 and 793 only
        # Uncomment the appropriate puts lines to enable additional residues
        #puts $outfile1  "$length_vector1"
        #puts $outfile2  "$length_vector2"
        #puts $outfile3  "$length_vector3"
        #puts $outfile4  "$length_vector4"
        #puts $outfile5  "$length_vector5"
        #puts $outfile6  "$length_vector6"
        #puts $outfile7  "$length_vector7"
        #puts $outfile8  "$length_vector8"
        #puts $outfile9  "$length_vector9"
        #puts $outfile10 "$length_vector10"
        #puts $outfile11 "$length_vector11"
        puts $outfile12 "$length_vector12"
        puts $outfile13 "$length_vector13"
    }

    # Close output files for this replica
    close $outfile12
    close $outfile13

    # Clear trajectory frames before loading the next replica
    animate delete all
}

exit
