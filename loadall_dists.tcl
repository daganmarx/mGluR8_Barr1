# loadall_dists.tcl
#
# Description:
#   Batch analysis script for VMD. Iterates over all trajectory replicas and
#   computes inter-subunit Ca-Ca distances across a stretch of TM6 residues
#   (784-804) between mGluR chain A (PROA) and chain B (PROD).
#
#   For each frame, the distances for all residues in the range are written
#   sequentially to a single output file, producing a 2D data block with
#   shape (nframes * 21 residues) x 1. The residue index is implicit (order
#   is fixed: 784, 785, ..., 804).
#
#   Output (per replica):
#       TM6_dist_rep<k>.txt — CA-CA distances (Å), 21 values per frame
#
# Requirements:
#   - VMD with Orient and tempoUserVMD packages installed
#   - PSF file: now_0000.psf
#   - DCD trajectory files: trajectory_now_<k>.dcd (k = 0 to 71)
#
# Segment IDs used:
#   PROA : mGluR chain A — TM6 residues 784-804
#   PROD : mGluR chain B — TM6 residues 784-804
#
# Usage:
#   vmd -dispdev none -e loadall_dists.tcl

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

    # Output file: all TM6 inter-subunit distances for this replica
    set outfile1 [open TM6_dist_rep$k.txt w]

    # --- Per-frame, per-residue analysis ---
    for {set i 0} {$i < $nf} {incr i} {

        # Loop over TM6 residues 784-804
        for {set j 784} {$j <= 804} {incr j 1} {

            # Select the CA of residue j in each subunit at the current frame
            set atomsel1 [atomselect top "segid PROA and resid $j and name CA" frame $i]
            set atomsel2 [atomselect top "segid PROD and resid $j and name CA" frame $i]

            # Get XYZ coordinates
            set coord1 [lindex [$atomsel1 get {x y z}] 0]
            set coord2 [lindex [$atomsel2 get {x y z}] 0]

            # Compute inter-subunit distance vector and its magnitude
            set v1              [vecsub $coord1 $coord2]
            set length_vector1  [veclength $v1]

            # Write distance to output (one value per line, residues in order)
            puts $outfile1 $length_vector1
        }
    }

    # Close output file for this replica
    close $outfile1

    # Clear trajectory frames before loading the next replica
    animate delete all
}

exit
