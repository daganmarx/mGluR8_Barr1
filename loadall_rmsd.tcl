# loadall_rmsd.tcl
#
# Description:
#   Batch analysis script for VMD. Iterates over all trajectory replicas and
#   computes per-frame RMSD for four structural regions of the mGluR dimer,
#   each aligned to its own frame-0 reference before measuring.
#
#   For each measurement, the full system is first superposed onto the
#   alignment selection's frame-0 reference, then the RMSD of a (potentially
#   different) selection is measured. This cross-chain RMSD design means:
#     - Aligning on chain A and measuring chain B RMSD captures asymmetric
#       drift between subunits relative to a common reference frame.
#
#   Metrics computed (all using backbone heavy atoms, no hydrogens):
#     1. Chain B full backbone RMSD, aligned to chain A frame 0
#        Written to: allchainA_RMSD_rep<k>.dat
#
#     2. Chain B full backbone RMSD, aligned to chain B frame 0
#        Written to: allchainB_RMSD_rep<k>.dat
#
#     3. Chain B structured-region RMSD (helix+sheet), aligned to chain A frame 0
#        Written to: struc_chainA_RMSD_rep<k>.dat
#
#     4. Chain B structured-region RMSD (helix+sheet), aligned to chain B frame 0
#        Written to: struc_chainB_RMSD_rep<k>.dat
#
# Requirements:
#   - VMD with Orient and tempoUserVMD packages installed
#   - PSF file: now_0000.psf
#   - DCD trajectory files: trajectory_now_<k>.dcd (k = 0 to 71)
#
# Segment IDs used:
#   PROA : mGluR chain A — used as alignment reference
#   PROD : mGluR chain B — used as RMSD target
#
# Note:
#   Arrestin RMSD is implemented but commented out. Uncomment the Barr/outfile5
#   blocks to enable it.
#
# Usage:
#   vmd -dispdev none -e loadall_rmsd.tcl

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

    # --- Alignment and RMSD selections ---
    # Full backbone (no H) for each chain — used for alignment and RMSD
    set LBDa [atomselect top "segid PROA and backbone and noh"]
    set LBDb [atomselect top "segid PROD and backbone and noh"]

    # Secondary-structure backbone only (helix + sheet, no H) for each chain
    set TMDa [atomselect top "segid PROA and (helix or sheet) and backbone and noh"]
    set TMDb [atomselect top "segid PROD and (helix or sheet) and backbone and noh"]

    # Uncomment to add arrestin RMSD:
    #set Barr [atomselect top "segid PROF and backbone and noh"]

    # Open output files for this replica
    set outfile1 [open allchainA_RMSD_rep$k.dat w]   ;# Full LBD: aligned to A, measure B
    set outfile2 [open allchainB_RMSD_rep$k.dat w]   ;# Full LBD: aligned to B, measure B
    set outfile3 [open struc_chainA_RMSD_rep$k.dat w] ;# Structured: aligned to A, measure B
    set outfile4 [open struc_chainB_RMSD_rep$k.dat w] ;# Structured: aligned to B, measure B
    #set outfile5 [open allBarr_RMSD_rep$k.dat w]

    set nf [molinfo top get numframes]  ;# Total number of loaded frames

    # Frame-0 reference selections (frozen at frame 0 for RMSD baseline)
    set frame0_LBDa [atomselect top "segid PROA and backbone and noh" frame 0]
    set frame0_LBDb [atomselect top "segid PROD and backbone and noh" frame 0]
    set frame0_TMDa [atomselect top "segid PROA and (helix or sheet) and backbone and noh" frame 0]
    set frame0_TMDb [atomselect top "segid PROD and (helix or sheet) and backbone and noh" frame 0]
    #set frame0_Barr [atomselect top "segid PROF and backbone and noh" frame 0]

    # --- Per-frame RMSD calculation ---
    for { set i 1 } { $i <= $nf } { incr i } {

        # Selection of all atoms at current frame (used for system-wide alignment moves)
        set all [atomselect top all frame $i]

        # --- Metric 1: Full backbone, align chain A -> measure chain B RMSD ---
        $LBDa frame $i
        $all move [measure fit $LBDa $frame0_LBDa]   ;# Fit whole system to chain A ref
        puts $outfile1 "[measure rmsd $LBDb $frame0_LBDb]"

        # --- Metric 2: Full backbone, align chain B -> measure chain B RMSD ---
        $LBDb frame $i
        $all move [measure fit $LBDb $frame0_LBDb]   ;# Fit whole system to chain B ref
        puts $outfile2 "[measure rmsd $LBDb $frame0_LBDb]"

        # --- Metric 3: Structured regions, align chain A -> measure chain B RMSD ---
        $TMDa frame $i
        $all move [measure fit $TMDa $frame0_TMDa]
        puts $outfile3 "[measure rmsd $TMDb $frame0_TMDb]"

        # --- Metric 4: Structured regions, align chain B -> measure chain B RMSD ---
        $TMDb frame $i
        $all move [measure fit $TMDb $frame0_TMDb]
        puts $outfile4 "[measure rmsd $TMDb $frame0_TMDb]"

        # Uncomment for arrestin RMSD:
        #$Barr frame $i
        #$all move [measure fit $Barr $frame0_Barr]
        #puts $outfile5 "[measure rmsd $Barr $frame0_Barr]"
    }

    # Close output files for this replica
    close $outfile1
    close $outfile2
    close $outfile3
    close $outfile4
    #close $outfile5

    # Clear trajectory frames before loading the next replica
    animate delete all
}

exit
