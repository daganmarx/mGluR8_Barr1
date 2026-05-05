# loadall_ECDdist.tcl
#
# Description:
#   Batch analysis script for VMD. Iterates over all trajectory replicas and
#   computes five distance metrics per frame related to the extracellular domain
#   (ECD/LBD) and cysteine-rich domain (CRD) of the mGluR dimer.
#
#   Metrics computed:
#     1. LBD closure distance, chain A: CA distance between res 155 and res 284
#        (reports on the open/closed state of the ligand-binding domain lobe)
#        Written to: LBDa_closure_dist_A155_E284_rep<k>.txt
#
#     2. LBD closure distance, chain B: same metric for chain B (PROD)
#        Written to: LBDb_closure_dist_A155_E284_rep<k>.txt
#
#     3. Inter-subunit LBD distance at res 229: CA-CA distance between the
#        equivalent residue on chain A and chain B
#        Written to: LBD_229_dist_rep<k>.txt
#
#     4. Inter-subunit CRD distance at res 553: CA-CA distance between
#        chain A and chain B
#        Written to: CRD_553_dist_rep<k>.txt
#
#     5. Inter-subunit lower lobe COM distance: distance between the mass-
#        weighted COMs of the lower lobe of the LBD (residues 195-333 and
#        469-511, CA only) in chain A vs chain B
#        Written to: LBDcom_dist_rep<k>.txt
#
# Requirements:
#   - VMD with Orient and tempoUserVMD packages installed
#   - PSF file: now_0000.psf
#   - DCD trajectory files: trajectory_now_<k>.dcd (k = 0 to 71)
#
# Segment IDs used:
#   PROA : mGluR chain A
#   PROD : mGluR chain B
#
# Usage:
#   vmd -dispdev none -e loadall_ECDdist.tcl

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

    # Open output files for this replica
    set outfile1 [open LBDa_closure_dist_A155_E284_rep$k.txt w]  ;# LBD closure, chain A
    set outfile2 [open LBDb_closure_dist_A155_E284_rep$k.txt w]  ;# LBD closure, chain B
    set outfile3 [open LBD_229_dist_rep$k.txt w]                 ;# Inter-subunit dist at res 229
    set outfile4 [open CRD_553_dist_rep$k.txt w]                 ;# Inter-subunit dist at res 553
    set outfile5 [open LBDcom_dist_rep$k.txt w]                  ;# Lower lobe COM distance

    # CA atom selections for closure distance residues
    set res155a [atomselect top "segid PROA and resid 155 and name CA"]  ;# LBD lobe A, chain A
    set res284a [atomselect top "segid PROA and resid 284 and name CA"]  ;# LBD lobe B, chain A
    set res155b [atomselect top "segid PROD and resid 155 and name CA"]  ;# LBD lobe A, chain B
    set res284b [atomselect top "segid PROD and resid 284 and name CA"]  ;# LBD lobe B, chain B

    # CA atom selections for inter-subunit distance residues
    set res229a [atomselect top "segid PROA and resid 229 and name CA"]
    set res229b [atomselect top "segid PROD and resid 229 and name CA"]
    set res553a [atomselect top "segid PROA and resid 553 and name CA"]
    set res553b [atomselect top "segid PROD and resid 553 and name CA"]

    # Lower lobe of the LBD: includes residues 195-333 and 469-511 (CA only)
    # Used for inter-subunit COM distance calculation
    set LowLobeA [atomselect top \
        "segid PROA and resid 469 to 511 and name CA \
         or segid PROA and resid 195 to 333 and name CA"]
    set LowLobeB [atomselect top \
        "segid PROD and resid 469 to 511 and name CA \
         or segid PROD and resid 195 to 333 and name CA"]

    # --- Per-frame analysis ---
    for {set i 0} {$i < $nf} {incr i} {

        # Update all selections to the current frame
        $res155a frame $i
        $res155b frame $i
        $res284a frame $i
        $res284b frame $i
        $res229a frame $i
        $res229b frame $i
        $res553a frame $i
        $res553b frame $i
        $LowLobeA frame $i
        $LowLobeB frame $i

        # Get mass-weighted centers (single-atom selections return COM = atom position)
        set coord1  [measure center $res155a weight mass]
        set coord2  [measure center $res284a weight mass]
        set coord3  [measure center $res155b weight mass]
        set coord4  [measure center $res284b weight mass]
        set coord5  [measure center $res229a weight mass]
        set coord6  [measure center $res229b weight mass]
        set coord7  [measure center $res553a weight mass]
        set coord8  [measure center $res553b weight mass]
        set coord9  [measure center $LowLobeA weight mass]
        set coord10 [measure center $LowLobeB weight mass]

        # Compute distances using vecdist (returns scalar directly)
        set v1 [vecdist $coord2 $coord1]   ;# LBD closure, chain A: res284-res155
        set v2 [vecdist $coord4 $coord3]   ;# LBD closure, chain B: res284-res155
        set v3 [vecdist $coord6 $coord5]   ;# Inter-subunit distance at res 229
        set v4 [vecdist $coord8 $coord7]   ;# Inter-subunit distance at res 553
        set v5 [vecdist $coord10 $coord9]  ;# Inter-subunit lower lobe COM distance

        # Write distances to output files
        puts $outfile1 $v1
        puts $outfile2 $v2
        puts $outfile3 $v3
        puts $outfile4 $v4
        puts $outfile5 $v5
    }

    # Close output files for this replica
    close $outfile1
    close $outfile2
    close $outfile3
    close $outfile4
    close $outfile5

    # Clear trajectory frames before loading the next replica
    animate delete all
}

exit
