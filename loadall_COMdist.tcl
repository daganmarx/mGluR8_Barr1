# loadall_COMdist.tcl
#
# Description:
#   Batch analysis script for VMD. Iterates over all trajectory replicas and
#   computes three center-of-mass (COM) distance metrics per frame:
#
#   1. CRD inter-subunit distance: distance between the COMs of the cysteine-
#      rich domain (CRD) of mGluR chain A (PROA) and chain B (PROD).
#      Written to: CRD_com_rep<k>.txt
#
#   2. TMD inter-subunit distance: distance between the COMs of the
#      transmembrane domain (TMD) of chain A and chain B.
#      Written to: TMD_com_rep<k>.txt
#
#   3. Finger-loop to TMD distance: distance between the COM of the arrestin
#      finger loop (FL, residues 62-74 of PROF) and the TMD COM of chain A.
#      Written to: FL_com_rep<k>.txt
#
# Requirements:
#   - VMD with Orient and tempoUserVMD packages installed
#   - PSF file: now_0000.psf
#   - DCD trajectory files: trajectory_now_<k>.dcd (k = 0 to 71)
#
# Segment IDs used:
#   PROA : mGluR chain A
#       - residues 500-579: CRD
#       - residues 580-850: TMD
#   PROD : mGluR chain B
#       - residues 500-579: CRD
#       - residues 580-850: TMD
#   PROF : Arrestin-1 (Barr1)
#       - residues 62-74: finger loop
#
# Note:
#   PBC wrapping uses lipid as the reference (rather than PROA), which centers
#   the system on the membrane rather than the receptor.
#
# Usage:
#   vmd -dispdev none -e loadall_COMdist.tcl

# Load the structural topology (PSF) file
mol load psf now_0000.psf

# Load required VMD packages
package require Orient        ;# Principal axis utilities
package require tempoUserVMD  ;# Provides dopbc for PBC correction
namespace import Orient::orient

# --- Main loop over trajectory replicas ---
for { set k 0 } { $k <= 71 } { incr k 1 } {

    # Load trajectory with PBC correction, centering on lipid bilayer
    dopbc -file trajectory_now_$k.dcd -frames 0:1:1000 \
          -ref "lipid" -wrapby "fragment"

    # Open output files for this replica
    set outfile1 [open CRD_com_rep$k.txt w]   ;# CRD inter-subunit COM distance
    set outfile2 [open TMD_com_rep$k.txt w]   ;# TMD inter-subunit COM distance
    set outfile3 [open FL_com_rep$k.txt w]    ;# Finger loop to TMD-A COM distance

    set nf [molinfo top get numframes]  ;# Total number of loaded frames

    # Define atom selections for mGluR domains and arrestin finger loop
    set CRDa [atomselect top "protein and segid PROA and resid 500 to 579"]  ;# CRD, chain A
    set CRDb [atomselect top "protein and segid PROD and resid 500 to 579"]  ;# CRD, chain B
    set TMDa [atomselect top "protein and segid PROA and resid 580 to 850"]  ;# TMD, chain A
    set TMDb [atomselect top "protein and segid PROD and resid 580 to 850"]  ;# TMD, chain B
    set FL   [atomselect top "protein and segid PROF and resid 62 to 74"]    ;# Arrestin finger loop

    # --- Per-frame analysis ---
    for { set i 1 } { $i <= $nf } { incr i 1 } {

        # Update all selections to the current frame
        $CRDa frame $i
        $CRDb frame $i
        $TMDa frame $i
        $TMDb frame $i
        $FL   frame $i

        # Compute mass-weighted centers of mass for each domain
        set CRDa_com [measure center $CRDa weight mass]
        set CRDb_com [measure center $CRDb weight mass]
        set TMDa_com [measure center $TMDa weight mass]
        set TMDb_com [measure center $TMDb weight mass]
        set FL_com   [measure center $FL   weight mass]

        # Compute inter-domain difference vectors
        set CRD_com_dist     [vecsub $CRDa_com $CRDb_com]   ;# CRD A-B separation
        set TMD_com_dist     [vecsub $TMDa_com $TMDb_com]   ;# TMD A-B separation
        set FL_TMDa_com_dist [vecsub $FL_com   $TMDa_com]   ;# FL to TMD-A distance

        # Compute scalar distances (vector magnitudes)
        set lenCRDcom    [veclength $CRD_com_dist]
        set lenTMDcom    [veclength $TMD_com_dist]
        set lenFLTMDcom  [veclength $FL_TMDa_com_dist]

        # Write distances to output files
        puts $outfile1 $lenCRDcom
        puts $outfile2 $lenTMDcom
        puts $outfile3 $lenFLTMDcom
    }

    # Close output files for this replica
    close $outfile1
    close $outfile2
    close $outfile3

    # Clear trajectory frames before loading the next replica
    animate delete all
}

exit
