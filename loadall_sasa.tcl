# loadall_sasa.tcl
#
# Description:
#   Batch analysis script for VMD. Iterates over a specified trajectory replica
#   and computes the solvent-accessible surface area (SASA) buried by contact
#   between the arrestin-1 N-lobe loops and the mGluR transmembrane domain
#   (chain A).
#
#   Buried SASA is computed as:
#       buried_sa = SASA(loops alone) - SASA(loops in context of chainA)
#
#   A larger buried SASA indicates greater contact/occlusion of the loops
#   by the receptor TMD surface.
#
#   The probe radius used is 1.4 Å (standard water probe).
#   Frames are sampled every 2 steps (incr i 2) to reduce compute time.
#
#   Loops analyzed:
#     - Arrestin N-lobe loops: residues 238-248 and 305-317 of segid PROF
#
#   Contact partner:
#     - mGluR chain A TMD + arrestin (segid PROA res 580-850 + segid PROF)
#       (arrestin is included in the SASA surface probe to account for
#        self-occlusion)
#
#   Output (per replica):
#       Nloops_chainA_rep<k>.txt — buried SASA (Å²), one value per sampled frame
#
# Requirements:
#   - VMD with Orient and tempoUserVMD packages installed
#   - PSF file: now_0000.psf
#   - DCD trajectory files: trajectory_now_<k>.dcd
#
# Segment IDs used:
#   PROF : Arrestin-1 (Barr1)
#   PROA : mGluR chain A — TMD residues 580-850
#
# Note:
#   - The script is currently set to run only replica 39 (k = 39 to 39).
#     Change the loop bounds to analyze other replicas.
#   - Many additional SASA metrics (Nlobe, membrane, ICL regions, etc.) are
#     implemented but commented out. Uncomment to enable.
#   - SASA is computed every other frame (step 2) for performance; adjust
#     `incr i 2` to `incr i 1` for full-resolution analysis.
#
# Usage:
#   vmd -dispdev none -e loadall_sasa.tcl

# Load the structural topology (PSF) file
mol load psf now_0000.psf

# Load required VMD packages
package require Orient        ;# Principal axis utilities
package require tempoUserVMD  ;# Provides dopbc for PBC correction
namespace import Orient::orient

# --- Main loop over trajectory replicas (currently: replica 39 only) ---
for { set k 39 } { $k <= 39 } { incr k 1 } {

    # Load trajectory with PBC correction, centering on lipid bilayer
    dopbc -file trajectory_now_$k.dcd -frames 0:1:1000 \
          -ref "lipid" -wrapby "fragment"

    # --- Atom selections ---
    # Arrestin N-lobe loops: two loop regions used for SASA computation
    set Nloops [atomselect top \
        "(segid PROF and resid 238 to 248) or (segid PROF and resid 305 to 317)"]

    # Contact surface: mGluR chain A TMD combined with the full arrestin
    # (used as the reference surface for buried SASA calculation)
    set chainA [atomselect top "(segid PROA and resid 580 to 850) or (segid PROF)"]

    set nf [molinfo top get numframes]  ;# Total number of loaded frames

    # Full arrestin protein selection (used as the probe surface for isolated SASA)
    set prot [atomselect top "protein and segid PROF"]

    # Open output file for buried SASA
    set outfile3 [open Nloops_chainA_rep$k.txt w]

    # --- Per-frame SASA analysis (every other frame for performance) ---
    for {set i 0} {$i < $nf} {incr i 2} {
        puts $i  ;# Print frame number to terminal for progress tracking

        # Update selections to current frame
        $Nloops frame $i
        $chainA frame $i
        $prot   frame $i

        # SASA of the loops in isolation (probe = full arrestin protein surface)
        set Nloops_sasa [measure sasa 1.4 $prot -restrict $Nloops]

        # SASA of the loops in contact with the receptor TMD
        # (probe = chainA surface; buried area = difference)
        set Nloops_chainA_sasa [measure sasa 1.4 $chainA -restrict $Nloops]

        # Buried SASA = area hidden by receptor contact
        set buried_sa3 [expr $Nloops_sasa - $Nloops_chainA_sasa]

        # Write buried SASA to output file
        puts $outfile3 $buried_sa3
    }

    # Close output file for this replica
    close $outfile3

    # Clear trajectory frames before loading the next replica
    animate delete all
}

exit
