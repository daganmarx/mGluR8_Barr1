# loadall_Barractangle.tcl
#
# Description:
#   Batch analysis script for VMD. Iterates over a series of MD trajectory
#   replicas (DCD files) and computes the inter-lobe bending angle of
#   arrestin-1 (Barr1) for each frame.
#
#   The bending angle is defined as the angle between the first principal axes
#   of the arrestin N-lobe and C-lobe. This captures the relative orientation
#   (opening/closing or twisting) of the two lobes throughout the simulation.
#
#   An angle > 100 degrees is reflected to its supplement (180 - angle) to
#   enforce a consistent 0-90 degree convention and avoid discontinuities
#   from principal axis sign ambiguity.
#
#   Output (per replica):
#       Barr1_prinax1_angle_rep<k>.dat — inter-lobe angle (degrees), one value per frame
#
# Requirements:
#   - VMD with the Orient and tempoUserVMD packages installed
#   - A PSF file: now_0000.psf
#   - DCD trajectory files named: trajectory_now_<k>.dcd (k = 0 to 71)
#
# Segment IDs used:
#   PROF : Arrestin-1 (Barr1)
#       - residues 1-175:   N-lobe
#       - residues 176-360: C-lobe
#
# Note:
#   Additional output files for principal axes 2 and 3 are implemented but
#   commented out. Uncomment the relevant blocks to enable them.
#
# Usage:
#   Run headlessly from the command line:
#       vmd -dispdev none -e loadall_Barractangle.tcl
#   or source interactively within VMD:
#       source loadall_Barractangle.tcl

# Load the structural topology (PSF) file
mol load psf now_0000.psf

# Load required VMD packages
package require Orient        ;# Provides principal axis calculation utilities
package require tempoUserVMD  ;# Provides dopbc for PBC wrapping/correction
namespace import Orient::orient

# --- Main loop over trajectory replicas ---
for { set k 0 } { $k <= 71 } { incr k 1 } {

    # Load trajectory k with PBC correction.
    # Frames are wrapped by fragment, aligned to segid PROA (mGluR chain A).
    dopbc -file trajectory_now_$k.dcd -frames 0:1:1000 \
          -ref "segid PROA" -wrapby "fragment"

    # Open output file for principal axis 1 inter-lobe angle (this replica)
    set outfile1 [open Barr1_prinax1_angle_rep$k.dat w]
    # Uncomment below to also output axes 2 and 3:
    #set outfile2 [open Barr1_prinax2_angle_rep$k.dat w]
    #set outfile3 [open Barr1_prinax3_angle_rep$k.dat w]

    set nf [molinfo top get numframes]  ;# Total number of loaded frames

    # Define atom selections for arrestin lobes
    set Nlobe [atomselect top "segid PROF and resid 1 to 175"]    ;# N-lobe
    set Clobe [atomselect top "segid PROF and resid 176 to 360"]  ;# C-lobe

    # --- Per-frame analysis ---
    for { set i 1 } { $i <= $nf } { incr i 1 } {

        # Update lobe selections to the current frame
        $Nlobe frame $i
        $Clobe frame $i

        # Compute the first principal axis (longest inertia axis) for each lobe.
        # Index 0 = first principal axis (largest eigenvalue).
        set prinax1N [lindex [draw principalaxes $Nlobe] 0]
        set prinax1C [lindex [draw principalaxes $Clobe] 0]
        # Uncomment to compute additional principal axes:
        #set prinax2N [lindex [draw principalaxes $Nlobe] 1]
        #set prinax2C [lindex [draw principalaxes $Clobe] 1]
        #set prinax3N [lindex [draw principalaxes $Nlobe] 2]
        #set prinax3C [lindex [draw principalaxes $Clobe] 2]

        # Magnitudes of principal axis vectors (should be ~1 for unit vectors)
        set length2  [veclength $prinax1N]
        set length2b [veclength $prinax1C]

        # Dot product between N-lobe and C-lobe first principal axes
        set dot_prod_1 [vecdot $prinax1N $prinax1C]
        #set dot_prod_2 [vecdot $prinax2N $prinax2C]
        #set dot_prod_3 [vecdot $prinax3N $prinax3C]

        # Cosine of the inter-lobe angle (normalized dot product)
        set cos_1 [expr "$dot_prod_1 / $length2 / $length2b"]
        #set cos_2 [expr "$dot_prod_2 / $length3 / $length3b"]
        #set cos_3 [expr "$dot_prod_3 / $length4 / $length4b"]

        # Convert to angle in radians, then degrees
        set ang_rad_1  [tcl::mathfunc::acos "$cos_1"]
        set ang_deg_1  [expr "57.2957795 * $ang_rad_1"]
        #set ang_rad_2 [tcl::mathfunc::acos "$cos_2"]
        #set ang_deg_2 [expr "57.2957795 * $ang_rad_2"]
        #set ang_rad_3 [tcl::mathfunc::acos "$cos_3"]
        #set ang_deg_3 [expr "57.2957795 * $ang_rad_3"]

        # Reflect angles > 100 degrees to their supplement.
        # Principal axes have sign ambiguity (axis and its negative are equivalent),
        # so angles > 90 degrees are folded back into the 0-90 range.
        if {$ang_deg_1 > 100} {
            set ang_deg_1 [expr 180 - $ang_deg_1]
        }
        #if {$ang_deg_2 > 100} { set ang_deg_2 [expr 180 - $ang_deg_2] }
        #if {$ang_deg_3 > 100} { set ang_deg_3 [expr 180 - $ang_deg_3] }

        # Write result to output file
        puts $outfile1 $ang_deg_1
        #puts $outfile2 $ang_deg_2
        #puts $outfile3 $ang_deg_3
    }

    # Close output file for this replica
    close $outfile1
    #close $outfile2
    #close $outfile3

    # Clear all trajectory frames before loading the next replica
    animate delete all
}

exit
