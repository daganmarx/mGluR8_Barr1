# loadall_xyz.tcl
#
# Description:
#   Batch analysis script for VMD. Iterates over a series of MD trajectory
#   replicas (DCD files) and computes two angular metrics per frame:
#
#   1. XY-plane angle (adj_ang): The angle between the arrestin inter-lobe
#      vector (N-lobe COM -> C-lobe COM, projected onto XY plane) and the
#      principal axis of the mGluR chain A (also projected onto XY plane).
#      Reported as an adjusted angle (90 - raw angle), written to:
#          Barr1TMD_newxy_angle_rep<k>.dat
#
#   2. Z-tilt angle (ang_deg_z): The angle between the arrestin inter-lobe
#      vector and the Z-axis (i.e., tilt out of the membrane plane).
#      Adjusted so that 0 degrees = membrane-parallel orientation.
#      Written to:
#          Barr1chainA_z_angle_rep<k>.dat
#
# Requirements:
#   - VMD with the Orient and tempoUserVMD packages installed
#   - A PSF file: now_0000.psf
#   - DCD trajectory files named: trajectory_now_<k>.dcd (k = 20 to 71)
#
# Segment IDs used:
#   PROF : Arrestin-1 (Barr1) — residues 1-175 = N-lobe, 176-360 = C-lobe
#   PROA : mGluR chain A — residues 580-900 used for principal axis calculation
#
# Usage:
#   Run from within VMD:
#       vmd -dispdev none -e loadall_xyz.tcl
#   or source interactively:
#       source loadall_xyz.tcl

# Load the structural topology (PSF) file
mol load psf now_0000.psf

# Load required VMD packages
package require Orient        ;# Provides principal axis calculation utilities
package require tempoUserVMD  ;# Provides dopbc for PBC wrapping/correction
namespace import Orient::orient

# --- Main loop over trajectory replicas ---
for { set k 20 } { $k <= 71 } { incr k 1 } {

    # Load trajectory k, applying periodic boundary correction.
    # Frames are wrapped by fragment, aligned to segid PROA (mGluR chain A).
    dopbc -file trajectory_now_$k.dcd -frames 0:1:1000 \
          -ref "segid PROA" -wrapby "fragment"

    # Open output files for this replica
    set outfile1 [open Barr1TMD_newxy_angle_rep$k.dat w]  ;# XY-plane angle
    set outfile2 [open Barr1chainA_z_angle_rep$k.dat w]   ;# Z-tilt angle

    set nf [molinfo top get numframes]  ;# Total number of loaded frames

    # Define atom selections
    # Arrestin N-lobe: residues 1-175 of segid PROF
    set Nlobe  [atomselect top "segid PROF and resid 1 to 175"]
    # Arrestin C-lobe: residues 176-360 of segid PROF
    set Clobe  [atomselect top "segid PROF and resid 176 to 360"]
    # mGluR chain A: residues 580-900 of segid PROA
    set chainA [atomselect top "segid PROA and resid 580 to 900"]

    # --- Per-frame analysis ---
    for { set i 1 } { $i <= $nf } { incr i 1 } {

        # Update each selection to the current frame
        $Nlobe  frame $i
        $Clobe  frame $i
        $chainA frame $i

        # Compute centers of mass for arrestin lobes and mGluR chain A
        set Nlobe_com  [measure center $Nlobe  weight mass]
        set Clobe_com  [measure center $Clobe  weight mass]
        set chainA_com [measure center $chainA weight mass]

        # Arrestin inter-lobe vector: points from C-lobe COM to N-lobe COM
        set barr1 [vecsub $Nlobe_com $Clobe_com]

        # mGluR principal axis 2 (index 1): long axis of the receptor in XY plane
        set mGluR [lindex [draw principalaxes $chainA] 1]

        # Project both vectors onto the XY plane (set Z component to 0)
        set barr1_xy [lreplace $barr1 2 2]
        set mGluR_xy [lreplace $mGluR 2 2]

        # Reference Z-axis vector for tilt calculation
        set zvec "0 0 1"

        # Vector magnitudes for angle normalization
        set length2 [veclength $mGluR_xy]   ;# |mGluR XY projection|
        set length3 [veclength $barr1_xy]   ;# |arrestin XY projection|

        # --- Z-tilt angle: angle between arrestin vector and Z-axis ---
        set dot_prod_z [vecdot $barr1 $zvec]
        set cos_z      [expr "$dot_prod_z / $length3"]
        set ang_rad_z  [tcl::mathfunc::acos "$cos_z"]
        # Convert to degrees; subtract from 180 so that 0 = parallel to membrane
        set ang_deg_z  [expr "180 - 57.2957795 * $ang_rad_z"]

        # --- XY-plane angle: angle between arrestin and mGluR principal axis ---
        set dot_prod_xy [vecdot $barr1_xy $mGluR_xy]
        set cos_xy      [expr "$dot_prod_xy / $length3 / $length2"]
        set ang_rad_xy  [tcl::mathfunc::acos "$cos_xy"]
        # Convert to degrees
        set ang_deg_xy  [expr "57.2957795 * $ang_rad_xy"]
        # Adjust so that 0 = aligned with receptor axis
        set adj_ang     [expr "90 - $ang_deg_xy"]

        # Print to terminal and write to output files
        puts $adj_ang
        puts $ang_deg_z
        puts $outfile1 $adj_ang
        puts $outfile2 $ang_deg_z
    }

    # Close output files for this replica
    close $outfile1
    close $outfile2

    # Clear all trajectory frames before loading the next replica
    animate delete all
}

exit
