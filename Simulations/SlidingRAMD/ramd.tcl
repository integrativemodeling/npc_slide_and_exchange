#*****************************************************************************
#
# Random Acceleration Molecular Dynamics (RAMD)
#  - a 'tclForces' implementation for NAMDv2.5+
# 
# Author:  Harish Vashisth
# 
# (c) 2008 Board of Trustees of Drexel University
#
# Disclaimer:  This code for research purposes only.   
#
#*****************************************************************************
# modified by Jerry Karp, 10.2.14

# start
     
# make a list of all the ligand atoms to which force needs to be applied
# atoms indices are taken from PDB files

       set atoms {}
       for { set i 4047} { $i <= 4104  } { incr i } {
              lappend atoms $i
    }
       foreach atom $atoms {
        addatom $atom
    }

# make a group of ligand atoms on whose center of mass force is to be applied

   set groupid [addgroup $atoms]
   print "FSFG GROUP ID IS: $groupid "
  

# make a group of protein atoms
       set numatoms1 3994
       set atoms1 {}
       for { set i 1} { $i <= $numatoms1 } { incr i } {
              lappend atoms1 $i
       }
       set groupid1 [addgroup $atoms1]
       print "NTF2's GROUP ID IS: $groupid1 "


# some initializations   
# `forcefreq' and `rmin' are read from NAMD config. file and can be defined as below:
set forcefreq 20
set rmin      0.006

         set force     "0.0 0.0 0.0"

# a counter for force and time
         set forcecount $forcefreq
         set time      0

# old center of mass position
         set comold    "0.0 0.0 0.0"

         set totalwork 0
         set totaldist 0

# expulsion parameter which takes values either `0' or `1'
         set expulsion 0

# distance from initial position at which ligand can be considered exited from protein
# it is a system-specific parameter and should be different for problem at hand
         set R         65.0
         set dold      "0.0 0.0 0.0"
         set displace  0


#*******************************************************************************************   
# ACTUAL PROCEDURE WHICH IS CALLED AT EACH TIME STEP FOR FORCE APPLICATION FROM CONFIG. FILE
#*******************************************************************************************

proc calcforces { } {

      
# some global variables defintion 

         global atoms 
         global groupid force 
         global forcefreq
         global forcecount
         global comold time totalwork num rmin totaldist
         global expulsion  R
         global dold displace
         global atoms1 numatoms1 groupid1

# calculate and print the distance ligand must move in 'N' steps
# info: a check for correct parameter (`forcefreq' and `rmin') reading
        
         set vmin [expr "($rmin)/($forcefreq)"]
         if {$time == 0} {
            print " R-minimum (A):    $rmin"
            print " V-minimum (A/fs): $vmin"
	   }

# terminate NAMD in case of successful expulsion or if the ligand is stuck in protein
# conditions when expulsion is set to `1' are defined at the end and are system-specific

    if {$expulsion == 1} {
       set process [pid]
       print "Process to be killed is: $process"
       exec kill $process
      } 

# run further only if the ligand is not ejected

 if {$expulsion != 1} {

# info: print force vector details
        set f1 [lindex $force 0]
        set f2 [lindex $force 1]
        set f3 [lindex $force 2]
        set fmag [expr "sqrt((($f1)*($f1)+($f2)*($f2)+($f3)*($f3)))"]
        print "Force vector added at time $time is: $force"
        print "Force constant = $fmag"

# apply force on COM of the ligand
# force applied to the group is automatically distributed in proportion to all atoms by NAMD

        addforce $groupid $force
     
# get coordinates and masses of all the atoms in the groups and lists
# atom lists and groups are defined outside this `proc', at the top of this code

        loadcoords coords
        loadmasses masses  
        

# info: calculate COM position after each force application
    
        set comnew  "$coords($groupid)"
        print "Ligand-COM is at:        $time   $comnew"
        print "Center of mass new at $time = $comnew"
        print "Center of mass old at $time = $comold"

# info: keep track of the ligand's COM distance from COM of protein
# this is again system-specific
# users can define their own conditions to track the position of the ligand in their system
       
        set protein "$coords($groupid1)"
        print "Protein-COM is at:       $time   $protein"
        set d [vecsub $comnew $protein]
        set d1 [lindex $d 0]
        set d2 [lindex $d 1]
        set d3 [lindex $d 2]
        set r0 [expr "sqrt((($d1)*($d1)+($d2)*($d2)+($d3)*($d3)))"]
        print "Ligand is so far:        $time   $r0"



# compute the displacement of the ligand COM and not the distance
        set rem0 [expr "$time % $forcefreq"]
	if {$rem0 == 0 } {
           set dnew "$comnew"
         }
        if {$rem0 == 0 && $time != 0} {
        set disp [vecsub $dnew $dold]
        set e [lindex $disp 0]
        set f [lindex $disp 1]
        set g [lindex $disp 2]
        set displace [expr "sqrt((($e)*($e)+($f)*($f)+($g)*($g)))"]    
        print "Displacement in $time steps is =    $displace"  

        }

# calculate work done, distance moved after first force application
# work done is just for our information and not needed for expulsion trajectories
# this step is carried out only after time counter is `1'
    
   if {$time != 0} {
        set dist [vecsub $comnew $comold]
        set a [lindex $dist 0]
        set b [lindex $dist 1]
        set c [lindex $dist 2]
        set distance [expr "sqrt((($a)*($a)+($b)*($b)+($c)*($c)))"]    
        print "Distance moved =     $distance"
        set work [expr "(($f1)*($a)+($f2)*($b)+($f3)*($c))"]
        print "Time step and work done =       $time  $work"
        set totalwork [expr "$totalwork + $work"]
        print "Time step and ligand distance and accumlated work  = $time   $r0    $totalwork"
                 
      }
     
       print "TEST PRINT $forcecount = $forcefreq && $displace < $rmin"



# recalculate the random force direction at every `forcefreq' steps
# simultaneously check for minimum velocity threshold too
  
 if { $forcecount == $forcefreq && $displace < $rmin } {     

       print "NEW RANDOM FORCE CHOSEN in Kcal/(mol)(A) at $time"
# choice of `k' (kcal/mol.A) is system-specific
# different `k' values should be tried for successful expulsions
# we have used `k' = 20, 10, 15, and 5 (given in Table 1 as `f_0' values)
       set k   20.0
# random choice of Euler angles for x,y,z directions
       set pi  [expr "2.0*asin(1.0)"]
       set p   [expr "rand()"]
       set phi [expr "(($pi*$p)-($pi)/(2))"]
       set t   [expr "rand()"]
       set theta [expr "((2*$pi*$t)-($pi))"]
       set rx [expr "cos($theta)*cos($phi)"]
       set ry [expr "sin($theta)*cos($phi)"]
       set rz [expr "sin($phi)"]
       set r "$rx $ry $rz"
       set f [vecscale [expr "$k"] $r ]

# compute the force vector

       set force [vecscale [expr 1.0] $f]

# increment force counter
       set forcecount 0
       
    }


# print in case the ligand is already in solvent
# we define it as a distance from initial positions and set counter to `1'
# it is system-specific and other conditions can be used as well
        
       if {$r0 >= $R && $expulsion == 0 } {
          print "EJECTED AT TIME $time   $r0"
          set expulsion "1"
	 }  


# a distance threshold is defined if the ligand is assumed stuck in protein
# this is also system-specific and not necessarily has to be as given here
        
      if {$r0 <= 11.5} {
          print "Moving to center at $time"
          set expulsion "1"
	 } 


# some counter increments

        set comold "$comnew"
        if {$rem0 == 0 } {
           set dold "$dnew"
           set forcecount 0
         }
     
        incr time
     	incr forcecount
        return

   }

  }


#*****************
# END
#*****************




