
 proc _vcos { v1 v2} {
 #
 # returns cos theta between v1 and v2
 #
     set v1n [vecnorm $v1]
     set v2n [vecnorm $v2]
     return [vecdot $v1n $v2n]
 }

proc _cosNH {seltext t tt {molid top}} {
  set Nt [[atomselect $molid "$seltext and name N" frame $t] get "x y z"]   
  set Ht [[atomselect $molid "$seltext and name H" frame $t] get "x y z"] 
  set Ntt [[atomselect $molid "$seltext and name N" frame $tt] get "x y z"]   
  set Htt [[atomselect $molid "$seltext and name H" frame $tt] get "x y z"] 
  set dt [vecsub [lindex $Nt 0] [lindex $Ht 0]]
  set dtt [vecsub [lindex $Ntt 0] [lindex $Htt 0]]
  return [_vcos $dt $dtt]
}

proc _P2 {i dt {step 10} {i0 0} {i1 -1} {molid top} } {
    # Returns the mean second-rank Legendre polynomial of the
    # correlation between the N-H vectors (time correlation
    # function)
    #
    # Params:
    # i - resid number (not residue!)
    # dT - correlation over dT frames
    # step - skip rate when computing correlation average (in frames)
    # i0 - first from-frame from which to calculate correlation
    # i1 - last from-frame for which to calculate correlation, or -1 for end
    set seltext "resid $i"
    if { [[atomselect $molid "$seltext and name H"] num] == 0 } {
	return 0.0
    }
    set n [expr [molinfo top get numframes] - $dt]
    if {$i0 >= $n } {
	error "_P2: from frame $i0 plus delta-T $dt is beyond last frame [expr $n-1]"
    }
    if { ($i1 >= 0) && ($i1 < $n) } {
	set n $i1
    }
    set P2sum 0.0
    set nsteps 0
    for { set i_from $i0 } { $i_from < $n } {set i_from [expr $i_from+$step]} {
	set i_to [expr $i_from+$dt];
  	set c [_cosNH $seltext $i_from $i_to $molid]
	set P2c [expr 0.5 * (3 * $c * $c - 1.0)]
	## puts [format "%6.0f %6.0f %6.0f %8.2f" $dt $i_from $i_to $P2c]
	set P2sum [expr $P2sum + $P2c]
	incr nsteps
    }
    set P2 [expr $P2sum / $nsteps]
    return $P2
}

proc _P2_multiple {i dt {step 10} {f0 0} {f1 -1} {molids "top"}} {
    # Returns the mean second-rank Legendre polynomial of the
    # correlation between the N-H vectors (time correlation
    # function) averaged over a list of mol_ids
    #
    # Params:
    # i - resid number (not residue!)
    # dT - correlation over dT frames
    # stepsize - skip rate when computing correlation average
    # mol_ids list of molecular ids
    # f0 - first frame from which to calculate correlation
    # f1 - last frame for which to calculate correlation, or -1 for end
    puts "residue (not resid) $i dt $dt step $step f0 $f0 f1 $f1 molids $molids"
    puts [[atomselect [lindex $molids 0] "residue $i and name CA"] get resname]
    set W 0.0 
    set sumP2 0.0
    for {set ii 0} {$ii < [llength $molids]} {incr ii} {
       set molid [lindex $molids $ii]
       set w [molinfo $molid get numframes]
       set curP2 [_P2 $i $dt $step $f0 $f1 $molid]
       set sumP2 [expr $sumP2 + $w * $curP2]
       set W [expr $W + $w]
    }
    return [expr $sumP2 / $W]
}

proc plot_P2 { dT {stepsize 1} {f0 0} {f1 -1} {from_res 1} {to_res -1} {molids "top"} {p -1} } {
    # dT - correlation over dT frames
    # stepsize - skip rate when computing correlation average
    # f0 - first frame from which to calculate correlation
    # f1 - last frame for which to calculate correlation, or -1 for end
    # from_resid - first resid (not residue!) in range
    # to_resid - last resid (not residue!) in range (or -1 for last residue)
    # p - plot object or -1 for new plot
     
    set x [list]
    set y [list]
    set CAs [atomselect top "protein and name CA"]
    set ress [$CAs get resid]
    for {set i 0} { $i < [$CAs num] } {incr i} {
	set res [lindex $ress $i]
	if { $res < $from_res || 
	     $res > $to_res && 
	     $to_res != -1 } {
	    continue 
	}
	set P2 [_P2_multiple $res $dT $stepsize $f0 $f1 $molids]
	if { $P2 == 0.0 } {
		continue
	}
	set resid [expr $res+1]
	lappend x $resid
	lappend y $P2
    }
    puts $x
    set legend "$f0 to $f1"
    puts "legend=$legend"
    if {$p == -1} {
	set p [multiplot -x $x -y $y —legend $legend -xlabel "resid" -marker circle -ymin 0.0 -ymax 1.0 -plot]
    } else {
	$p add $x $y —legend $legend -plot
    }
    return $p
}


# calculate P2 for residue (not resid) k for time interval $i*$ns_per_frame with 
# sampling rate $stepf frames
proc _do_k_i { fp k i {ns_per_frame 0.06} {stepf 100} } {
	set t [expr $i*$ns_per_frame]
	set p2 [_P2 $k $i $stepf 0 -1];
	lappend x $t; lappend y $p2;  
	puts $fp [format "%.3f %.4f" $t $p2];
}

# stepf - skip rate in frames
proc get_resid_P2 { resid  {ns_per_frame 0.06} {stepf 100}} {
	set x [list]; set y [list]; 
	set fname [format "%d.txt" $resid]
	set fp [open $fname w]
	for {set i 0} {$i<100} {incr i} {
		_do_k_i $fp $resid $i $ns_per_frame $stepf;
	}; 
	flush $fp;
	for {set i 100} {$i<1000} {set i [expr $i+10]} {
		_do_k_i $fp $resid $i $ns_per_frame $stepf;
	}
	flush $fp;	
	for {set i 1000} {$i<4000} {set i [expr $i+100]} {
		_do_k_i $fp $resid $i $ns_per_frame $stepf;
	}
	close $fp;	
}




proc get_resid_P2_intervals { resid {ns_per_frame 0.06} {stepf 25} } {
    # Computes P2(dT) time autocorrelation function for different
    # values of dT, for each of nine partially overlapping time
    # segments (to later estimate std-errors).
    #
    # Params:
    #   resid        - the residue id
    #   ns_per_frame - the nanoseconds per frame rate of the trajectory
    #   stepf - the frame skip rate when computing the P@ value for
    #   each value of dT. A high <stepf> value can speed the run
    #   considerably, but should be kept low enough so as not to
    #   sacrifice accuracy (this can be tested by robustness of
    #   results)
    #
    # Output is writting to file "<resid>.txt", first column is dt,
    # columns 2-10 are P2s for each time segment
    set fname [format "%d.txt" $resid]
    set fp [open $fname w]
    set n [	molinfo top get numframes];
    set k 5; # number of intervals is k*2-1, with segments of size n/k, in floor(n/k/2.0) increments
    set interval [expr round($n/$k)];
    set max_dF [expr $interval-10*$stepf]; 
    ## set max_dF [expr round(250.0/$ns_per_frame)]
    # df = delta-frames
    for {set dF 0} {$dF<$max_dF} {set dF [expr ceil(($dF+0.1)*1.03)]} { 
	set dNS [expr $dF*$ns_per_frame];
	puts -nonewline $fp [format {%9.2f } $dNS];
	puts -nonewline    "[format {%9.2f } $dNS] ";
	for {set from 0} {$from<=[expr $n-$interval]} {set from [expr $from+floor($interval/2.0)]} {
	    set to [expr $from+$interval];
	    if { [expr $from+$dF]>=$n } {
		set from_final [expr $n-$dF-1];
	    } else {
		set from_final $from
	    }
	    set p2 [_P2 $resid $dF $stepf $from_final $to];
	    puts -nonewline $fp [format {%5.3f } $p2];
	}
	puts $fp ""; 
	puts "";
    }
    close $fp
}

