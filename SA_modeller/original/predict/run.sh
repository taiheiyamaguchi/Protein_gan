#!/bin/tcsh
set nsteps=10000
set log_frequency=100
set trajectory_frequency=100
set restraint_in="rest.dat"
set log_out=`printf "anneal.log"`
set trajectory_out=`printf "tra.pdb"`
set final_out=`printf "final.pdb"`
../anneal/anneal \
	$nsteps \
	$log_out \
	$log_frequency \
	$restraint_in \
	$trajectory_out \
	$trajectory_frequency \
	$final_out
end
