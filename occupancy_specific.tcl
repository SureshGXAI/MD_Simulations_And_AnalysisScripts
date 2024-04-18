set DIR [lindex $argv 0]
set trj $DIR-solvate_output.xtc
set prmtop $DIR-solvate.prmtop
mol new $trj type xtc first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all
mol addfile $prmtop type parm7 first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all

set top [atomselect top "all"]
set prot [atomselect top "protein and resid 50"]
set lig [atomselect top "protein and resid 67"]

source /usr/local/lib/vmd/plugins/noarch/tcl/hbonds1.2/hbonds.tcl 
hbonds -sel1 $prot -sel2 $lig -type unique -writefile yes -dist 3.8 -ang 120 -polar yes -detailout $DIR-details_occupancy_specific.dat

