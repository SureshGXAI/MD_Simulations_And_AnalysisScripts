set DIR [lindex $argv 0]
set ligname [lindex $argv 1]
set trj $DIR-solvate_output.xtc
set prmtop $DIR-solvate.prmtop
mol new $trj type xtc first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all
mol addfile $prmtop type parm7 first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all

set top [atomselect top "all"]
set prot [atomselect top "protein"]
set lig [atomselect top "resname $ligname"]

source /usr/local/lib/vmd/plugins/noarch/tcl/hbonds1.2/hbonds.tcl 
hbonds -sel1 $prot -sel2 $lig -type unique -writefile yes -dist 3.5 -ang 120 -polar yes -detailout $DIR-details_occupancy.dat

