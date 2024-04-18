set var [lindex $argv 0]
set ligname [lindex $argv 1]

mol new ../$var-solvate.prmtop 
mol addfile ../$var-solvate_output.xtc type xtc first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all

package require pbctools
pbc wrap -all -compound res -center bb -centersel "protein or nucleic"

set frame0 [atomselect top "noh and (protein or resname $ligname)" frame 0]
set lig [atomselect top "noh and resname $ligname"]
set prot [atomselect top "noh and protein"]

set nf [molinfo top get numframes]

for {set i 1 } {$i < $nf } { incr i } {
    set sel [atomselect top "noh and (protein or resname $ligname)" frame $i]
    set all [atomselect top all frame $i]
    $all move [measure fit $sel $frame0]
    puts [open rmsd_protein.dat w] "[measure rmsd $prot $frame0]"
    puts [open rmsd_ligand.dat w] "[measure rmsd $lig $frame0]"

}

quit
