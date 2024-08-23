# This script identifies surface residues, specifically methionines on an Fc fragment of an antibody. These methionines are readily oxidized and the oxidized form is less stable. This script can also be modified to show the aromatic residues on the surface.
#
# Liu et al. (2008) Structure and stability changes of human IgG1 Fc as a consequence of methionine oxidation. Biochemistry, 47, 5088. https://doi.org/10.1021/bi702238b

# Instructions:

# open PyMOL
# drag this file onto the PyMOL viewing window, which will run all the commands below.

# To view the surface aromatic residues, uncomment (remove the starting '## ') the two lines containing 'Surf_Aromatics'.

# clear any PyMOL settings and set working directory
reinitialize
cd ~/Desktop

# the python definition below defines the command to find surface residues
python
'''
http://pymolwiki.org/index.php/FindSurfaceResidues
'''

from __future__ import print_function
from pymol import cmd


def findSurfaceAtoms(selection="all", cutoff=2.5, quiet=1):
    """
DESCRIPTION
    Finds those atoms on the surface of a protein
    that have at least 'cutoff' exposed A**2 surface area.
USAGE
    findSurfaceAtoms [ selection, [ cutoff ]]
SEE ALSO
    findSurfaceResidues
    """
    cutoff, quiet = float(cutoff), int(quiet)

    tmpObj = cmd.get_unused_name("_tmp")
    cmd.create(tmpObj, "(" + selection + ") and polymer", zoom=0)

    cmd.set("dot_solvent", 1, tmpObj)
    cmd.get_area(selection=tmpObj, load_b=1)

    # threshold on what one considers an "exposed" atom (in A**2):
    cmd.remove(tmpObj + " and b < " + str(cutoff))

    selName = cmd.get_unused_name("exposed_atm_")
    cmd.select(selName, "(" + selection + ") in " + tmpObj)

    cmd.delete(tmpObj)

    if not quiet:
        print("Exposed atoms are selected in: " + selName)

    return selName


def findSurfaceResidues(selection="all", cutoff=2.5, doShow=0, quiet=1):
    """
DESCRIPTION

    Finds those residues on the surface of a protein
    that have at least 'cutoff' exposed A**2 surface area.

USAGE

    findSurfaceResidues [ selection, [ cutoff, [ doShow ]]]

ARGUMENTS

    selection = string: object or selection in which to find exposed
    residues {default: all}

    cutoff = float: cutoff of what is exposed or not {default: 2.5 Ang**2}

RETURNS

    (list: (chain, resv ) )
        A Python list of residue numbers corresponding
        to those residues w/more exposure than the cutoff.

    """
    cutoff, doShow, quiet = float(cutoff), int(doShow), int(quiet)

    selName = findSurfaceAtoms(selection, cutoff, quiet)

    exposed = set()
    cmd.iterate(selName, "exposed.add((chain,resv))", space=locals())

    selNameRes = cmd.get_unused_name("exposed_res_")
    cmd.select(selNameRes, "byres " + selName)

    if not quiet:
        print("Exposed residues are selected in: " + selNameRes)

    if doShow:
        cmd.show_as("spheres", "(" + selection + ") and polymer")
        cmd.color("white", selection)
        cmd.color("yellow", selNameRes)
        cmd.color("red", selName)

    return sorted(exposed)

cmd.extend("findSurfaceAtoms", findSurfaceAtoms)
cmd.extend("findSurfaceResidues", findSurfaceResidues)
python end

# download the wt structure
fetch 5JII, type = pdb, name = 'WT', async = 0

# remove all except the protein in chain A
# the full Fc protein is a dimer (chain A and chain B), but chain B is removed for clarity
remove not chain A
remove not polymer.protein

# run the new command defined by the python script above
# the default is not to show the residues yet
# all the surface residues are now in an object named 'exposed_res_01'
findSurfaceResidues

# show only the methionines from all the surface residues
select Surf_Met, exposed_res_01 and resn Met
show sticks, Surf_Met

# use the commands below to show surface aromatic residues
## select Surf_Aromatics, exposed_res_01 and (resn Phe or resn Tyr or resn Trp)
## show sticks, Surf_Aromatics

# hide hydrogens and orient the view
hide (hydro)

set_view (\
     0.731381595,   -0.593263984,    0.336333036,\
    -0.554532349,   -0.230288655,    0.799664259,\
    -0.396957517,   -0.771367610,   -0.497412682,\
     0.000000000,   -0.000000000, -142.723129272,\
    40.463905334,   48.197937012,   39.846763611,\
   112.524002075,  172.922256470,  -20.000000000 )

deselect
