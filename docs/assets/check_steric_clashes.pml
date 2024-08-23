## Check a protein for steric clashes ##

# This script identifies a small steric clash created by the T4 lysozyme A129F substitution. First, it shows the wt protein (pdb code: 1l63) with no steric clash near Ala129. Next, modify the script to show A129F (pdb code: 1QTC), which shows a little steric clash near Phe129, but less than might be expected because the surrounding region adjusts for the addition of the large phenyl group.

# T4 lysozyme A129F is nevertheless destabilized relative to the wild type. The melting temperature of the A129F variant (94 °C) is 19 °C lower than that for the wild type (113 °C). The measured free energy of unfolding is 1.2 kcal/mol lower at 40 °C for the A129F variant indicating that it is less stable.

# Liu et al. (2000) The introduction of strain and its effects on the structure and stability of T4 lysozyme. J. Mol. Biol., 295, 127. https://doi.org/10.1006/jmbi.1999.3300

# Instructions:

# open PyMOL
# drag this file onto the PyMOL viewing window, which will run all the command below.
# this script shows the wt protein and should show no steric clashes.

# To view the A129F protein, remove the '## ' from two commands below (one starting with fetch, the other with util) and place the '## ' in front of the corresponding commands for the wt protein. Then rerun the script. It should show a tiny bit more steric clashing.

# clear any PyMOL settings and set working directory
reinitialize
cd ~/Desktop

# the python definition below defines the command to show bumping interactions
python
'''
http://pymolwiki.org/index.php/show_bumps
(c) 2011 Thomas Holder, MPI for Developmental Biology
License: BSD-2-Clause
'''
from pymol import cmd

def show_bumps(selection='(all)', name='bump_check', quiet=1):
    '''
DESCRIPTION
    Visualize VDW clashes
ARGUMENTS

    selection = string: atom selection {default: all}

    name = string: name of CGO object to create {default: bump_check}
    '''
    cmd.delete(name)
    cmd.create(name, selection, zoom=0)
    cmd.set('sculpt_vdw_vis_mode', 1, name)
    cmd.set('sculpt_field_mask', 0x020)  # cSculptVDW
    for state in range(1, 1 + cmd.count_states('%' + name)):
        cmd.sculpt_activate(name, state)
        strain = cmd.sculpt_iterate(name, state, cycles=0)
        if not int(quiet):
            print('VDW Strain in state %d: %f' % (state, strain))
    cmd.show_as('cgo', name)

cmd.extend('show_bumps', show_bumps)
python end

# download the wt structure
fetch 1l63, type = pdb, name = 'WT', async = 0

# download the A129F variant structure
## fetch 1QTC, type = pdb, name = 'A129F', async = 0

# remove water molecules for clarity
# show side chains as sticks colored by atom
remove resn HOH
set cartoon_side_chain_helper, on
show sticks, all

# use green carbons for wt
util.cbag wt

# use purple carbons for A129F
## util.cbab A129F

# show steric clashes as disks between atoms
# creates an object named 'bump_check'
show_bumps

# color residue 129 orange
util.cbao resi 96

# zoom to residue 129
set_view (\
     0.753239810,    0.359877706,   -0.550557375,\
    -0.643795133,    0.231922507,   -0.729203105,\
    -0.134736434,    0.903710246,    0.406377584,\
     0.000000000,    0.000000000,  -66.521865845,\
    29.566997528,    6.366400242,   -6.990599632,\
    52.446342468,   80.597389221,  -20.000000000 )

# you should see no bumping interactions (red disks) near Ala129 for wt,
# but slightly more bumping interactions near Phe129 for A129F
