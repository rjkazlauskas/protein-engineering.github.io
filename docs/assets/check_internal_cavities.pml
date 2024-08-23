## Check a protein for internal cavities ##

# This script identifies an internal pocket created by the T4 lysozyme L99A substitution by comparing the wild type structure (pdb code: 1l63) and the L99A variant (pdb code: 1l90). 
# Eriksson _et al._ (1992) Response of a protein structure to cavity-creating mutations and its relation to the hydrophobic effect. _Science_, **255**, 178. https://doi.org/10.1126/science.1553543

# Instructions: 
# open PyMOL
# drag this file onto the PyMOL viewing window

# clear any PyMOL settings and set working directory
reinitialize
cd ~/Desktop

# download the wild-type (1l63) and Leu99Ala (1l90) structures from the protein data bank
# both proteins also contain the Cys54Thr/Cys97Ala substitutions to remove a disulfide 
# crosslink. The stability and activity of this cysteine-free mutant is essentially
# identical with the wild-type enzyme.

fetch 1l63, type = pdb, name = 'WT', async = 0
fetch 1l90, type = pdb, name = 'Leu99Ala', async = 0

# remove water molecules for clarity
# show side chains as sticks colored by atom
# use green carbons for wild-type structure
# use purple carbons for Leu99Ala structure
remove resn HOH
set cartoon_side_chain_helper, on
show sticks, all
util.cbag WT
util.cbab Leu99Ala

# display surface of pockets and cavities only
set surface_cavity_mode, 1

# set transparency of surface to 40% and color to yellow
set transparency, 0.25
set surface_color, yellow
show surface

# show position 99 as orange sticks
util.cbao resi 99

# adjust to zoom in and orient view to position 99 and the associated cavity
set_view (\
     0.280132115,   -0.950214684,   -0.136434376,\
     0.416004360,   -0.007920738,    0.909332395,\
    -0.865141094,   -0.311487436,    0.393072486,\
     0.000000000,    0.000000000,  -74.134986877,\
    30.149997711,    5.736000538,    5.680692673,\
    58.448589325,   89.821388245,  -20.000000000 )
        
# toggle visibility of WT and Leu99Ala to compare cavities
# note the close overlap of the two structures
# the internal surface is dark and hard to see, type the command 'ray' to see it in yellow