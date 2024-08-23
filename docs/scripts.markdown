---
layout: page
title: Scripts
permalink: /scripts/
---
#### [Check internal cavities]({% link assets/check_internal_cavities.pml %})

This script identifies an internal pocket created by the T4 lysozyme L99A substitution by comparing the wild type structure (pdb code: 1l63) and the L99A variant (pdb code: 1l90). 

Eriksson _et al._ (1992) Response of a protein structure to cavity-creating mutations and its relation to the hydrophobic effect. _Science_, **255**, 178.

#### [Check steric clashes]({% link assets/check_steric_clashes.pml %})

This script identifies a small steric clash created by the T4 lysozyme A129F
substitution. First, it shows the wt protein (pdb code: 1l63) with no steric clash near Ala129. Next, modify the script to show A129F (pdb code: 1QTC), which shows a little steric clash near Phe129, but less than might be expected because the surrounding region
adjusts for the addition of the large phenyl group.

Liu _et al._ (2000) The introduction of strain and its effects on the structure and stability of T4 lysozyme. _J. Mol. Biol._, **295**, 127.

#### [Check surface residues]({% link assets/check_surface_residues.pml %})

This script identifies surface residues, specifically methionines on an Fc fragment of an antibody. These methionines are readily oxidized and the oxidized form is less stable. This script can also be modified to show the aromatic residues on the surface.

Liu _et al._ (2008) Structure and stability changes of human IgG1 Fc as a consequence of methionine oxidation. _Biochemistry_, **47**, 5088.
