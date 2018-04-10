Pairwise Distance Pipeline
Built: April 10, 2018
Author: Braden Katzman
Language: MATLAB
----------------------------------

This is a small pipeline to calculate the "neighbors" of the pharyngeal mass at a given time point.

Neighbors are defined as nuclei that are within a certain distance, d = , of the pharynx.

The given time point chosen to calculate the neighbors is t = 300. This is approximately 30 minutes after
the bilaterally symmetric, "two-sheet" stage is recognizable. After t = 300, the pharyngeal mass starts to 
compress and form a spherical structure in the anterior, dorsal end of the embryo.

The goal is to identify those cells that are "neighbors" of the pharynx at the end of the "two-sheet" stage
and identify which of these then migrate away from the pharynx toward the embryo shell.