# Creation of Dihedral Escher-like Tilings Based on As-Rigid-As-Possible Deformation

This code is a C++ implementation of the algorithms developed in the following papers. The source code will be available after the paper is published. 

[1] Yuichi Nagata and Shinji Imahori, Creation of Dihedral Escher-like Tilings Based on As-Rigid-As-Possible Deformation, ****, ****, ****, 2024.

Outline: An Escher-like tiling is a tiling consisting of one or a few artistic shapes of tile. The proposed algorithm generates Escher-like tilings consisting of two distinct shapes (dihedral Escher-like tilings) that are as similar as possible to the two goal shapes specified by the user (see the figure below). Within the algorithm, the two tile shapes are optimized to minimize a distance function that measures the difference from the two goal shapes under the constraints that they can tile the plane. The distance between the tile and goal shapes is evaluated based on the as-rigid-as-possible (ARAP) deformation energy. This allows for the generation of satisfactory tile shapes by deforming the goal shapes in physically plausible ways even when large deformations of the goal shapes are indispensable to form possible tile shapes.

<img src="Images/fig.png" width="90%">
