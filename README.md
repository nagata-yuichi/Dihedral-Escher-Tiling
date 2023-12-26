# Creation of Dihedral Escher-like Tilings Based on As-Rigid-As-Possible Deformation

This code is a C++ implementation of the algorithms developed in the following papers. 

[1] Yuichi Nagata and Shinji Imahori, Creation of Dihedral Escher-like Tilings Based on As-Rigid-As-Possible Deformation, ****, ****, ****, 2024. [[link]](https://dl.acm.org/doi/10.1145/3638048)

Outline: An Escher-like tiling is a tiling consisting of one or a few artistic shapes of tile. The proposed algorithm generates Escher-like tilings consisting of two distinct shapes (dihedral Escher-like tilings) that are as similar as possible to the two goal shapes specified by the user. Within the algorithm, the two tile shapes are optimized to minimize a distance function that measures the difference from the two goal shapes under the constraints that they can tile the plane. The distance between the tile and goal shapes is evaluated based on the as-rigid-as-possible (ARAP) deformation energy. This allows for the generation of satisfactory tile shapes by deforming the goal shapes in physically plausible ways even when large deformations of the goal shapes are indispensable to form possible tile shapes.

<img src="Images/fig.png" width="100%">

# Environment
The source code is implemented in C++ with the Eigen library. The author ran the program on Ubuntu 20.04, where the following preparation may be required to compile the source code. 
- Make the Eigen library [[link]](https://eigen.tuxfamily.org) available on your PC. 
- Install the X11 library on your PC. `$ sudo apt install libx11-dev`

# Directories
- Program: Source codes for generating and displaying tilings
- Data: Files of goal shapes (mesh representation)
- Make_Polygon: A program for creating goal shapes (polygon)
- Make_Mesh: Source codes for constructing goal shapes (mesh)

# Compile
Execute the following in the directory Program.

```
$ chmod 700 build_I_get_conf.exe
$ ./build_I_get_conf.exe
```
Then, an executable file "jikken_I_get_conf" is generated.

```
$ chmod 700 build_I_conf.exe
$ ./build_I_conf.exe
```
Then, an executable file "jikken_I_conf" is generated.

```
$ chmod 700 build_display.exe
$ ./build_display.exe
```
Then, an executable file "display" is generated.

Note:
You may need to modify the code that includes Eigen libray (see error messages). 

# How to perform the EST
(1) Selecting promising template configurations

execution:
```
$ ./jikken_I_get_conf <string1> <string2> <string3> <integer1> <integer2>
```
&nbsp; \<string1\> : file name of the first goal shape (mesh representation)  
&nbsp; \<string2\> : file name of the second goal shape (mesh representation)  
&nbsp; \<string3\> : file name to which results are written  
&nbsp; \<interger1\> : the number of template configurations stored (1000000 in the paper)  
&nbsp; \<interger2\> : the number of threads for parallel execution  

example:
```
$ ./jikken_I_get_conf seahorse_60_36.dat pegasus_60_39.dat AAA 1000000 20
```
Then, the specified number of template configurations are stored in a file (AAA_seahorse_60_36_pegasus_60_39.conf in this example).

(2) Solve the optimization problem for the selected template configurations

execution:
```
$ ./jikken_I_conf <string1> <string2> <string3> <integer1> <double1> <string4> <interger2>
```
&nbsp; \<string1\>: file name of the first goal shape (mesh representation)  
&nbsp; \<string2\>: file name of the second goal shape (mesh representation)  
&nbsp; \<string3\>: file name to which results are written  
&nbsp; \<interger1\>: the number of top solutions stored (10 in the paper)  
&nbsp; \<double1\>: the value of alpha for the inner points (0.5 in the paper)  
&nbsp; \<interger2\>: the number of threads for parallel execution  

example:
```
$ ./jikken_I_conf seahorse_60_36.dat pegasus_60_39.dat BBB 100 0.5 AAA_seahorse_60_36_pegasus_60_39.conf 20    
```
Then, the specified number of top solutions are stored in a file (BBB_seahorse_60_36_pegasus_60_39.tile in this example).
