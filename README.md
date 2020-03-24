# networkMFPT

Matlab code for analytically calculating diffusive mean first passage times (MFPT) on a network

----------------

This code is associated with
Brown, Aidan I., Laura M. Westrate, and Elena F. Koslover. "Impact of global structure on diffusive exploration of organelle networks." Scientific Reports 10.1 (2020): 1-13.

Please cite the article above if using the code.

Any questions should be directed to Elena Koslover (ekoslover@physics.ucsd.edu).

--------------


You need to provide a network structure file. Containing the following.
- lines of format "NODE ID x y" or "NODE ID x y z" for 2D or 3D networks, respectively
- lines of format "EDGE ID N1 N2 len" where ID is an integer for the edge ID, N1 and N2 are IDs of nodes connected by the edge. "len" is the edge length (optional, can be omitted). If no edge length is supplied for a given edge, the straight-line distance between the two points is calculated.

An example network structure file is provided in example.net

Run example.m to calculate the mean first passage time with two target nodes (node 3 and 6). This the the MFPT to hit the first of the targets for a diffusing particle.

The resulting MFPT values for each starting node are given in the array MFPTs.
