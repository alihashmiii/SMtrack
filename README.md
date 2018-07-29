Single Particle Tracking
=======

The package implements a robust single molecule tracking scheme based upon a modified version of Lineage Mapper (Chalfoun et al, Scientific Reports 2016). A novel implementation of a per-frame particle jump-distance has been incorporated, which relies on a mean distance (computed via a Delaunay Mesh) and its minimization such that a maximum of 0 or 1 particle association is obtained between consecutive frames.

The implementation is expected to work successfully with numerous particle detection strategies. Another novel aspect is the inclusion of
"subpixel particle localization" which relies on a two-dimensional Gaussian-Fit.

`Possible improvements:`

###### 1. add a better particle detection scheme



##### superposition of all particle trajectories for one video:

![alt-text](https://github.com/alihashmiii/SMtrack/blob/master/for%20readme/input.gif) ![alt-text](https://github.com/alihashmiii/SMtrack/blob/master/for%20readme/particle%20Trajectories.png)
