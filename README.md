Single Particle Tracking
=======

`Status: In progress`

SMTrack package implements a robust single molecule tracking scheme based upon a modified version of Lineage Mapper (Chalfoun et al, Scientific Reports 2016). A novel implementation of a per-frame particle jump-distance has been incorporated, which relies on a mean distance (computed via a Delaunay Mesh) and its minimization such that a maximum of 0 or 1 particle association is obtained between consecutive frames.

The implementation is expected to work successfully with numerous particle detection strategies. Another novel aspect is the inclusion of
"subpixel particle localization" which relies on a two-dimensional Gaussian-Fit.

`To Do:`

###### 0. Rework on the 2D Gaussian-Fits (high priority and important for subpixel localization)


`Possible improvements (Near Future):`

###### 1. Add a better particle detection scheme (example, FogBank by NIST) ... can considerably improve the results 

###### 2. Not necessary but perhaps important for low density cases: Kalman predictor to handle occlusions when a particle motion-model exists (e.g. constant velocity etc ..) 




##### Input:

![alt-text](https://github.com/alihashmiii/SMtrack/blob/master/for%20readme/input.gif) 


##### for computing per-frame jump distance:

![alt-text](https://github.com/alihashmiii/SMtrack/blob/master/for%20readme/overlayed%20jump%20distance.png)


##### superposition of all particle trajectories for the above movie (subpixel localization = False):

![alt-text](https://github.com/alihashmiii/SMtrack/blob/master/for%20readme/segoutput.png)


##### a Gaussian is fit to the data for achieving subpixel localization:

![alt-text](https://github.com/alihashmiii/SMtrack/blob/master/for%20readme/gaussianFit.png)


##### superposition of all particle trajectories for the above movie (subpixel localization = True):

![alt-text](https://github.com/alihashmiii/SMtrack/blob/master/for%20readme/subpixoutput.png)

