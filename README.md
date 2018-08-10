Single Particle Tracking
=======

`Currently Active`

SMTrack package implements a robust single molecule tracking scheme based upon a modified version of Lineage Mapper (Chalfoun et al, Scientific Reports 2016). A novel implementation of a per-frame particle jump-distance has been incorporated, which relies on a mean distance (computed via a Delaunay Mesh) and its minimization such that a maximum of 0 or 1 particle association is obtained between consecutive frames. There is also the option to specify a custom search distance.

The implementation is expected to work successfully with numerous particle detection strategies. Another novel aspect is the inclusion of "subpixel particle localization" which relies on a two-dimensional Gaussian-Fit.

Note: the current detection scheme utilizes a Laplacian of Gaussian (LoG) filter with thresholding. Another variant uses the maxima of the distance transform of the thresholded LoG-filtered image as seeds to WatershedComponents (gradient descent version or rainfall algorithm (Osma-Ruiz)) to segment for individual spots. 2D Gaussians can be fit to either the thresholded image or the watershed components.

`Possible improvements (Near Future):`

###### 1. Add a different detection scheme (example, FogBank by NIST)

###### 2. Not necessary but perhaps important for low density cases: Kalman predictor to handle occlusions when a particle motion-model exists (e.g. constant velocity etc ..) 




##### Input:

![alt-text](https://github.com/alihashmiii/SMtrack/blob/master/for%20readme/input.gif) 


##### for computing per-frame jump distance:

![alt-text](https://github.com/alihashmiii/SMtrack/blob/master/for%20readme/overlayed%20jump%20distance.png)


##### superposition of all particle trajectories for the above movie `subpixel localization = False`:

![alt-text](https://github.com/alihashmiii/SMtrack/blob/master/for%20readme/segment_output.png)


##### a Gaussian (rainbow) is fit to the data (yellow) for achieving subpixel localization:

![alt-text](https://github.com/alihashmiii/SMtrack/blob/master/for%20readme/gaussianFit.png)


##### superposition of all particle trajectories for the above movie `subpixel localization = True`:

![alt-text](https://github.com/alihashmiii/SMtrack/blob/master/for%20readme/subpixoutput.png)

