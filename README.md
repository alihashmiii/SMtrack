Single Particle Tracking 
=======

single molecule tracking algorithm based upon a modified version of Lineage Mapper. The package implements a robust single molecule tracking scheme. The procedure is somewhat similar to the underlying algorithm proposed in Lineage Mapper (Chalfoun et al, Sci. Reports 2016). A novel implementation of a per-frame particle jump distance has been incorporated, which relies on a mean distance(computed via a
Delaunay Mesh) and its minimization such that a maximum of 0 or 1 particle association is obtained between consecutive frames. The
implementation is expected to work successfully with numerous particle detection strategies. Another novel aspect is the inclusion of " subpixel particle localization" which can be achieved using a two-dimensional Gaussian-Fit

