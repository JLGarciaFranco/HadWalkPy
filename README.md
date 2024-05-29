This is a repository for code to compute the divergent and non-divergent wind components using available horizontal U and V wind components. 
The main features of the code are:
  1. Computation of the irrotational and non-divergent components of the wind (Helmholtz decomposition)
  2. Local Hadley-Walker circulation characterized by streamfunction and vertical velocities. 

The code is built into a single class under the file: decompose_wind.py. 

An example code with comments is found under the file: example.py

The decomposition is based on the study of Schwendike et al., (2014), which decomposes a full 3D circulation into 2 orthogonal 2D circulations, which is also called the $\psi$ method originally proposed by Keyser et al., (1989).

*References*

Keyser, D., B. D. Schmidt, and D. G. Duffy (1989), A technique for representing three-dimensional vertical circulations in baroclinic disturbance, Mon. Weather Rev., 117, 2463–2494. 
Schwendike, J., Govekar, P., Reeder, M. J., Wardle, R., Berry, G. J., & Jakob, C. (2014). Local partitioning of the overturning circulation in the tropics and the connection to the Hadley and Walker circulations. Journal of Geophysical Research: Atmospheres, 119(3), 1322-1339.
