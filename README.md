# STPGenerator
Generate problem instances for the Euclidean Steiner Tree Problem. 

Usage:

    ./STPGenerator <pointset type> [options]
    
    Pointset type must be one of:
    	cube : points distributed randomly in a unit cube
    	sphere : points distributed randomly in a unit sphere
    	edge : points distributed randomly on diagonal of a hypercube
    	spheresurface : points distributed evenly on a unit sphere surface
    	sausage : the corners of a sequence of regular simplices that are face-to-face
    	grid : the vertices of a cubic grid of width floor(n^(1/d))
    
    Options can be any list of:
    	-n <int> : number of points (standard=10)
    	-d <int> : dimension of points (standard=2)
    	-s <int> : seed for generating random points (standard=time())
    	-name <string> : name of point set (standard="")
