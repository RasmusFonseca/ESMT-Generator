# ESMT-Generator
Generate problem instances for the Euclidean Steiner Minimal Tree problem. 

# Usage:

```
./STPGenerator <pointset type> [options]

Pointset type must be one of:
	cube : points distributed randomly in a unit cube
	sphere : points distributed randomly in a unit sphere
	edge : points distributed randomly on diagonal of a hypercube
	spheresurface : points distributed evenly on a unit sphere surface
	sausage : the corners of a sequence of regular simplices that are face-to-face
	grid : the vertices of a cubic grid of width floor(n^(1/d))

Options can be any list of:
	-n <int> : number of points (default=10)
	-d <int> : dimension of points (default=2)
	-s <int> : seed for generating random points (default=time())
	-name <string> : name of point set (default="")
```

For example, running `./STPGenerator spheresurface -n 5 -d 2 -name Circle` will output the following to stdout:
```
33D32945 STP File, STP Format Version 1.0

SECTION Comments
Name "Circle"
Problem "Euclidean Steiner Tree Problem"
Remark "Generated with STPGenerator. 2 dimensions."
END

SECTION Graph
Nodes 5
END

SECTION Coordinates
DD 0.2394268165 0.9709144141
DD -0.9850239519 0.1724175578
DD -0.4381797079 0.8988873921
DD 0.1382872742 -0.9903921596
DD -0.6461654211 -0.7631973851
END
```
