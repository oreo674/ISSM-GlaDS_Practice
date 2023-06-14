// gmsh mesh  
// L, W                                                   
L = 100.0e3 ;               
W = 20.0e3 ; 


// The anisotropic 2D mesh generator can be selected with:

Mesh.Algorithm = 2 ;

// One can force a 4 step Laplacian smoothing of the mesh with:

//Mesh.Smoothing = 4 ;

lc =  500.0 ;

Point(1) = {0.0,0.0,0.0,lc};
Point(2) = { L ,0.0,0.0,lc};
Point(3) = { L , W ,0.0,lc};
Point(4) = {0.0, W ,0.0,lc};

Line(1) = {1,2} ;
Line(2) = {2,3} ;
Line(3) = {3,4} ;
Line(4) = {4,1} ;
 
Line Loop(5) = {1,2,3,4};

Plane Surface(10) = {5};

Physical Line(1) = {1} ;
Physical Line(2) = {2} ;
Physical Line(3) = {3} ;
Physical Line(4) = {4} ;

Physical Surface(5) = {10} ;


