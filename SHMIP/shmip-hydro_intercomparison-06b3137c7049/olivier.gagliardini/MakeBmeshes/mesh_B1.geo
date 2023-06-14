L = 100.0e3 ; 
W = 20.0e3 ; 
Mesh.Algorithm = 2 ; 
lc =  500.0 ; 
Point(1) = {0.0,0.0,0.0,lc}; 
Point(2) = { L ,0.0,0.0,lc}; 
Translate {0.0, W , 0.0} {Duplicata{Point{2} ; } }
Translate {0.0, W , 0.0} {Duplicata{Point{1} ; } }
Point(5)={ 5.9000000e+04, 8.0000000e+03,0.0,lc}; 
Line(1) = {1,2} ; 
Line(2) = {2,3} ; 
Line(3) = {4,3} ; 
Line(4) = {4,1} ; 
Line Loop(5) = {1,2,-3,4}; 
Plane Surface(10) = {5}; 
Point{5} In Surface{10}; 
Physical Point(1) = {5}; 
Physical Line(2) = {1}; 
Physical Line(3) = {2}; 
Physical Line(4) = {3}; 
Physical Line(5) = {4}; 
Physical Surface(6) = {10}; 
