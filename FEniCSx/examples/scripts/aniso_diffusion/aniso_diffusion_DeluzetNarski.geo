// Gmsh unit square [0,1]x[0,1]
Mesh.MshFileVersion = 2.2;
h = 0.1;  //refinement factor h
Point(1) = {0.0,0.0, 0.0, h};
Point(2) = {0.0, 1.0, 0.0, h};
Point(3) = {1.0, 1.0, 0.0, h};
Point(4) = {1.0, 0.0, 0.0, h};
Line(7) = {1,2};
Line(8) = {2,3};
Line(9) = {3,4};
Line(10) = {4,1};
Line Loop(5) = {7,8,9,10};
Plane Surface(6) = {5};
Physical Curve(14) = {7};
Physical Curve(13) = {8};
Physical Curve(12) = {9};
Physical Curve(11) = {10};
Physical Surface(15) ={6};
