// Non-conformal unit square [0,1]x[0,1]
Mesh.MshFileVersion = 2.2;

h = 1.1;
delta = 0.008730192;
Point(1)  = {1.0, 0.0, 0.0, h};
Point(2)  = {0.5, 0.0, 0.0, h};
Point(3)  = {1.0, 0.5+delta, 0.0, h};
Point(4)  = {0.5, 0.5-delta, 0.0, h};
Point(5)  = {1.0, 1.0, 0.0, h};
Point(6)  = {0.5, 1.0, 0.0, h};
Point(7)  = {0.5, 0.0, 0.0, h};
Point(8)  = {0.0, 0.0, 0.0, h};
Point(9)  = {0.5, 0.5+delta, 0.0, h};
Point(10) = {0.0, 0.5-delta, 0.0, h};
Point(11) = {0.5, 1.0, 0.0, h};
Point(12) = {0.0, 1.0, 0.0, h};
Line(1) = {1,2};
Line(2) = {3,4};
Line(3) = {1,3};
Line(4) = {2,4};
Line(5) = {5,6};
Line(6) = {3,5};
Line(7) = {4,6};
Line(8) = {7,8};
Line(9) = {9,10};
Line(10) = {7,9};
Line(11) = {8,10};
Line(12) = {11,12};
Line(13) = {9,11};
Line(14) = {10,12};

Line Loop(15) = {11,-9, -10, 8};
Line Loop(16) = {14,-12, -13, 9};
Line Loop(17) = {7,-5,-6, 2};
Line Loop(18) = {4,-2, -3, 1};

Mesh.RecombineAll = 1;
// Mesh.CharacteristicLengthFactor=1.0;

Plane Surface(1) = {15};
Plane Surface(2) = {16};
Plane Surface(3) = {17};
Plane Surface(4) = {18};

Physical Line(1) = {11,14};
Physical Line(2) = {12,5};
Physical Line(3) = {6,3};
Physical Line(4) = {1,8};
Physical Surface(1) ={1,2,3,4};