// Constructs an (approximately) annular mesh with quad elements
// Loosely based on https://github.com/ExCALIBUR-NEPTUNE/nektar-diffusion/blob/master/example/Torus/domain.geo

// ================================== Params ==================================
// Circle radii
r_inner = 1.0;
r_outer = 2.0;

// Radial resolution
res_r = 10;
// Target azimuthal resolution; actual res will be 3*floor(res_theta/3)
res_theta = 30;

//=============================== gmsh entities ===============================
// N.B. Constructs the mesh in 3 sections of 2*PI/3 each (to satisfy gmsh 'Circle' constraints) then recombines them at the end

// Centre
Point(0) = {0,0,0,1.0};

// Inner circle points at theta = 0, 2*PI/3, 4PI/3
Point(1) = {r_inner,0,0,1.0};
Point(2) = {r_inner*Cos(2*Pi/3),r_inner*Sin(2*Pi/3),0,1.0};
Point(3) = {r_inner*Cos(4*Pi/3),r_inner*Sin(4*Pi/3),0,1.0};

// Outer circle points at theta = 0, 2*PI/3, 4PI/3
Point(11) = {r_outer,0,0,1.0};
Point(12) = {r_outer*Cos(2*Pi/3),r_outer*Sin(2*Pi/3),0,1.0};
Point(13) = {r_outer*Cos(4*Pi/3),r_outer*Sin(4*Pi/3),0,1.0};

// Radials at theta = 0, 2*PI/3, 4PI/3
Line(21) = {1,11};
Line(22) = {2,12};
Line(23) = {3,13};

// Inner boundary segments
Circle(31) = {1, 0, 2};
Circle(32) = {2, 0, 3};
Circle(33) = {3, 0, 1};

// Outer boundary segments
Circle(41) = {11, 0, 12};
Circle(42) = {12, 0, 13};
Circle(43) = {13, 0, 11};

// Define planes for each segment
Line Loop(51) = {21, 41, -22, -31};
Plane Surface(61) = {51};
Line Loop(52) = {22, 42, -23, -32};
Plane Surface(62) = {52};
Line Loop(53) = {23, 43, -21, -33};
Plane Surface(63) = {53};

// Split radial lines
Transfinite Line {21,22,23} = res_r+1;
// Split inner and outer circles
Transfinite Line {31,32,33,41,42,43} = res_theta/3+1;

// Mesh the three sections, then recombine
Transfinite Surface {61} = {1, 11, 12, 2};
Transfinite Surface {62} = {2, 12, 13, 3};
Transfinite Surface {63} = {3, 13, 11, 1};
Recombine Surface {61,62,63};

//=================================== Labels ==================================
// Whole domain
Physical Surface(0) = {61, 62, 63};
// Inner boundary
Physical Line(1) = {31,32,33};
//Physical Line("inner") = {31,32,33};
// Outer boundary
Physical Line(2) = {41,42,43};
//Physical Line("outer") = {41,42,43};