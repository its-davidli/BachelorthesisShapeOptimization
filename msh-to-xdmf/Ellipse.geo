// Gmsh project created on Wed Jun  4 16:19:45 2025
SetFactory("OpenCASCADE");
//+
Ellipse(1) = {0, 0, 0, 2, 1, 0, 2*Pi};
//+
Curve Loop(1) = {1};
//+
Curve Loop(2) = {1};
//+
Curve Loop(3) = {1};
//+
Curve Loop(4) = {1};
//+
Curve Loop(5) = {1};
//+
Curve Loop(6) = {1};
//+
Curve Loop(7) = {1};
//+
Surface(1) = {7};
//+
Curve Loop(9) = {1};
//+
Plane Surface(2) = {9};
//+
Physical Curve(10) = {1};
//+
Physical Surface(11) = {1};
//+
MeshSize { PointsOf{ Surface{1}; } } = 0.03;
