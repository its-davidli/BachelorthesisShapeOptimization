// Gmsh project created on Wed Jun  4 16:19:45 2025
SetFactory("OpenCASCADE");
//+
Ellipse(1) = {0, 0, 0, 4, 1, 0, 2*Pi};
//+
Curve Loop(1) = {1};
//+
Surface(1) = {1};
//+
Plane Surface(2) = {1};
//+
Physical Curve(10) = {1};
//+
Physical Surface(11) = {1};
//+
MeshSize { PointsOf{ Surface{1}; } } = 0.05;

