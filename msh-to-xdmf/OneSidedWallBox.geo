// Gmsh project created on Wed Jun  4 16:19:45 2025
SetFactory("OpenCASCADE");
//+
//+
Box(1) = {0, 0, 0, 1, 1, 5};
//+
Physical Surface(1) = {1,2};
Physical Surface(2) = {1,2,3,4,5,6};
Physical Volume(1) = {1};
//+
MeshSize { PointsOf{ Volume{1}; } } = 0.2;