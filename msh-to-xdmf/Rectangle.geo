// Gmsh project created on Wed May 28 18:32:26 2025
SetFactory("OpenCASCADE");
//+
Rectangle(1) = {-1, -0.5, 0, 2, 1, 0.1};
//+
Physical Curve(1) = {1, 2, 3, 4,5,6,7,8};
//+
Physical Surface(1) = {1};

MeshSize { PointsOf{ Surface{1}; } } = 0.018;
