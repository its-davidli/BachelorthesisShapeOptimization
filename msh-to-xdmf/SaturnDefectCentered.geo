// Gmsh project created on Wed Jun  4 16:19:45 2025
SetFactory("OpenCASCADE");
//+
//+
Box(1) = {-0.5, -0.5, -0.5, 1.0, 1.0, 1.0};
//+
Sphere(2) = {0, 0, 0, 0.05, -Pi/2, Pi/2, 2*Pi};
Physical Surface(1) = {1,2,7};
Physical Surface(2) = {1,2,3,4,5,6};
Physical Surface(3) = {7};
Physical Surface(4) = {1,2}; 
Volume(3) = {1,2};
Physical Volume(1) = {3};
Delete { Volume{1,2}; }
MeshSize { PointsOf{ Volume{3}; } } = 0.075;


