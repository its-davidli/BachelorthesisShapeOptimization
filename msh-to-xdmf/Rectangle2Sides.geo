// Gmsh project created on Wed May 28 18:32:26 2025
SetFactory("OpenCASCADE");
//+

//+
Point(1) = {-2, -1.75, 0, 1.0};
//+
Point(2) = {2, -1.75, 0, 1.0};
//+
Point(3) = {2, 1.75, 0, 1.0};
//+
Point(4) = {-2, 1.75, 0, 1.0};
//+
Point(5) = {0, 1.8, 0, 1.0};
Point(6) = {0, -1.8, 0, 1.0};

Point(11) = {-2, -1.25,0,1.0};
Point(12) = {2, -1.25,0,1.0};
Point(13) = {2, 1.25,0,1.0};
Point(14) = {-2, 1.25,0,1.0};
Line(1) = {11,14};
BSpline(2) = {14,4,5,3, 13};
Line(3) = {13,12};
BSpline(4) = {12,2,6,1,11};

Curve Loop(1) = {1,2,3,4};
Plane Surface(1) = {1};

Physical Curve(1) = {2,4};
Physical Curve(2) = {1,3};
Physical Surface(1) = {1};
Mesh.MeshSizeMin = 0.03;
Mesh.MeshSizeMax = 0.03;