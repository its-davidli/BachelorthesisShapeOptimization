// Gmsh project created on Wed Jun  4 16:19:45 2025
SetFactory("OpenCASCADE");
//+
Ellipse(1) = {0, 0, 0, 2, 1, 0, Pi};
Line(2) = {2,1};
Extrude {{1, 0, 0}, {1, 0, 0}, 2*Pi} {Curve{1}; }
//Curve Loop(5) = {1,2,3};

Surface Loop(1) = {1};
Volume(1) = {1};

Physical Surface(1) = {2};
Physical Volume(1) = {1};

MeshSize { PointsOf{ Surface{1}; } } = 0.2;

