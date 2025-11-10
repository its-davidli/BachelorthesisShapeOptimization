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
//MeshSize { PointsOf{ Surface{1,2,3,4,5,6}; } } = 0.1;
MeshSize { PointsOf{ Volume{1}; } } = 0.15;

Point(100) = {0.5, 0.5, 2.5, 0.1};

// distance to the center
Field[1] = Distance;
Field[1].PointsList = {100};

// small elements at center, grading to 0.1 near the surfaces
Field[2] = Threshold;
Field[2].IField = 1;
Field[2].LcMin = 0.7;
Field[2].LcMax = 0.1;
Field[2].DistMin = 0.0;
Field[2].DistMax = 0.7;

Background Field = 2;