// Gmsh project created on Wed Jun  4 16:19:45 2025
SetFactory("OpenCASCADE");
//+
//+
Box(1) = {-0.1, -0.1, -0.2, 0.2, 0.2, 0.4};
//+
Sphere(2) = {0, 0, 0.02, 0.025, -Pi/2, Pi/2, 2*Pi};
Physical Surface(1) = {7};
Physical Surface(2) = {5,6};
Physical Surface(3) = {1,2,3,4};
Volume(3) = {1,2};
Physical Volume(1) = {3};
Delete { Volume{1,2}; }
//MeshSize { PointsOf{ Volume{3}; } } = 0.01;

lc = 0.05;
Point(100) = {0, 0, 0.02, lc};

// Say we would like to obtain mesh elements with size lc/30 near curve 2 and
// point 5, and size lc elsewhere. To achieve this, we can use two fields:
// "Distance", and "Threshold". We first define a Distance field (`Field[1]') on
// points 100. This field returns the distance to point 100.
Field[1] = Distance;
Field[1].PointsList = {100};


// We then define a `Threshold' field, which uses the return value of the
// `Distance' field 1 in order to define a simple change in element size
// depending on the computed distances
//
// SizeMax -                     /------------------
//                              /
//                             /
//                            /
// SizeMin -o----------------/
//          |                |    |
//        Point         DistMin  DistMax
Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = lc / 7;
Field[2].SizeMax = lc;
Field[2].DistMin = 0.075;
Field[2].DistMax = 0.15;

Background Field = 2;

