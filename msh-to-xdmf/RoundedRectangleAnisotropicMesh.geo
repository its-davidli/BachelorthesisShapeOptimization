// Gmsh project created on Wed May 28 18:32:26 2025
SetFactory("OpenCASCADE");
//+
Rectangle(1) = {-1, -0.5, 0, 2, 1, 0.2};
//+
Physical Curve(9) = {6, 7, 8, 1, 2, 3, 4, 5};
//+
Physical Surface(10) = {1};

lc = 0.05;
Point(100) = {0, 0, 0, lc};

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
Field[2].SizeMin = lc / 5;
Field[2].SizeMax = lc;
Field[2].DistMin = 0.2;
Field[2].DistMax = 0.3;

Background Field = 2;



