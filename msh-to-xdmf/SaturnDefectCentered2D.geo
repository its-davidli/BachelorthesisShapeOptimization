// Gmsh project created on Wed Jun  4 16:19:45 2025
SetFactory("OpenCASCADE");
//+
Rectangle(1) = {-0.5, -0.5, 0, 1.0, 1.0, 0};
//+
Circle(5) = {0, 0.025, 0, 0.025, 0, 2*Pi};
Line Loop(6) = {5};
Plane Surface(2) = {6};

Physical Curve(1) = {1, 4, 5};
Physical Curve(2) = {1, 2, 3, 4};
Physical Curve(3) = {5};

BooleanDifference(3) = { Surface{1}; Delete; }{ Surface{2}; Delete; };
Physical Surface(1) = {3};

lc = 0.05;
Point(100) = {0, 0.025, 0, lc};

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
Field[2].SizeMin = lc / 10;
Field[2].SizeMax = lc;
Field[2].DistMin = 2;
Field[2].DistMax = 2.3;

Background Field = 2;

//+
