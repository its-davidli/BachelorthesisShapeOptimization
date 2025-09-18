// Gmsh project created on Tue Aug  5 17:39:13 2025
SetFactory("OpenCASCADE");
//+

lc = 0.5;
Rectangle(1) = {-5, -5, 0, 10, 10, 0};
Delete {
    Surface{1};
}
//+
Physical Curve(1) = {1,2,3,4};
//+
Circle(5) = {0, 0, 0, 1, 0, 2*Pi};
//+
Physical Curve(2) = {5};

Line Loop(3) = {1,2,3,4};
//+
Line Loop(4) = {5};
//+
Plane Surface(2) = {3, 4};
//Plane Surface(3) = {4};
//+
Physical Surface(1) = {2};
//Physical Surface(2) = {3};
Mesh.MeshSizeMin = 0.1;
Mesh.MeshSizeMax = 0.1;

//Point(100) = {0, 0, 0, lc};
// Say we would like to obtain mesh elements with size lc/30 near curve 2 and
// point 5, and size lc elsewhere. To achieve this, we can use two fields:
// "Distance", and "Threshold". We first define a Distance field (`Field[1]') on
// points 5 and on curve 2. This field returns the distance to point 5 and to
// (100 equidistant points on) curve 2.
//Field[1] = Distance;
//Field[1].SurfacesList = {3};
//Field[1].Sampling = 100;


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
//Field[2] = Threshold;
//Field[2].InField = 1;
//Field[2].SizeMin = lc / 15;
//Field[2].SizeMax = lc;
//Field[2].DistMin = 1;
//Field[2].DistMax = 3;
//
//Background Field = 2;
//
//Mesh.MeshSizeFromPoints = 0;
//Mesh.MeshSizeFromCurvature = 0;
//Mesh.MeshSizeExtendFromBoundary = 0;
//Mesh.Algorithm = 5;

