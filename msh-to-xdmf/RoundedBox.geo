// Gmsh project created on Wed Jun  4 16:19:45 2025
//SetFactory("OpenCASCADE");
//+
//+
lc = 0.1;
Point(1) = {0.1, 0, 0, lc};
Point(2) = {0.9, 0, 0, lc};
Point(3) = {1, 0.1, 0, lc};
Point(4) = {1, 0.9, 0, lc};
Point(5) = {0.9, 1, 0, lc};
Point(6) = {0.1, 1, 0, lc};
Point(7) = {0, 0.9, 0, lc};
Point(8) = {0, 0.1, 0, lc};
// Center Points Arcs
Point(9) = {0.1, 0.1, 0, lc};
Point(10) = {0.9, 0.1, 0, lc};
Point(11) = {0.9, 0.9,0,lc};
Point(12) = {0.1, 0.9, 0,lc};

// Corner arcs (radius 0.1)
Circle(1) = {8, 9, 1};
Circle(2) = {2, 10, 3};
Circle(3) = {4, 11, 5};
Circle(4) = {6, 12, 7};

// Straight edges
Line(5) = {1,2};
Line(6) = {3,4};
Line(7) = {5,6};
Line(8) = {7,8};

// Line loop
Line Loop(1) = {1,5,2,6,3,7,4,8};
Plane Surface(1) = {1};
Extrude {{0,0,5},{0, 0, 1}, {0.5, 0.5, 0}, 0} {
  Surface{1}; 
}

Physical Surface(1) = {49, 33};
Physical Surface(2) = {21,25,29,33,37,41,45,49,50};
Physical Volume(1) = {1};

