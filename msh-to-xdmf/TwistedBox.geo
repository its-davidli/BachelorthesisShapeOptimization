// Gmsh project created on Wed Jun  4 16:19:45 2025
//SetFactory("OpenCASCADE");
//+
//+
lc = 0.15;
Point(1) = {0,0,0,lc};
Point(2) = {1,0,0,lc};
Point(3) = {1,1,0,lc};
Point(4) = {0,1,0,lc};
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
Line Loop(1) = {1,2,3,4};
Plane Surface(1) = {1};
Extrude {{0,0,5},{0, 0, 1}, {0.5, 0.5, 0}, Pi/4} {
  Surface{1}; 
}

Physical Surface(1) = {17, 25};
Physical Surface(2) = {1, 13, 17, 21, 25, 26};
Physical Volume(1) = {1};
