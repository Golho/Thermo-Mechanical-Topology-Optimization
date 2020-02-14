// Gmsh project created on Wed Feb 19 18:39:42 2020
SetFactory("OpenCASCADE");
//+
Circle(1) = {0, -0, 0, 0.1, e, 2*Pi};
//+
Rectangle(1) = {0, 0, 0, 0.05, 0.05, e};
//+
Curve Loop(2) = {1};
//+
Curve Loop(3) = {5, 2, 3, 4};
//+
Plane Surface(2) = {2, 3};
//+
Physical Surface("circle with rectangle") = {2};
