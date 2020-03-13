// Gmsh project created on Mon Mar  9 17:23:50 2020
SetFactory("OpenCASCADE");
//+
Rectangle(1) = {0, 0, 0, 0.1, 0.1, 0};
//+
Point(5) = {0.05, 0.05, 0, 1.0};
//+
Physical Point("inlet") = {5};
//+
Physical Point("outlet") = {4, 1, 2, 3};
//+
Physical Surface("solid") = {1};
//+
Point(6) = {0.05, -0, -0, 1.0};
//+
Point(7) = {0.1, 0.05, 0, 1.0};
//+
Point(8) = {0.05, 0.1, 0, 1.0};
//+
Point(9) = {0, 0.05, -0, 1.0};
//+
Line(5) = {6, 5};
//+
Line(6) = {5, 9};
//+
Line(7) = {5, 8};
//+
Line(8) = {5, 7};
//+
Recombine Surface {1};
