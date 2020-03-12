// Gmsh project created on Mon Mar  9 17:42:13 2020
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {0.5, 0, 0, 1.0};
//+
Point(3) = {1, 0, 0, 1.0};
//+
Point(4) = {0, 0.5, -0, 1.0};
//+
Point(5) = {0, 1, -0, 1.0};
//+
Point(6) = {0.5, 0.5, -0, 1.0};
//+
Point(7) = {1, 0.5, 0, 1.0};
//+
Point(8) = {1, 1, 0, 1.0};
//+
Point(9) = {0.5, 1, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 7};
//+
Line(4) = {7, 8};
//+
Line(5) = {8, 9};
//+
Line(6) = {9, 5};
//+
Line(7) = {5, 4};
//+
Line(8) = {4, 1};
//+
Line(9) = {2, 6};
//+
Line(10) = {6, 4};
//+
Line(11) = {9, 6};
//+
Line(12) = {6, 7};
//+
Curve Loop(1) = {7, 10, 11, 6};
//+
Surface(1) = {1};
//+
Curve Loop(3) = {5, 11, 12, 4};
//+
Surface(2) = {3};
//+
Curve Loop(5) = {8, 1, 9, 10};
//+
Surface(3) = {5};
//+
Curve Loop(7) = {2, 3, -12, -9};
//+
Surface(4) = {7};
//+
Transfinite Surface {1} = {6, 9, 5, 4};
//+
Transfinite Surface {2} = {7, 8, 9, 6};
//+
Transfinite Surface {4} = {3, 7, 6, 2};
//+
Transfinite Surface {3} = {2, 6, 4, 1};
//+
Recombine Surface {1, 3, 4, 2};
//+
Physical Surface("solid") = {1, 3, 4, 2};
//+
Physical Point("inlet") = {6};
//+
Physical Point("outlet") = {5, 1, 3, 8};
