// Gmsh project created on Sat Apr 11 02:17:06 2020
SetFactory("OpenCASCADE");

//+
Rectangle(1) = {0, 0, 0, 0.05, 0.05, 0};
//+
Extrude {0, 0, 0.05} {
  Surface{1}; 
}
//+
Physical Point("inlet") = {1};
//+
Physical Point("outlet") = {7};
//+
Physical Volume("solid") = {1};
//+
Transfinite Volume{1} = {1, 5, 8, 4, 2, 6, 7, 3};
//+
Transfinite Surface {5} = {1, 5, 8, 4};
//+
Transfinite Surface {6} = {5, 6, 7, 8};
//+
Transfinite Surface {4} = {4, 8, 7, 3};
//+
Transfinite Surface {1} = {4, 3, 2, 1};
//+
Transfinite Surface {3} = {3, 7, 6, 2};
//+
Transfinite Surface {2} = {2, 6, 5, 1};
