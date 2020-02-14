[p e t] = gmsh2pdetoolbox("quad1500.msh", 2, 3);

fem = heatFEM(p, e, t, 2)

fem.
