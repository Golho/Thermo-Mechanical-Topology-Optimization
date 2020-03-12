gmsh = gmshParser('meshes/long_quad_coarse.msh');
timeSteps = 5;
tFinal = 20;

% Create boundary conditions
prescribed = struct(...
    'physicalName', 'outlet', ...
    'type', 'Dirichlet', ...
    'value', 0, ...
    'timeSteps', 1:timeSteps ...
);

flux = struct(...
    'physicalName', 'inlet', ...
    'type', 'Neumann', ...
    'value', 10, ...
    'timeSteps', 1:timeSteps ...
);

% Create body conditions
body = struct(...
    'type', 'main', ...
    'physicalName', 'solid' ...
);

fem = OptHeatFEM(gmsh, tFinal, timeSteps);
fem.addBoundaryCondition(prescribed);
fem.addBoundaryCondition(flux);
fem.addBodyCondition(body);

topOpt = HeatComplianceProblem(fem, Elements.QUA_4);
designPar = 0.5*ones(size(fem.mainDensities));
g(1) = topOpt.objective(designPar)

[Ex, Ey] = fem.getElemTemp(0);
figure
elfield2(Ex, Ey, designPar);

% Numerical gradient
dgdphi = numGrad(@topOpt.objective, designPar, 1e-6)

% Analytical gradient
der_g = topOpt.gradObjective(designPar)