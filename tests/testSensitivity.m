clear; close all;
%%
gmsh = gmshParser('meshes/square_coarse.msh');
timeSteps = 50;
tFinal = 1000;

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
    'value', 0.1/0.001, ...
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

fem.setMaterial(struct('D', 10*eye(2), 'density', 1, 'heatCapacity', 1e6), 0.001);

material_1 = struct('kappa', 0.1, 'cp', 5e5);
material_2 = struct('kappa', 10, 'cp',  1e6);

options = struct(...
    'p_kappa', 3, ...
    'p_cp', 3, ...
    'filter', true, ...
    'filterRadius', 0.02, ...
    'material_1', material_1, ...
    'material_2', material_2 ...
);

%%
topOpt = HeatComplianceProblem(fem, Elements.QUA_4, options, 0.4);
designPar = 0.5*ones(size(fem.mainDensities));
g(1) = topOpt.objective(designPar)

[Ex, Ey] = fem.getElemTemp(0);
figure
elfield2(Ex, Ey, designPar);

% Numerical gradient
dgdphi = numGrad(@topOpt.objective, designPar, 1e-6)

% Analytical gradient
der_g = topOpt.gradObjective(designPar)

%%
topOpt = MaxTemperatureProblem(fem, Elements.QUA_4, options, 0.4);
designPar = 0.5*ones(size(fem.mainDensities));
topOpt.normalize(designPar);
g(1) = topOpt.objective(designPar)

[Ex, Ey, ed] = fem.getElemTemp(fem.timeSteps-1);
figure
elfield2(Ex, Ey, ed);

% Numerical gradient
dgdphi = numGrad(@topOpt.objective, designPar, 1e-7)

% Analytical gradient
der_g = topOpt.gradObjective(designPar)