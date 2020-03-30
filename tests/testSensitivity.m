clear; close all;
timeSteps = 5;
tFinal = 1000;
volumeFraction = 0.4;
radius = 0.025;

material_1 = Material(1, 1e6, 0.1*eye(3));
material_2 = Material(1, 1e6, 10*eye(3));
%% Structured mesh
isnear = @(x, a) abs(x-a) < 1e-3; 
mesh = StructuredMesh([5, 0.1], [5, 0.1]);
globalCoord = mesh.coordinates();

centerNodes = find( (globalCoord(1, :) >= 0.046 & globalCoord(1, :) <= 0.054) & ...
                    (globalCoord(2, :) >= 0.046 & globalCoord(2, :) <= 0.054));
                
cornerNodes = find( (isnear(globalCoord(1, :), 0) | isnear(globalCoord(1, :), 0.1)) & ...
                    (isnear(globalCoord(2, :), 0) | isnear(globalCoord(2, :), 0.1)));

% Create boundary conditions
prescribed = struct(...
    'nodes', cornerNodes, ...
    'type', 'Dirichlet', ...
    'value', 0, ...
    'timeSteps', 1:timeSteps ...
);

fluxCorner = struct(...
    'nodes', centerNodes, ...
    'type', 'Neumann', ...
    'value', 1/length(centerNodes), ...
    'timeSteps', 1:timeSteps ...
);

% Create body conditions
body = struct(...
    'type', 'main' ...
);

fem = OptHeatFEMStructured(mesh, tFinal, timeSteps);

fem.addBoundaryCondition(prescribed);
fem.addBoundaryCondition(fluxCorner);
fem.addBodyCondition(body);

fem.setMaterial(material_2);

options = struct(...
    'p_kappa', 3, ...
    'p_cp', 3, ...
    'filter', true, ...
    'filterRadius', radius, ...
    'filterWeightFunction', @(dx, dy, dz) max(radius-sqrt(dx.^2+dy.^2+dz.^2), 0), ...
    'material_1', material_1, ...
    'material_2', material_2 ...
);
%%
topOpt = HeatComplianceProblem(fem, options, 0.4);
designPar = 0.5*ones(size(fem.designPar));
g(1) = topOpt.objective(designPar)

% Numerical gradient
dgdphi = numGrad(@topOpt.objective, designPar, 1e-8);

% Analytical gradient
der_g = topOpt.gradObjective(designPar);
norm(dgdphi - der_g) / norm(dgdphi)
assert(norm(dgdphi - der_g) / norm(dgdphi) < 1e-5, "Sensitivities does not match");
%%
topOpt = MaxTemperatureProblem(fem, options, 0.4);
designPar = 0.5*ones(size(fem.designPar));
topOpt.normalize(designPar);
g(1) = topOpt.objective(designPar)

% Numerical gradient
dgdphi = numGrad(@topOpt.objective, designPar, 1e-7);

% Analytical gradient
der_g = topOpt.gradObjective(designPar);
norm(dgdphi - der_g) / norm(dgdphi)
assert(norm(dgdphi - der_g) / norm(dgdphi) < 1e-5, "Sensitivities does not match");
%% GMSH mesh
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

fem = OptHeatFEM(gmsh, Elements.QUA_4, tFinal, timeSteps);
fem.addBoundaryCondition(prescribed);
fem.addBoundaryCondition(flux);
fem.addBodyCondition(body);

fem.setMaterial(material_2);

options = struct(...
    'p_kappa', 3, ...
    'p_cp', 3, ...
    'filter', true, ...
    'filterRadius', radius, ...
    'filterWeightFunction', @(dx, dy, dz) max(radius-sqrt(dx.^2+dy.^2+dz.^2), 0), ...
    'material_1', material_1, ...
    'material_2', material_2 ...
);

%%
topOpt = HeatComplianceProblem(fem, options, 0.4);
designPar = 0.5*ones(size(fem.designPar));
g(1) = topOpt.objective(designPar)

% Numerical gradient
dgdphi = numGrad(@topOpt.objective, designPar, 1e-8);

% Analytical gradient
der_g = topOpt.gradObjective(designPar);
norm(dgdphi - der_g) / norm(dgdphi)
assert(norm(dgdphi - der_g) / norm(dgdphi) < 1e-5, "Sensitivities does not match");
%%
topOpt = MaxTemperatureProblem(fem, options, 0.4);
designPar = 0.5*ones(size(fem.designPar));
topOpt.normalize(designPar);
g(1) = topOpt.objective(designPar)

% Numerical gradient
dgdphi = numGrad(@topOpt.objective, designPar, 1e-7);

% Analytical gradient
der_g = topOpt.gradObjective(designPar);
norm(dgdphi - der_g) / norm(dgdphi)
assert(norm(dgdphi - der_g) / norm(dgdphi) < 1e-5, "Sensitivities does not match");