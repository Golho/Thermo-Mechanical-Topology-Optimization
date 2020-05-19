clear; close all;
timeSteps = 1;
volumeFraction = 0.4;
radius = 0.025;

material_1 = Material(1, 1e6, 0.1*eye(3), 1e1, 0.25);
material_2 = Material(1, 1e6, 10*eye(3), 1e9, 0.25);
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
    'components', [1, 1], ...
    'timeSteps', 1:timeSteps ...
);

fluxCorner = struct(...
    'nodes', centerNodes, ...
    'type', 'Neumann', ...
    'value', 1/length(centerNodes), ...
    'components', [0, 1], ...
    'timeSteps', 1:timeSteps ...
);

% Create body conditions
body = struct(...
    'type', 'main' ...
);

fem = OptMechFEMStructured(mesh, timeSteps, "plane stress");

fem.addBoundaryCondition(prescribed);
fem.addBoundaryCondition(fluxCorner);
fem.addBodyCondition(body);

fem.setMaterial(material_2);
[E, EDer, alpha, alphaDer] = MechSIMP(material_1, material_2, 3, 3);
fem.addInterpFuncs(E, EDer, alpha, alphaDer);

options = struct(...
    "heavisideFilter", false, ...
    'designFilter', true, ...
    'filterRadius', radius, ...
    'filterWeightFunction', @(dx, dy, dz) max(radius-sqrt(dx.^2+dy.^2+dz.^2), 0), ...
    'material_1', material_1, ...
    'material_2', material_2 ...
);
%%
topOpt = MechComplianceProblem(fem, options, 0.4);
designPar = 0.5*ones(size(fem.designPar));
g(1) = topOpt.objective(designPar)

% Numerical gradient
dgdphi = numGrad(@topOpt.objective, designPar, 1e-8);

% Analytical gradient
der_g = topOpt.gradObjective(designPar);
norm(dgdphi - der_g) / norm(dgdphi)
assert(norm(dgdphi - der_g) / norm(dgdphi) < 1e-5, "Sensitivities does not match");

%%
leftCenterNode = find(globalCoord(1, :) == 0 & ...
                     (globalCoord(2, :) >= 0.046 & globalCoord(2, :) <= 0.054));

output = struct(...
    'nodes', leftCenterNode, ...
    'type', 'dummy', ...
    'name', 'josse', ...
    'value', -1, ...
    'components', [1, 0], ...
    'timeSteps', 1:timeSteps ...
);

spring = struct(...
    'nodes', leftCenterNode, ...
    'type', 'Robin', ...
    'value', 1e5, ...
    'components', [1, 0], ...
    'timeSteps', 1:timeSteps ...
);

fem.addBoundaryCondition(output);
fem.addBoundaryCondition(spring);

topOpt = FlexibilityProblem(fem, options, 0.4, 0.0005);
designPar = 0.5*ones(size(fem.designPar));
g(1) = topOpt.objective(designPar)

% Numerical gradient
dgdphi = numGrad(@topOpt.objective, designPar, 1e-8);

% Analytical gradient
der_g = topOpt.gradObjective(designPar);
norm(dgdphi - der_g) / norm(dgdphi)
assert(norm(dgdphi - der_g) / norm(dgdphi) < 1e-5, "Sensitivities does not match");
%%
g(1) = topOpt.constraint2(designPar)

% Numerical gradient
dgdphi = numGrad(@topOpt.constraint2, designPar, 1e-6);

% Analytical gradient
der_g = topOpt.gradConstraint2(designPar);
norm(dgdphi - der_g) / norm(dgdphi)
assert(norm(dgdphi - der_g) / norm(dgdphi) < 1e-5, "Sensitivities does not match");
%%
topOpt = FlexibilityProblem2(fem, options, 0.4);
designPar = 0.5*ones(size(fem.designPar));
g(1) = topOpt.objective(designPar)

% Numerical gradient
dgdphi = numGrad(@topOpt.objective, designPar, 1e-8);

% Analytical gradient
der_g = topOpt.gradObjective(designPar);
norm(dgdphi - der_g) / norm(dgdphi)
assert(norm(dgdphi - der_g) / norm(dgdphi) < 1e-5, "Sensitivities does not match");