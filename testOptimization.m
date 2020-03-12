gmsh = gmshParser('meshes/square_new.msh');
timeSteps = 20;
tFinal = 0.1;

% Create boundary conditions
% prescribed = struct(...
%     'physicalName', 'outlet', ...
%     'type', 'Dirichlet', ...
%     'value', 0, ...
%     'timeSteps', 1:timeSteps ...
% );

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
%fem.addBoundaryCondition(prescribed);
fem.addBoundaryCondition(flux);
fem.addBodyCondition(body);

topOpt = HeatComplianceProblem(fem, Elements.QUA_4);
initial = 0.4*ones(size(fem.mainDensities));
%%
opt.algorithm = NLOPT_LD_SLSQP;

opt.min_objective = @topOpt.nlopt_objective;
opt.lower_bounds = zeros(size(fem.mainDensities));
opt.upper_bounds = ones(size(fem.mainDensities));
opt.fc = {@topOpt.nlopt_constraint1};
opt.maxtime = 20;
opt.verbose = 1;
opt.xtol_abs = 1e-7*ones(size(fem.mainDensities));

[optDesign, fmin, returnCode] = nlopt_optimize(opt, initial);
fprintf('Return code: %d\n', returnCode);
%%
% termCriteria = @(x1, x2) norm(x1 - x2) < 1e-3;
% 
% logger = log4m.getLogger('testOpt.txt');
% 
% optimizer = OptimizerMMA(topOpt, initial, termCriteria, logger);
% optDesign = optimizer.loop();

%%
[Ex, Ey, ed] = fem.getElemTemp(timeSteps-1);

figure
title("Optimal density")
elfield2(Ex, Ey, optDesign);
colorbar

figure
title("Terminal temperature distribution");
elfield2(Ex, Ey, ed);
colorbar

lowerNodes = find(fem.nodeCoordinates(:, 2) <= 0.01);
[~, I] = sort(fem.nodeCoordinates(lowerNodes, 1));
lowerNodes = lowerNodes(I);
lowerCoords = fem.nodeCoordinates(lowerNodes, 1);
lowerTemps = fem.temperatures(lowerNodes, end);

figure
lowerPlot = plot(lowerCoords, lowerTemps);