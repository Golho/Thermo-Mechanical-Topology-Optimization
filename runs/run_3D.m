clear;
close all;
jobManager = JobManager();
%%
isnear = @(x, a) abs(x-a) < 1e-3;
mesh = StructuredMesh([41, 0.05], [41, 0.05], [41, 0.05]);
tFinal = 400;
timeSteps = 20;
radius = 0.005;
volumeFraction = 0.4;

globalCoord = mesh.coordinates();

centerNodes = find( globalCoord(1, :) == 0 & ...
    globalCoord(2, :) == 0 & ...
    globalCoord(3, :) == 0);

cornerNodes = find( globalCoord(1, :) == 0.05 & ...
    globalCoord(2, :) == 0.05 & ...
    globalCoord(3, :) == 0.05);

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

material_1 = Material(1, 5e5, 0.1*eye(3));
material_2 = Material(1, 1e6, 10*eye(3));

options = struct(...
    'p_kappa', 3, ...
    'p_cp', 3, ...
    'filter', true, ...
    'filterRadius', radius, ...
    'filterWeightFunction', @(dx, dy, dz) max(radius-sqrt(dx.^2+dy.^2+dz.^2), 0), ...
    'material_1', material_1, ...
    'material_2', material_2 ...
    );

% tempFig = figure(1);
% colorbar
% designFig = figure(2);
% colorbar
% [tempPlot, designPlot] = plotIntermediate(fem, initial, tempFig, designFig);
%
% intermediateFunc = @(femModel, designPar) plotIntermediate(femModel, designPar, tempFig, designFig, tempPlot, designPlot);

opt.maxtime = 5*60;
opt.verbose = 1;
opt.ftol_rel = 1e-6;
%opt.xtol_abs = 1e-7*ones(size(fem.mainDensities));
%%
opt.algorithm = NLOPT_LD_MMA;
fem = OptHeatFEMStructured(mesh, tFinal, timeSteps);
fem.addBoundaryCondition(prescribed);
fem.addBoundaryCondition(fluxCorner);
fem.addBodyCondition(body);

fem.setMaterial(material_2);


topOpt_1 = MaxTemperatureProblem(fem, options, volumeFraction);
initial = volumeFraction*ones(size(fem.designPar));
topOpt_1.normalize(initial);

job = Job(topOpt_1, initial, opt);
jobManager.add(job);
%%
jobManager.runAll();
%%
jobManager.plotAll();
%%
saveAnswer = questdlg("Would you like to save all jobs?", "Yes", "No");
switch saveAnswer
    case "Yes"
        jobManager.saveAll();
    case "No"
        disp("Did not save the jobs");
end