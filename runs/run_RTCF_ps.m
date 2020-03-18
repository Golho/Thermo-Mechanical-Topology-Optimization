clear;
close all;
jobManager = JobManager();
%%
gmsh = gmshParser('meshes/square_01x01.msh');
timeSteps = 50;
tFinal = 1000;
volumeFraction = 0.2;

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

material_1 = struct('kappa', 0.1, 'cp', 1e6);
material_2 = struct('kappa', 10, 'cp',  1e6);

options = struct(...
    'p_kappa', 3, ...
    'p_cp', 3, ...
    'filter', true, ...
    'filterRadius', 0.005, ...
    'material_1', material_1, ...
    'material_2', material_2 ...
);

fem.assemble();
initial = volumeFraction*ones(size(fem.mainDensities));

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
opt.algorithm = NLOPT_LD_MMA;
%%
topOpt_1 = MaxTemperatureProblem(fem, Elements.QUA_4, options, volumeFraction);
topOpt_1.normalize(initial);

job = Job(topOpt_1, initial, opt);
jobManager.add(job);
%%
options.p_kappa = 2;
options.p_cp = 3;

topOpt_1 = MaxTemperatureProblem(copy(fem), Elements.QUA_4, options, volumeFraction);
topOpt_1.normalize(initial);

job = Job(topOpt_1, initial, opt);
jobManager.add(job);
%%
options.p_kappa = 2;
options.p_cp = 2;
topOpt_1 = MaxTemperatureProblem(copy(fem), Elements.QUA_4, options, volumeFraction);
topOpt_1.normalize(initial);

job = Job(topOpt_1, initial, opt);
jobManager.add(job);
%%
options.p_kappa = 3;
options.p_cp = 2;
topOpt_1 = MaxTemperatureProblem(copy(fem), Elements.QUA_4, options, volumeFraction);
topOpt_1.normalize(initial);

job = Job(topOpt_1, initial, opt);
jobManager.add(job);
%%
options.p_kappa = 2;
options.p_cp = 4;
topOpt_1 = MaxTemperatureProblem(copy(fem), Elements.QUA_4, options, volumeFraction);
topOpt_1.normalize(initial);

job = Job(topOpt_1, initial, opt);
jobManager.add(job);
%%
options.p_kappa = 4;
options.p_cp = 2;
topOpt_1 = MaxTemperatureProblem(copy(fem), Elements.QUA_4, options, volumeFraction);
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