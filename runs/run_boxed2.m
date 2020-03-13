clear;
close all;
jobManager = JobManager();

gmsh = gmshParser('meshes/square_01x01.msh');
timeSteps = 50;
volumeFraction = 0.4;
material_1 = struct('kappa', 0.1, 'cp', 1e6);
material_2 = struct('kappa', 10, 'cp',  1e6);

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

options = struct(...
    'p_kappa', 3, ...
    'p_cp', 3, ...
    'filter', true, ...
    'filterRadius', 0.003, ...
    'material_1', material_1, ...
    'material_2', material_2 ...
);

opt.maxtime = 5*60;
opt.verbose = 1;
opt.ftol_rel = 1e-6;
opt.algorithm = NLOPT_LD_MMA;
%%
tFinal = 1;

fem = OptHeatFEM(gmsh, tFinal, timeSteps);
fem.addBoundaryCondition(flux);
fem.addBodyCondition(body);

fem.setMaterial(struct('D', 10*eye(2), 'density', 1, 'heatCapacity', 1e6), 0.001);

topOpt = HeatComplianceProblem(fem, Elements.QUA_4, options, volumeFraction);
initial = volumeFraction*ones(size(fem.mainDensities));

job = Job(topOpt, initial, opt);
jobManager.add(job);

%%
tFinal = 10;

fem = OptHeatFEM(gmsh, tFinal, timeSteps);
fem.addBoundaryCondition(flux);
fem.addBodyCondition(body);

fem.setMaterial(struct('D', 10*eye(2), 'density', 1, 'heatCapacity', 1e6), 0.001);

topOpt = HeatComplianceProblem(fem, Elements.QUA_4, options, volumeFraction);
initial = volumeFraction*ones(size(fem.mainDensities));

job = Job(topOpt, initial, opt);
jobManager.add(job);
%%
tFinal = 100;

fem = OptHeatFEM(gmsh, tFinal, timeSteps);
fem.addBoundaryCondition(flux);
fem.addBodyCondition(body);

fem.setMaterial(struct('D', 10*eye(2), 'density', 1, 'heatCapacity', 1e6), 0.001);

topOpt = HeatComplianceProblem(fem, Elements.QUA_4, options, volumeFraction);
initial = volumeFraction*ones(size(fem.mainDensities));

job = Job(topOpt, initial, opt);
jobManager.add(job);
%%
tFinal = 1000;

fem = OptHeatFEM(gmsh, tFinal, timeSteps);
fem.addBoundaryCondition(flux);
fem.addBodyCondition(body);

fem.setMaterial(struct('D', 10*eye(2), 'density', 1, 'heatCapacity', 1e6), 0.001);

topOpt = HeatComplianceProblem(fem, Elements.QUA_4, options, volumeFraction);
initial = volumeFraction*ones(size(fem.mainDensities));

job = Job(topOpt, initial, opt);
jobManager.add(job);
%%
tFinal = 1000;

fem = OptHeatFEM(gmsh, tFinal, timeSteps);
fem.addBoundaryCondition(flux);
fem.addBodyCondition(body);

fem.setMaterial(struct('D', 10*eye(2), 'density', 1, 'heatCapacity', 1e6), 0.001);

topOpt = HeatComplianceProblem(fem, Elements.QUA_4, options, volumeFraction);
initial = volumeFraction*ones(size(fem.mainDensities));

job = Job(topOpt, initial, opt);
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