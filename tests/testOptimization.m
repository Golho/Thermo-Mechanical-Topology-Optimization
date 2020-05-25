clear;
close all;
jobManager = JobManager();

opt.maxtime = 5*60;
opt.verbose = 1;
opt.ftol_rel = 1e-6;
%opt.xtol_abs = 1e-7*ones(size(fem.mainDensities));
opt.algorithm = NLOPT_LD_MMA;

timeSteps = 20;
tFinal = 600;
volumeFraction = 0.4;
radius = 0.005;

material_1 = Material(1, 1e6, 0.1*eye(3));
material_2 = Material(1, 1e6, 10*eye(3));

options = struct(...
    'heavisideFilter', false, ...
    'designFilter', true, ...
    'filterRadius', radius, ...
    'filterWeightFunction', @(dx, dy, dz) max(radius-sqrt(dx.^2+dy.^2+dz.^2), 0), ...
    'material_1', material_1, ...
    'material_2', material_2 ...
);
%%
gmsh = gmshParser('meshes/square_01x01_12.msh');

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

fem = OptHeatFEM(gmsh, Elements.QUA_4, tFinal, timeSteps);
fem.addBoundaryCondition(prescribed);
fem.addBoundaryCondition(flux);
fem.addBodyCondition(body);

fem.setMaterial(material_2);

[kappaF, kappaFDer, cp, cpDer] = HeatSIMP(material_1, material_2, 3, 3);
fem.addInterpFuncs(kappaF, kappaFDer, cp, cpDer);

initial = volumeFraction*ones(size(fem.designPar));

% tempFig = figure(1);
% colorbar
% designFig = figure(2);
% colorbar
% [tempPlot, designPlot] = plotIntermediate(fem, initial, tempFig, designFig);
% 
% intermediateFunc = @(femModel, designPar) plotIntermediate(femModel, designPar, tempFig, designFig, tempPlot, designPlot);
%%
topOpt_1 = HeatComplianceProblem(fem, options, volumeFraction);

job = Job(topOpt_1, initial, opt);
jobManager.add(job);
%%
gmsh = gmshParser('meshes/cube_005x005x005_1.msh');
fem = OptHeatFEM(gmsh, Elements.HEX_8, tFinal, timeSteps);
fem.addBoundaryCondition(prescribed);
fem.addBoundaryCondition(flux);
fem.addBodyCondition(body);

fem.setMaterial(material_2);
[kappaF, kappaFDer, cp, cpDer] = HeatSIMP(material_1, material_2, 3, 3);
fem.addInterpFuncs(kappaF, kappaFDer, cp, cpDer);

initial = volumeFraction*ones(size(fem.designPar));

topOpt_1 = MaxTemperatureProblem(fem, options, volumeFraction);
topOpt_1.normalize(initial);
job = Job(topOpt_1, initial, opt);
jobManager.add(job);
%%
mesh = StructuredMesh([15, 0.1], [15, 0.1]);
globalCoord = mesh.coordinates();

centerNodes = find( (globalCoord(1, :) >= 0.049 & globalCoord(1, :) <= 0.051) & ...
                    (globalCoord(2, :) >= 0.049 & globalCoord(2, :) <= 0.051));
                
cornerNodes = find( (globalCoord(1, :) == 0 | globalCoord(1, :) == 0.1) & ...
                    (globalCoord(2, :) == 0 | globalCoord(2, :) == 0.1));

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
    'value', 10/length(centerNodes), ...
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
[kappaF, kappaFDer, cp, cpDer] = HeatSIMP(material_1, material_2, 3, 3);
fem.addInterpFuncs(kappaF, kappaFDer, cp, cpDer);

topOpt_2 = HeatComplianceProblem(fem, options, volumeFraction);

initial = volumeFraction*ones(size(fem.designPar));

job = Job(topOpt_2, initial, opt);
jobManager.add(job);
%%
mesh = StructuredMesh([7, 0.05], [7, 0.05], [7, 0.05]);
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
    'value', 10/length(centerNodes), ...
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
[kappaF, kappaFDer, cp, cpDer] = HeatSIMP(material_1, material_2, 3, 3);
fem.addInterpFuncs(kappaF, kappaFDer, cp, cpDer);

topOpt_2 = MaxTemperatureProblem(copy(fem), options, volumeFraction);
initial = volumeFraction*ones(size(fem.designPar));
topOpt_2.normalize(initial);
job = Job(topOpt_2, initial, opt);
jobManager.add(job);
%%
topOpt_2 = MaxTemperatureProblem(copy(fem), options, volumeFraction);
topOpt_2.normalize(initial);
job = Job(topOpt_2, initial, opt);
jobManager.add(job);
%%
jobManager.runAll();
%%
jobManager.plotAll();