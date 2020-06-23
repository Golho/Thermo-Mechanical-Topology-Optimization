clear;
close all;
jobManager = JobManager();
%%
isnear = @(x, a) abs(x-a) < 1e-3; 
mesh = StructuredMesh([81, 0.1], [81, 0.1]);
timeSteps = 50;
radius = 0.005;
volumeFraction = 0.2;

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

material_1 = Material(1, 5e5, 0.1);
material_2 = Material(1e3, 1e3, 10);
materials = [material_1, material_2];

options = struct(...
    'heavisideFilter', false, ...
    'designFilter', true, ...
    'filterRadius', radius, ...
    'filterWeightFunction', @(dx, dy, dz) max(radius-sqrt(dx.^2+dy.^2+dz.^2), 0), ...
    'materials', materials, ...
    'plot', false ...
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
for tFinal = [100, 200, 500, 2000, 6000, 10000]
    fem = OptHeatFEMStructured(numel(materials), mesh, tFinal, timeSteps);
    fem.addBoundaryCondition(prescribed);
    fem.addBoundaryCondition(fluxCorner);
    fem.addBodyCondition(body);

    fem.setMaterial(material_2);
    massLimit = volumeFraction * sum(fem.volumes * material_2.density);
    
    for p_kappa = [2, 3]
        for p_cp = [2, 3]
            [kappaF, kappaFDer, cp, cpDer] = HeatSIMP(materials, p_kappa, p_cp);
            fem.addInterpFuncs(kappaF, kappaFDer, cp, cpDer);

            topOpt_1 = MaxTemperatureProblem(fem, options, massLimit);
            initial = volumeFraction*ones(size(fem.designPar));
            topOpt_1.normalize(initial);

            job = Job(topOpt_1, initial, opt);
            jobManager.add(job);
        end
    end
end
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