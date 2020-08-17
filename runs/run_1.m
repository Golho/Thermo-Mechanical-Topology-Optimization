clear;
close all;
jobManager = JobManager();
% 2D max temp optimization, variation of time and penalization parameter
%%
width = 0.05;
height = 0.05;
mesh = StructuredMesh([101, width], [101, height]);
timeSteps = 50;
radius = 0.005;
volumeFraction = 0.2;

globalCoord = mesh.coordinates();

centerNode = find( globalCoord(1, :) == 0 & globalCoord(2, :) == 0);
                
cornerNode = find( globalCoord(1, :) == width & globalCoord(2, :) == height);

% Create boundary conditions
prescribed = struct(...
    'nodes', cornerNode, ...
    'type', 'Dirichlet', ...
    'value', 0, ...
    'timeSteps', 1:timeSteps ...
);

fluxCorner = struct(...
    'nodes', centerNode, ...
    'type', 'Neumann', ...
    'value', 1/length(centerNode), ...
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

opt.maxtime = 20*60;
opt.verbose = 1;
opt.ftol_rel = 1e-7;
%opt.xtol_abs = 1e-7*ones(size(fem.mainDensities));
%%
opt.algorithm = NLOPT_LD_MMA;
for tFinal = [100, 500, 2000, 10000]
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
jobManager.runAndSaveAll();
%%
jobManager.plotAll();