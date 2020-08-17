clear; close all;

global mechFEM
global materials
global mesh
global tempPrescribed
global body
global opt
global jobManager

jobManager = JobManager();
jobs = Job.empty();

opt.maxtime = 20*60;
opt.verbose = 1;
opt.ftol_rel = 1e-6;
%opt.xtol_abs = 1e-7*ones(size(fem.mainDensities));
opt.algorithm = NLOPT_LD_MMA;

tFinal = 0.2;
timeSteps = 50;
radius = 10e-6;
k = 1e7;
volumeFraction = 0.2;

void = Material(1e-3, 1e9, 1e-3, 100e3, 0.3, 0);
material_1 = Material(1, 1e6, 10, 100e9, 0.3, 1e-5);
material_2 = Material(1, 1e6, 10, 100e9, 0.3, 2e-5);
materials = [void, material_1, material_2];
%% Define the geometry
width = 400e-6;
height = 200e-6;
mesh = StructuredMesh([141, width], [71, height]);
globalCoord = mesh.coordinates();

% Define the node groups
topAndLeftNodes = find(globalCoord(1, :) == 0 | ...
    globalCoord(2, :) == height);
topRightNode = find(globalCoord(1, :) == width & globalCoord(2, :) == height);
topCornerExpanded = find(globalCoord(2, :) == height & ...
    globalCoord(1, :) <= width/25);
bottomNodes = find(globalCoord(2, :) == 0);
bottomRightNode = find(globalCoord(2, :) == 0 & globalCoord(1, :) == width);
rightCornerExpanded = find((globalCoord(1, :) >= width*9/10 & globalCoord(2, :) == 0) | ...
    (globalCoord(2, :) <= height/10 & globalCoord(1, :) == width));
rightCorner = find(globalCoord(1, :) >= width*9/10 & globalCoord(2, :) == 0);

% Create boundary conditions
fixed = struct(...
    'nodes', topAndLeftNodes, ...
    'type', 'Dirichlet', ...
    'value', 0, ...
    'components', [1, 1], ...
    'timeSteps', 1:timeSteps ...
    );

symmetry = struct(...
    'nodes', bottomNodes, ...
    'type', 'Dirichlet', ...
    'value', 0, ...
    'components', [0, 1], ...
    'timeSteps', 1:timeSteps ...
    );

output = struct( ...
    'nodes', bottomRightNode, ...
    'type', 'dummy', ...
    'name', 'josse', ...
    'value', 1e7, ...
    'components', [1, 0], ...
    'timeSteps', timeSteps ...
    );

spring = struct( ...
    'nodes', bottomRightNode, ...
    'type', 'Robin', ...
    'value', k, ...
    'components', [1, 0], ...
    'timeSteps', timeSteps ...
    );


tempPrescribed = struct(...
    'nodes', topCornerExpanded, ...
    'type', 'Dirichlet', ...
    'value', 100, ...
    'timeSteps', 1:timeSteps ...
    );

% Create body conditions
body = struct(...
    'type', 'main' ...
    );

% Create the FEM model and add the boundary/body conditions
mechFEM = OptMechFEMStructured(numel(materials), mesh, timeSteps, "plane stress");

mechFEM.addBoundaryCondition(fixed);
mechFEM.addBoundaryCondition(symmetry);
mechFEM.addBoundaryCondition(output);
mechFEM.addBoundaryCondition(spring);
mechFEM.addBodyCondition(body);



options = struct(...
    'heavisideFilter', false, ...
    'designFilter', true, ...
    'filterRadius', radius, ...
    'filterWeightFunction', @(dx, dy, dz) max(radius-sqrt(dx.^2+dy.^2+dz.^2), 0), ...
    'materials', materials, ...
    'plot', true ...
    );

%%
opt.maxeval = 70;
p_cp = 3;
p_kappa = 2;
p_E = 3;
p_alpha = 4;

for tFinal = [1e-3, 1e-2, 1e-1]
    addJob(tFinal, p_cp, p_kappa, p_E, p_alpha);
end

tFinal = 1e-2;
p_cp = 3;
p_kappa = 3;
p_E = 3;
p_alpha = 3;

addJob(tFinal, p_cp, p_kappa, p_E, p_alpha);

for p_cp = [2, 4]
    addJob(tFinal, p_cp, p_kappa, p_E, p_alpha);
end

p_cp = 3;

for p_kappa = [2, 4]
    addJob(tFinal, p_cp, p_kappa, p_E, p_alpha);
end

p_kappa = 3;

for p_E = [2, 4]
    addJob(tFinal, p_cp, p_kappa, p_E, p_alpha);
end

p_E = 3;

for p_alpha = [2, 4]
    addJob(tFinal, p_cp, p_kappa, p_E, p_alpha);
end

p_alpha = 3;


%%

%%
jobManager.runAndSaveAll();
%%
jobManager.plotAll();
%%
% saveAnswer = questdlg("Would you like to save all jobs?", "Yes", "No");
% switch saveAnswer
%     case "Yes"
%         jobManager.saveAll();
%     case "No"
%         disp("Did not save the jobs");
% end

function addJob(tFinal, p_cp, p_kappa, p_E, p_alpha)
    global mechFEM
    global materials
    global mesh
    global tempPrescribed
    global body
    global opt
    global jobManager
    
    tFinal = 0.2;
    timeSteps = 50;
    radius = 10e-6;
    k = 1e7;
    volumeFraction = 0.2;
    
    massLimit = volumeFraction * sum(mechFEM.volumes);%*material_1.density);
    
    heatFEM_i = OptThermoMechStructured(mechFEM, numel(materials), mesh, tFinal, timeSteps, 1);

    heatFEM_i.addBoundaryCondition(tempPrescribed);
    heatFEM_i.addBodyCondition(body);

    [E, EDer, alpha, alphaDer] = MechSIMP(materials, p_E, p_alpha);
    heatFEM_i.mechFEM.addInterpFuncs(E, EDer, alpha, alphaDer);

    [kappaF, kappaFDer, cp, cpDer] = HeatSIMP(materials, p_kappa, p_cp);
    heatFEM_i.addInterpFuncs(kappaF, kappaFDer, cp, cpDer);

    coupledFEM = heatFEM_i;
    i = 1;
    for beta = [1, 2, 3, 4]
        filters(1) = HeavisideFilter(beta, 0.3, heavisideUpdater(1, 1));
        filters(2) = HeavisideFilter(beta, 0.5, heavisideUpdater(1, 1));
        filters(3) = HeavisideFilter(beta, 0.7, heavisideUpdater(1, 1));
        options = struct(...
            'heavisideFilter', false, ...
            'designFilter', true, ...
            'filterRadius', radius, ...
            'filterWeightFunction', @(dx, dy, dz) max(radius-sqrt(dx.^2+dy.^2+dz.^2), 0), ...
            'materials', materials, ...
            'plot', true ...
            );
        topOpt = MaxDisplacementProblem_Coupled_Robust2(coupledFEM, options, massLimit, 0.5, filters);

        initial = volumeFraction * ones(size(heatFEM_i.designPar));
        %initial(1, :) = 0.1;
        initial(2, :) = 0.5;

        if i > 1
            opt.maxeval = 50;
            jobs(i) = Job(topOpt, initial, opt);
            jobs(i).linkJob(jobs(i-1));
        else
            opt.maxeval = 60;
            jobs(i) = Job(topOpt, initial, opt);
        end
        jobManager.add(jobs(i));
        i = i + 1;
    end
end