clear; close all;

jobManager = JobManager();

opt.maxtime = 20*60;
opt.verbose = 1;
opt.ftol_rel = 1e-6;
%opt.xtol_abs = 1e-7*ones(size(fem.mainDensities));
opt.algorithm = NLOPT_LD_MMA;

deltaTemp = 100;
tFinal = 1;
timeSteps = 20;
radius = 10e-6;
k = 65/10e-6;
volumeFraction = 0.25;

void = Material(1, 1e6, 0.01*eye(3), 100e3, 0.3, 0);
material_1 = Material(2000, 1e3, 10*eye(3), 100e9, 0.3, 4e-5);
material_2 = Material(1000, 1e3, 10*eye(3), 100e9, 0.3, 2e-5);

materials = [void, material_1, material_2];

%%
width = 400e-6;
height = 200e-6;
mesh = StructuredMesh([81, width], [41, height]);
globalCoord = mesh.coordinates();


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
    'value', 1e4, ...
    'components', [1, 0], ...
    'timeSteps', timeSteps ...
);

flux = struct(...
    'nodes', topCornerExpanded, ...
    'type', 'Neumann', ...
    'value', 1, ...
    'timeSteps', 1:timeSteps ...
);

% Create body conditions
body = struct(...
    'type', 'main' ...
);

mechFEM = OptMechFEMStructured(numel(materials), mesh, timeSteps, "plane stress");

mechFEM.addBoundaryCondition(fixed);
mechFEM.addBoundaryCondition(symmetry);
mechFEM.addBoundaryCondition(output);
mechFEM.addBodyCondition(body);

mechFEM.setTemperatures(deltaTemp*ones(size(mechFEM.temperatureChanges)));

mechFEM.setMaterial(material_2);
massLimit = volumeFraction * sum(mechFEM.volumes*material_2.density);

options = struct(...
    'heavisideFilter', false, ...
    'designFilter', true, ...
    'filterRadius', radius, ...
    'filterWeightFunction', @(dx, dy, dz) max(radius-sqrt(dx.^2+dy.^2+dz.^2), 0), ...
    'materials', materials, ...
    'plot', true ...
);

%%
p_kappa = 3;
p_cp = 1;
for p_E = [3]
    for p_alpha = [3]
        mechFEM_i = copy(mechFEM);

        heatFEM_i = OptThermoMechStructured(mechFEM_i, numel(materials), mesh, tFinal, timeSteps, 1);

        heatFEM_i.addBoundaryCondition(flux);
        heatFEM_i.addBodyCondition(body);

        heatFEM_i.setMaterial(material_2);

        [E, EDer, alpha, alphaDer] = MechSIMP(materials, p_E, p_alpha);
        mechFEM_i.addInterpFuncs(E, EDer, alpha, alphaDer);

        [kappaF, kappaFDer, cp, cpDer] = HeatSIMP(materials, p_kappa, p_cp);
        heatFEM_i.addInterpFuncs(kappaF, kappaFDer, cp, cpDer);

        coupledFEM = heatFEM_i;
        topOpt = ThermallyActuatedProblem(coupledFEM, options, massLimit);
        initial = volumeFraction * ones(size(heatFEM_i.designPar));

        job = Job(topOpt, initial, opt);
        jobManager.add(job);
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