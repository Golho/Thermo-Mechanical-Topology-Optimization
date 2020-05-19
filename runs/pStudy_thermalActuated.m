clear; close all;

jobManager = JobManager();

opt.maxtime = 5*60;
opt.verbose = 1;
opt.ftol_rel = 1e-7;
%opt.xtol_abs = 1e-7*ones(size(fem.mainDensities));
opt.algorithm = NLOPT_LD_MMA;

timeSteps = 50;
radius = 0.005;
P = 1;
k = 5e5;
tFinal = 2000;

material_1 = Material(1, 1.767e6, 1e-2*eye(3), 3e2, 0.45, 0*ones(3, 1));
material_2 = Material(1, 1.767e6, 0.22*eye(3), 1.1e9, 0.45, 8e-5*ones(3, 1));

%%
width = 0.1;
height = 0.1;
mesh = StructuredMesh([31, width], [31, height]);
globalCoord = mesh.coordinates();


topAndLeftNodes = find(globalCoord(1, :) == 0 | ...
    globalCoord(2, :) == height);
topRightNode = find(globalCoord(1, :) == height & globalCoord(2, :) == width);
topCornerExpanded = find(globalCoord(2, :) == height & ...
    globalCoord(1, :) <= width/25);
bottomNodes = find(globalCoord(2, :) == 0);
rightCornerExpanded = find((globalCoord(1, :) >= width*9/10 & globalCoord(2, :) == 0) | ...
                           (globalCoord(2, :) <= height/10 & globalCoord(1, :) == width));
rightCorner = find(globalCoord(1, :) >= width*9/10 & globalCoord(2, :) == 0);

% Create boundary conditions
fixed = struct(...
    'nodes', rightCornerExpanded, ...
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
    'nodes', topRightNode, ...
    'type', 'dummy', ...
    'name', 'josse', ...
    'value', 1e4, ...
    'components', [1, 0], ...
    'timeSteps', timeSteps ...
);

spring = struct( ...
    'nodes', topRightNode, ...
    'type', 'Robin', ...
    'value', k, ...
    'components', [1, 0], ...
    'timeSteps', 1:timeSteps ...
);

tempPrescribed = struct( ...
    'nodes', bottomNodes, ...
    'type', 'Dirichlet', ...
    'value', 100, ...
    'timeSteps', 1:timeSteps ...
);


% Create body conditions
body = struct(...
    'type', 'main' ...
);

mechFEM = OptMechFEMStructured(mesh, timeSteps, "plane stress");

mechFEM.addBoundaryCondition(fixed);
mechFEM.addBoundaryCondition(symmetry);
mechFEM.addBoundaryCondition(output);
mechFEM.addBoundaryCondition(spring);
mechFEM.addBodyCondition(body);

mechFEM.setMaterial(material_2);

options = struct(...
    'heavisideFilter', false, ...
    'designFilter', true, ...
    'filterRadius', radius, ...
    'filterWeightFunction', @(dx, dy, dz) max(radius-sqrt(dx.^2+dy.^2+dz.^2), 0), ...
    'material_1', material_1, ...
    'material_2', material_2 ...
);

%%
p_kappa = 3;
p_cp = 1;
for p_E = [3]
    for p_alpha = [0.25, 0.5, 1, 2, 3]
        mechFEM_i = copy(mechFEM);

        heatFEM_i = OptThermoMechStructured(mechFEM_i, mesh, tFinal, timeSteps, 1);

        heatFEM_i.addBoundaryCondition(tempPrescribed);
        heatFEM_i.addBodyCondition(body);

        heatFEM_i.setMaterial(material_2);

        [E, EDer, alpha, alphaDer] = MechSIMP(material_1, material_2, p_E, p_alpha);
        mechFEM_i.addInterpFuncs(E, EDer, alpha, alphaDer);

        [kappaF, kappaFDer, cp, cpDer] = HeatSIMP(material_1, material_2, p_kappa, p_cp);
        heatFEM_i.addInterpFuncs(kappaF, kappaFDer, cp, cpDer);

        coupledFEM = heatFEM_i;
        topOpt = ThermallyActuatedProblem2(coupledFEM, options);
        initial = 0.5*ones(size(heatFEM_i.designPar));

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