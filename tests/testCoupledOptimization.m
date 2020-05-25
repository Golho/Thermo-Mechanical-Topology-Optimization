clear; close all;

jobManager = JobManager();

opt.maxtime = 5*60;
opt.verbose = 1;
opt.ftol_rel = 1e-6;
%opt.xtol_abs = 1e-7*ones(size(fem.mainDensities));
opt.algorithm = NLOPT_LD_MMA;

timeSteps = 20;
volumeFraction = 0.2;
radius = 1e-5;
P = 1;
k = 500e3;
tFinal = 0.2;

material_1 = Material(1, 1e6, 1e-3*eye(3), 1e3, 0.31, 0*ones(3, 1));
material_2 = Material(1, 1e6, 10*eye(3), 1e9, 0.31, 1.5e-5*ones(3, 1));

%%
width = 5e-4;
height = 2.5e-4;
mesh = StructuredMesh([17, width], [9, height]);
globalCoord = mesh.coordinates();

topAndLeftNodes = find(globalCoord(1, :) == 0 | ...
    globalCoord(2, :) == height);
topCornerExpanded = find(globalCoord(2, :) == height & ...
    globalCoord(1, :) <= width/25);
bottomNodes = find(globalCoord(2, :) == 0);
rightCornerExpanded = find(globalCoord(1, :) == width & globalCoord(2, :) < height/25);

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
    'nodes', rightCornerExpanded, ...
    'type', 'dummy', ...
    'name', 'josse', ...
    'value', 1e8, ...
    'components', [1, 0], ...
    'timeSteps', timeSteps ...
);

spring = struct( ...
    'nodes', rightCornerExpanded, ...
    'type', 'Robin', ...
    'value', k, ...
    'components', [1, 0], ...
    'timeSteps', 1:timeSteps ...
);

tempPrescribed = struct( ...
    'nodes', rightCornerExpanded, ...
    'type', 'Dirichlet', ...
    'value', 0, ...
    'timeSteps', 1:timeSteps ...
);

heatInput = struct(...
    'nodes', topCornerExpanded, ...
    'type', 'Neumann', ...
    'value', 1, ...
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
p_E = 1;
p_kappa = 1;
p_cp = 1;
for p_alpha = [1]
    mechFEM_i = copy(mechFEM);
    
    heatFEM_i = OptThermoMechStructured(mechFEM_i, mesh, tFinal, timeSteps, 1);

    heatFEM_i.addBoundaryCondition(tempPrescribed);
    heatFEM_i.addBoundaryCondition(heatInput);
    heatFEM_i.addBodyCondition(body);

    heatFEM_i.setMaterial(material_2);
    
    [E, EDer, alpha, alphaDer] = MechSIMP(material_1, material_2, p_E, p_alpha);
    mechFEM_i.addInterpFuncs(E, EDer, alpha, alphaDer);

    [kappaF, kappaFDer, cp, cpDer] = HeatSIMP(material_1, material_2, p_kappa, p_cp);
    heatFEM_i.addInterpFuncs(kappaF, kappaFDer, cp, cpDer);

    coupledFEM = heatFEM_i;
    topOpt = ThermallyActuatedProblem(coupledFEM, options, volumeFraction);
    initial = volumeFraction*ones(size(heatFEM_i.designPar));

    job = Job(topOpt, initial, opt);
    jobManager.add(job);
end
%%
jobManager.runAll();
%%
jobManager.plotAll();