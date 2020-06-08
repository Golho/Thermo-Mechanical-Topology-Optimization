clear; close all;

jobManager = JobManager();

opt.maxtime = 40*60;
opt.verbose = 1;
opt.ftol_rel = 1e-8;
%opt.xtol_abs = 1e-7*ones(size(fem.mainDensities));
opt.algorithm = NLOPT_LD_MMA;

deltaTemp = 100;
timeSteps = 1;
radius = 0.01;
k = 65e6;
volumeFraction = 0.25;

void = Material(1, 1e7, 0.01*eye(3), 1e3, 0.4, 0);
copper = Material(8900, 390, 402*eye(3), 130e9, 0.355, 23e-5);
aluminium = Material(2700, 240, 88*eye(3), 70e9, 0.335, 16e-5);
materials = [void, aluminium, copper];
%%
width = 0.4;
height = 0.4;
mesh = StructuredMesh([21, width], [21, height]);
globalCoord = mesh.coordinates();


leftNodes = find(globalCoord(1, :) == 0);
centerRightNodes = find(globalCoord(1, :) == width & ...
                        globalCoord(2, :) >= height*5/11 & ...
                        globalCoord(2, :) <= height*6/11);

% Create boundary conditions
fixed = struct(...
    'nodes', leftNodes, ...
    'type', 'Dirichlet', ...
    'value', 0, ...
    'components', [1, 1], ...
    'timeSteps', 1:timeSteps ...
);

output = struct( ...
    'nodes', centerRightNodes, ...
    'type', 'dummy', ...
    'name', 'josse', ...
    'value', 1e3, ...
    'components', [0, 1], ...
    'timeSteps', timeSteps ...
);


% Create body conditions
body = struct(...
    'type', 'main' ...
);

mechFEM = OptMechFEMStructured(numel(materials), mesh, timeSteps, "plane stress");

mechFEM.addBoundaryCondition(fixed);
mechFEM.addBoundaryCondition(output);
mechFEM.addBodyCondition(body);

mechFEM.setTemperatures(deltaTemp*ones(size(mechFEM.temperatureChanges)));

mechFEM.setMaterial(aluminium);

options = struct(...
    'heavisideFilter', false, ...
    'designFilter', true, ...
    'filterRadius', radius, ...
    'filterWeightFunction', @(dx, dy, dz) max(radius-sqrt(dx.^2+dy.^2+dz.^2), 0), ...
    'materials', materials, ...
    'plot', false ...
);

%%
massLimit = volumeFraction * sum(mechFEM.volumes*copper.density);
initial = zeros(size(mechFEM.designPar));
initial(1, :) = volumeFraction;

initial(2, 1:(size(initial, 2)/2)) = 1;
initial(2, (size(initial, 2)/2):end) = 0;


p_kappa = 3;
p_cp = 1;
for p_E = [3]
    for p_alpha = [3]
        for k = 65*[1e4]
            for penalty = [0]
                spring = struct( ...
                    'nodes', centerRightNodes, ...
                    'type', 'Robin', ...
                    'value', k, ...
                    'components', [0, 1], ...
                    'timeSteps', 1:timeSteps ...
                );

                mechFEM_i = copy(mechFEM);

                mechFEM_i.addBoundaryCondition(spring);

                [E, EDer, alpha, alphaDer] = MechSIMP(materials, p_E, p_alpha);
                mechFEM_i.addInterpFuncs(E, EDer, alpha, alphaDer);

                topOpt = ThermallyActuatedProblemUniform(mechFEM_i, options, massLimit, penalty);

                job = Job(topOpt, initial, opt);
                jobManager.add(job);
            end
        end
    end
end
%%
jobManager.runAll();
%%
jobManager.plotAll();
%%
saveAnswer = questdlg("Would you like to save all jobs?", "Yes?");
switch saveAnswer
    case "Yes"
        jobManager.saveAll();
    case "No"
        disp("Did not save the jobs");
end