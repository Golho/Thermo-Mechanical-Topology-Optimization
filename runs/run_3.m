clear; close all;

jobManager = JobManager();

deltaTemp = 100;
timeSteps = 1;
radius = 10e-6;
k = 65/10e-7;
volumeFraction = 0.2;

void = Material(0, 1, 1, 100e3, 0.3, 0);
material_1 = Material(1, 1, 1, 100e9, 0.3, 1e-5);
material_2 = Material(1, 1, 1, 100e9, 0.3, 2e-5);
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
    'value', 1e6, ...
    'components', [1, 0], ...
    'timeSteps', timeSteps ...
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

options = struct(...
    'heavisideFilter', false, ...HeavisideFilter(1, 0.5, heavisideUpdater(1, 1)), ...
    'designFilter', true, ...
    'filterRadius', radius, ...
    'filterWeightFunction', @(dx, dy, dz) max(radius-sqrt(dx.^2+dy.^2+dz.^2), 0), ...
    'materials', materials, ...
    'plot', true ...
    );

%%
massLimit = volumeFraction * sum(mechFEM.volumes*material_2.density);
initial = zeros(size(mechFEM.designPar));
initial(1, :) = 1;

opt.maxtime = 20*60;
opt.verbose = 1;
opt.ftol_rel = 1e-8;
%opt.xtol_abs = 1e-7*ones(size(fem.mainDensities));
opt.algorithm = NLOPT_LD_MMA;

for p_E = [4]
    for p_alpha = [2]
        for k = [1e7]
            i = 1;
            for beta = [1, 2, 4, 8]
                filters(1) = HeavisideFilter(beta, 0.3, heavisideUpdater(1, 1));
                filters(2) = HeavisideFilter(beta, 0.5, heavisideUpdater(1, 1));
                filters(3) = HeavisideFilter(beta, 0.7, heavisideUpdater(1, 1));
                spring = struct( ...
                    'nodes', bottomRightNode, ...
                    'type', 'Robin', ...
                    'value', k, ...
                    'components', [1, 0], ...
                    'timeSteps', 1:timeSteps ...
                    );
                
                mechFEM.addBoundaryCondition(spring);
                
                [E, EDer, alpha, alphaDer] = MechSIMP(materials, p_E, p_alpha);
                mechFEM.addInterpFuncs(E, EDer, alpha, alphaDer);
                
                %                     topOpt = MaxDisplacementProblem(mechFEM, options, massLimit);
                %                     initial = volumeFraction*ones(size(mechFEM.designPar));
                %
                %                     job = Job(topOpt, initial, opt);
                %                     jobManager.add(job);
                %
                topOpt = MaxDisplacementProblem_Robust(mechFEM, options, massLimit, filters);
                initial = 0.1*ones(size(mechFEM.designPar));
                initial(2, :) = 0.5;
                
                
                jobs(i) = Job(topOpt, initial, opt);
                if i > 1
                    opt.maxeval = 50;
                    jobs(i) = Job(topOpt, initial, opt);
                    jobs(i).linkJob(jobs(i-1));
                else
                    opt.maxeval = 150;
                    jobs(i) = Job(topOpt, initial, opt);
                end
                jobManager.add(jobs(i));
                i = i + 1;
            end
        end
    end
end
%%
jobManager.runAll();
%%
jobManager.plotAll();