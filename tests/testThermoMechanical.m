%% GMSH meshes
clear; close all;
mesh = StructuredMesh([31, 1], [31, 1]);
globalCoord = mesh.coordinates();

centerNodes = find( (globalCoord(1, :) >= 0.49 & globalCoord(1, :) <= 0.51) & ...
                    (globalCoord(2, :) >= 0.49 & globalCoord(2, :) <= 0.51));
                
cornerNodes = find( (globalCoord(1, :) == 0 | globalCoord(1, :) == 1) & ...
                    (globalCoord(2, :) == 0 | globalCoord(2, :) == 1));

timeSteps = 10;
tFinal = 1000;

k = 400;
cp = 385;
rho = 8900;
alpha = 16.4e-6*ones(3, 1);
m = Material(rho, cp, k*eye(3), 110e9, 0.34, alpha);

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
    'value', 10000/length(centerNodes), ...
    'timeSteps', 1:timeSteps ...
);

% Create body conditions
body = struct(...
    'type', 'main' ...
);

femTemp = HeatFEMStructured(mesh, tFinal, timeSteps);

femTemp.addBoundaryCondition(prescribed);
femTemp.addBoundaryCondition(fluxCorner);
femTemp.addBodyCondition(body);

femTemp.setMaterial(m);

femTemp.assemble()
femTemp.solve();
%%
bottomNodes = find(globalCoord(2, :) == 0);
bottomLeftNode = find(globalCoord(1, bottomNodes) == 0);
topNodes = find(globalCoord(2, :) == 1);

% Create boundary conditions
fixed = struct(...
    'nodes', bottomNodes, ...
    'type', 'Dirichlet', ...
    'value', 0, ...
    'components', [0, 1], ...
    'timeSteps', 1:timeSteps ...
);

pulled = struct(...
    'nodes', bottomLeftNode, ...
    'type', 'Dirichlet', ...
    'value', 0, ...
    'components', [1, 0], ...
    'timeSteps', 1:timeSteps ...
);

% Create body conditions
body = struct(...
    'type', 'main' ...
);

fem = MechFEMStructured(mesh, timeSteps, "plane stress");

fem.addBoundaryCondition(fixed);
fem.addBoundaryCondition(pulled);
fem.addBodyCondition(body);

fem.setMaterial(m);
fem.setTemperatures(femTemp.temperatures);

fem.assemble()
fem.solve();
%% Start the plotting
figure(1)
sgtitle("Timestep 0");
subplot(2, 2, 1)
colorbar
eTemp = femTemp.getElemTemp(0);
tempPlot = elfield2(femTemp.Ex, femTemp.Ey, eTemp);
% Analytical solution for a semi-infinite body
subplot(2, 2, 2)
colorbar
eDisp = fem.getElemDisp(0);
dispPlot = eldisplace2(fem.Ex, fem.Ey, eDisp, 1000);

% Stresses
subplot(2, 2, 3)
colorbar
eStress = fem.getElemStress(0);
vmStress = vonMises(eStress);
stressPlot = elfield2(fem.Ex, fem.Ey, vmStress);

for t = 1:timeSteps-1
    sgtitle("Timestep " + t);
    eTemp = femTemp.getElemTemp(t);
    tempPlot.CData = eTemp;
    eDisp = fem.getElemDisp(t);
    dispPlot.CData = eDisp;
    eStress = fem.getElemStress(t);
    vmStress = vonMises(eStress);
    stressPlot.CData = vmStress;
    drawnow;
end