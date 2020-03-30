%% GMSH meshes
gmsh = gmshParser('meshes/long_quad.msh');
timeSteps = 100;
tFinal = 10;

k = 1;
cp = 1;
rho = 1;
q_0 = 1;
m = Material(rho, cp, k*eye(3));

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
    'value', q_0, ...
    'timeSteps', 1:timeSteps ...
);

% Create body conditions
body = struct(...
    'type', 'main', ...
    'physicalName', 'solid' ...
);

fem = HeatFEM(gmsh, Elements.QUA_4, tFinal, timeSteps);

fem.addBoundaryCondition(prescribed);
fem.addBoundaryCondition(flux);
fem.addBodyCondition(body);

fem.setMaterial(m);

fem.assemble()
fem.solve();
%% Start the plotting
% Analytical solution for a semi-infinite body
T_0 = 0;
a = k/(cp*rho);
analytical = @(x, t) 2*q_0/k*sqrt(a*t/pi)*exp(-x.^2/(4*a*t)) - q_0*x/k ...
    .*(1 - erf(x/(2*sqrt(a*t)))) + T_0;

figure(1)
colorbar
ed = fem.getElemTemp(0);
geoPlot = elfield2(fem.Ex, fem.Ey, ed);

lowerNodes = find(abs(fem.nodeCoordinates(2, :) - 0) <= 1e-6);
[~, I] = sort(fem.nodeCoordinates(1, lowerNodes));
lowerNodes = lowerNodes(I);
lowerCoords = fem.nodeCoordinates(1, lowerNodes);
lowerTemps = fem.temperatures(lowerNodes, 1);
figure(2)
lowerPlot = plot(lowerCoords, lowerTemps);
hold on
analPlot = plot(lowerCoords, analytical(lowerCoords, 0));

figure(1)
for t = 1:(timeSteps-1)
    figure(1)
    title("Frame " + t);
    ed = fem.getElemTemp(t);
    geoPlot.CData = ed;
    lowerPlot.YData = fem.temperatures(lowerNodes, t+1);
    analPlot.YData = analytical(lowerCoords, t*tFinal/(timeSteps-1));
    drawnow
end
%%
gmsh = gmshParser('meshes/square_01x01.msh');

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
    'value', 1, ...
    'timeSteps', 1:timeSteps ...
);

% Create body conditions
body = struct(...
    'type', 'main', ...
    'physicalName', 'solid' ...
);

fem = HeatFEM(gmsh, Elements.QUA_4, tFinal, timeSteps);

fem.addBoundaryCondition(prescribed);
fem.addBoundaryCondition(flux);
fem.addBodyCondition(body);

fem.setMaterial(m);

fem.assemble()
fem.solve();
%% Start the plotting
figure(3)
colorbar
ed = fem.getElemTemp(0);
geoPlot = elfield2(fem.Ex, fem.Ey, ed);

for t = 1:(timeSteps-1)
    title("Frame " + t);
    ed = fem.getElemTemp(t);
    geoPlot.CData = ed(1:4, :);
    drawnow
end
%%
gmsh = gmshParser('meshes/cube_005x005x005_1.msh');
fem = HeatFEM(gmsh, Elements.HEX_8, tFinal, timeSteps);

fem.addBoundaryCondition(prescribed);
fem.addBoundaryCondition(flux);
fem.addBodyCondition(body);

fem.setMaterial(m);

fem.assemble()
fem.solve();
%% Structured meshes
mesh = StructuredMesh([2000, 20], [2, 1]);
globalCoord = mesh.coordinates();

leftNodes = find(globalCoord(1, :) == 0);
                
rightNodes = find(globalCoord(1, :) == 20);

timeSteps = 100;
tFinal = 10;

q_0 = 1;
k = 1;
cp = 1;
rho = 1;

% Create boundary conditions
prescribed = struct(...
    'nodes', rightNodes, ...
    'type', 'Dirichlet', ...
    'value', 0, ...
    'timeSteps', 1:timeSteps ...
);

fluxCorner = struct(...
    'nodes', leftNodes, ...
    'type', 'Neumann', ...
    'value', q_0/2, ...
    'timeSteps', 1:timeSteps ...
);

% Create body conditions
body = struct(...
    'type', 'main' ...
);

fem = HeatFEMStructured(mesh, tFinal, timeSteps);

fem.addBoundaryCondition(prescribed);
fem.addBoundaryCondition(fluxCorner);
fem.addBodyCondition(body);
m = Material(rho, cp, k*eye(3));
fem.setMaterial(m);

fem.assemble()
fem.solve();
%% Start the plotting
% Analytical solution for a semi-infinite body
T_0 = 0;
a = k/(cp*rho);
analytical = @(x, t) 2*q_0/k*sqrt(a*t/pi)*exp(-x.^2/(4*a*t)) - q_0*x/k ...
    .*(1 - erf(x/(2*sqrt(a*t)))) + T_0;

figure(1)
colorbar

ed = fem.getElemTemp(0);
geoPlot = elfield2(fem.Ex(1:4, :), fem.Ey(1:4, :), ed(1:4, :));

lowerNodes = find(globalCoord(2, :) == 0 & globalCoord(3, :) == 0);
[~, I] = sort(globalCoord(1, lowerNodes));
lowerNodes = lowerNodes(I);
lowerCoords = globalCoord(1, lowerNodes);
lowerTemps = fem.temperatures(lowerNodes, 1);

figure(2)
lowerPlot = plot(lowerCoords, lowerTemps);
hold on
analPlot = plot(lowerCoords, analytical(lowerCoords, 0));
figure(1)
for t = 1:(timeSteps-1)
    figure(1)
    title("Frame " + t);
    ed = fem.getElemTemp(t);
    geoPlot.CData = ed(1:4, :);
    lowerPlot.YData = fem.temperatures(lowerNodes, t+1);
    analPlot.YData = analytical(lowerCoords, t*tFinal/(timeSteps-1));
    drawnow
end
%% 
clear;
mesh = StructuredMesh([45, 0.1], [45, 0.1]);
globalCoord = mesh.coordinates();

centerNodes = find( (globalCoord(1, :) >= 0.05 & globalCoord(1, :) <= 0.05) & ...
                    (globalCoord(2, :) >= 0.05 & globalCoord(2, :) <= 0.05));
                
cornerNodes = find( (globalCoord(1, :) == 0 | globalCoord(1, :) == 0.1) & ...
                    (globalCoord(2, :) == 0 | globalCoord(2, :) == 0.1));

timeSteps = 100;
tFinal = 0.1;

k = 1;
cp = 1;
rho = 1;

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

fem = HeatFEMStructured(mesh, tFinal, timeSteps);

fem.addBoundaryCondition(prescribed);
fem.addBoundaryCondition(fluxCorner);
fem.addBodyCondition(body);

m = Material(rho, cp, k*eye(3));
fem.setMaterial(m);

fem.assemble()
fem.solve();
%% Start the plotting
figure(3)
colorbar
ed = fem.getElemTemp(0);
geoPlot = elfield2(fem.Ex(1:4, :), fem.Ey(1:4, :), ed(1:4, :));

for t = 1:(timeSteps-1)
    title("Frame " + t);
    ed = fem.getElemTemp(t);
    geoPlot.CData = ed(1:4, :);
    drawnow
end