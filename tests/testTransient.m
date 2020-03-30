gmsh = gmshParser('meshes/square_01x01_002.msh');
timeSteps = 100;
tFinal = 0.1;

k = 1;
cp = 1;
rho = 1;

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
    'value', 10, ...
    'timeSteps', 1:timeSteps ...
);

% Create body conditions
body = struct(...
    'type', 'main', ...
    'physicalName', 'solid' ...
);

fem = HeatFEM(gmsh, tFinal, timeSteps);

fem.addBoundaryCondition(prescribed);
fem.addBoundaryCondition(flux);
fem.addBodyCondition(body);

fem.setMaterial(struct('D', k*eye(2), 'density', rho, 'heatCapacity', cp), 1);

fem.assemble()
fem.solve();
%% Start the plotting
% Analytical solution for a semi-infinite body
q_0 = 10;
T_0 = 0;
a = k/(cp*rho);
analytical = @(x, t) 2*q_0/k*sqrt(a*t/pi)*exp(-x.^2/(4*a*t)) - q_0*x/k ...
    .*(1 - erf(x/(2*sqrt(a*t)))) + T_0;

figure(1)
colorbar
[Ex, Ey, ed] = fem.getElemTemp(Elements.QUA_4, 0);

lowerNodes = find(abs(fem.nodeCoordinates(:, 2) - 0) <= 1e-6);
[~, I] = sort(fem.nodeCoordinates(lowerNodes, 1));
lowerNodes = lowerNodes(I);
lowerCoords = fem.nodeCoordinates(lowerNodes, 1);
lowerTemps = fem.temperatures(lowerNodes, 1);

geoPlot = elfield2(Ex, Ey, ed);
figure(2)
lowerPlot = plot(lowerCoords, lowerTemps);
hold on
analPlot = plot(lowerCoords, analytical(lowerCoords, 0));
figure(1)
for t = 1:(timeSteps-1)
    figure(1)
    title("Frame " + t);
    [Ex, Ey, ed] = fem.getElemTemp(Elements.QUA_4, t);
    geoPlot.CData = ed';
    lowerPlot.YData = fem.temperatures(lowerNodes, t+1);
    analPlot.YData = analytical(lowerCoords, t*tFinal/(timeSteps-1));
    drawnow
end

