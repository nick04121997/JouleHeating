%% create PDE model
model = createpde();
global sigma; 
sigma = 1E6;

%% PDE model setup
model = createpde();

% Generates a square
geom_descrip = [3;4;0;1;1;0;0;0;1;1;];
dl = decsg(geom_descrip);
geometryFromEdges(model,dl);
pdegplot(model,'EdgeLabels','on');

% Apply boundary conditions
applyBoundaryCondition(model,'dirichlet','Edge',1,'r',-V);
applyBoundaryCondition(model,'dirichlet','Edge',3,'r',V);
applyBoundaryCondition(model,'neumann','Edge',[2,4],'g',0);

delta_c = @(y) 1E-8*(1+y);

% Apply PDE coefficient
conductance_handle = @(location,state) delta_c(location.y);
specifyCoefficients(model,'m',0,'d',0,'c',conductance_handle,'a',0,'f',0,'face',1);

% Create mesh
generateMesh(model);

% Solve the PDE
result = solvepde(model);
voltage = result.NodalSolution;

x_mesh = result.Mesh.Nodes(1,:)';
y_mesh = result.Mesh.Nodes(2,:)';

voltage = result.NodalSolution;
figure(2);
pdeplot(model,'XYData',voltage,'ZData',voltage);
title('voltage');

figure(3);
pdemesh(model)

delta_mesh = delta_c(y_mesh);

[e_x, e_y] = evaluateGradient(result,x_mesh,y_mesh);

q_cal = delta_mesh.*(e_x.^2+e_y.^2).*sigma;

q_cal_fit = fit([x_mesh,y_mesh],q_cal,'linearinterp');







