%% create PDE model
global sigma; 
sigma = 1E6;
V = 115;

model = createpde();

% Bottom busbar co-ordinates
x_bottom = 0.0254*[0 3.622 7.411 11.107 14.724 18.279 21.772 25.214 28.611 32.037 35.463];
y_bottom = 0.0254*[0 0.815 3.101 5.532 8.081 10.715 13.431 16.211 19.046 22.421 25.795];

% Top busbar co-ordinates
x_top = 0.0254*[0 2.712 5.378 7.994 10.566 13.09 15.567 18.005 20.404 22.764 25.09];
y_top = 0.0254*[20.902 22.085 23.37 24.754 26.217 27.76 29.38 31.057 32.79 34.576 36.404];

% Defining geometry
geom_descrip = [2;22;x_bottom'; flip(x_top');y_bottom';flip(y_top')];
dl = decsg(geom_descrip);
geometryFromEdges(model,dl);

% Plot geometry
figure(1);
pdegplot(model,'EdgeLabels','on');

% Apply boundary conditions
applyBoundaryCondition(model,'dirichlet','Edge',[1:1:10],'r',-V);
applyBoundaryCondition(model,'dirichlet','Edge',[12:1:21],'r',V);
applyBoundaryCondition(model,'neumann','Edge',[11,22],'g',0);

% Create mesh
generateMesh(model);

delta_c = @(x,y) 1E-8*(1+x+y);

% Apply PDE coefficient
conductance_handle = @(location,state) delta_c(location.x,location.y);
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

delta_mesh = delta_c(x_mesh,y_mesh);

[e_x, e_y] = evaluateGradient(result,x_mesh,y_mesh);

q_cal = delta_mesh.*(e_x.^2+e_y.^2).*sigma;

q_cal_fit = fit([x_mesh,y_mesh],q_cal,'linearinterp');

figure(4)
plot(q_cal_fit, [x_mesh,y_mesh], q_cal);

delta_fit = fit([x_mesh,y_mesh], delta_mesh,'linearinterp');

figure(5)
plot(delta_fit, [x_mesh,y_mesh], delta_mesh);







