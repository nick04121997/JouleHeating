%% Clearing MATLAB workspace and command window, and closing all figures
clear;
close all;
clc;

%% System parameters
% tol = error tolerance for terminating iteration
% sigma = conductivity of material
% V = voltage applied to the busbars 
tol = .5;
global sigma V; 
sigma = 1E6;
V = 115;
delta_data = dlmread('thickness_init.csv',',',1,0);
heat_data = dlmread('heat_desired.csv',',',1,0);
x = delta_data(:,1);
y = delta_data(:,2);
delta = delta_data(:,3);
q_des = heat_data(:,3);

% Continuous function of desired q
q_des_c = fit([x,y],q_des,'linearinterp');

% Continuous function of delta
delta_c = fit([x,y],delta,'linearinterp');


%% PDE model setup
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

% Determine co-ordinate points for mesh nodes
x_mesh = model.Mesh.Nodes(1,:)';
y_mesh = model.Mesh.Nodes(2,:)';

% Apply PDE coefficient 
conductance_handle = @(location,state) delta_c(location.x,location.y);
specifyCoefficients(model,'m',0,'d',0,'c',conductance_handle,'a',0,'f',0,'face',1);

% Solve PDE
result = solvepde(model);

% Plot solution
voltage = result.NodalSolution;
figure(2);
pdeplot(model,'XYData',voltage,'ZData',voltage);
title('voltage');

% Query delta continuous function and heating desired continuous function
% at mesh node co-ordinate points 
delta_mesh = delta_c(x_mesh,y_mesh);
q_mesh = q_des_c(x_mesh,y_mesh);

% Calculate joule heating
[e_x, e_y] = evaluateGradient(result,x_mesh,y_mesh);
q_cal = delta_mesh.*(e_x.^2+e_y.^2).^2*sigma;

% Calculate error
err = abs(q_cal-q_mesh);
err = mean(err);

%% Begin iteration
i = 0;

while (abs(err) > tol)
    fprintf('Iteration number %g \n', i)
    fprintf('The error is %g \n \n', err)
    
    % Resisty updating 
    switch (mod(i,4))
        case 0
            delta_mesh = (q_cal./q_mesh).*delta_mesh;
        case 1
            delta_mesh = (q_mesh./q_cal).*delta_mesh;
        case 2
            delta_mesh = ((q_cal./q_mesh).^2).*delta_mesh;
        otherwise
            delta_mesh = ((q_mesh./q_cal).^2).*delta_mesh;
    end
    
    % Continuous fit for thickness
    delta_c = fit([x_mesh,y_mesh],delta_mesh,'linearinterp');
    
    % Apply PDE coefficient
    conductance_handle = @(location,state) delta_c(location.x,location.y);
    specifyCoefficients(model,'m',0,'d',0,'c',conductance_handle,'a',0,'f',0,'face',1);
    
    % Create mesh
    generateMesh(model);
    
    % Solve PDE
    result = solvepde(model);
    
    % Plot solution
    voltage = result.NodalSolution;
    pdeplot(model,'XYData',voltage,'ZData',voltage);
    title('voltage');
    
    % Calcuate joule heating 
    [e_x, e_y] = evaluateGradient(result,x_mesh,y_mesh);
    q_cal = delta_mesh.*(e_x.^2+e_y.^2)*sigma;
    
    % Calculate error
    err = abs(q_cal-q_mesh);
    err = mean(err);
    
    % Update iteration count
    i = i + 1;
end






