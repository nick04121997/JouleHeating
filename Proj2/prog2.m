%%%%
% Script for running the MATLAB iterative method for thickness prediction
% of a square geometry
%%%%

%% Clearing MATLAB workspace and command window, and closing all figures
clear;
close all;
clc;

%% System parameters
% tol = error tolerance for terminating iteration
% sigma = conductivity of material
% V = voltage applied to the busbars 
tol = 0.5;
global sigma V;
sigma = 1E6;
V = 115;

%% Importing data and preforming interpolation
% delta_data = initial thickness data
% heat_data = desired heating data
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

% Generates a square
geom_descrip = [3;4;0;1;1;0;0;0;1;1;];
dl = decsg(geom_descrip);
geometryFromEdges(model,dl);
pdegplot(model,'EdgeLabels','on');

% Apply boundary conditions
applyBoundaryCondition(model,'dirichlet','Edge',1,'r',-V);
applyBoundaryCondition(model,'dirichlet','Edge',3,'r',V);
applyBoundaryCondition(model,'neumann','Edge',[2,4],'g',0);

% Apply PDE coefficient
conductance_handle = @(location,state) delta_c(location.x,location.y);
specifyCoefficients(model,'m',0,'d',0,'c',conductance_handle,'a',0,'f',0,'face',1);

% Create mesh
generateMesh(model);

% Solve the PDE
result = solvepde(model);
voltage = result.NodalSolution;

% Plot solution
pdeplot(model,'XYData',voltage,'ZData',voltage);
title('voltage');

% Calculate joule heating
[e_x, e_y] = evaluateGradient(result,x,y);
q_cal = delta.*(e_x.^2+e_y.^2).^2*sigma;
err = abs(q_cal-q_des);
err = mean(err);

%% Begin iteration
i = 0;

while (err > tol)   
    fprintf('Iteration number %g \n', i)
    fprintf('The error is %g \n \n', err)
    
    % Resisty updating
    switch (mod(i,4))
        case 0
            delta = (q_cal./q_des).*delta;
        case 1
            delta = (q_des./q_cal).*delta;
        case 2
            delta = ((q_cal./q_des).^2).*delta;
        otherwise
            delta = ((q_des./q_cal).^2).*delta;
    end
    
    % Continuous fit for thickness
    delta_c = fit([x,y],delta,'linearinterp');
    
    % Apply PDE coefficients
    conductance_handle = @(location,state) delta_c(location.x,location.y);
    specifyCoefficients(model,'m',0,'d',0,'c',conductance_handle,'a',0,'f',0,'face',1);
    
    % Create mesh 
    generateMesh(model);
    
    % Solve the PDE
    result = solvepde(model);
    
    % Plot solution
    voltage = result.NodalSolution;
    pdeplot(model,'XYData',voltage,'ZData',voltage);
    title('voltage');
    
    % Caculate joule heating 
    [e_x, e_y] = evaluateGradient(result,x,y);
    q_cal = delta.*(e_x.^2+e_y.^2)*sigma;
    
    % Calcuate error
    err = abs(q_cal-q_des);
    err = mean(err);
    
    % Update iteration count
    i = i + 1;
end






