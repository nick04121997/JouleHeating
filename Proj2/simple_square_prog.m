%%%%
% simple_square_prog.m
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

%% Create solution movie
% Help visualize thickness updating and convergence
v = VideoWriter('solution.avi');
v.FrameRate = 2;
open(v);

%% Importing data and preforming interpolation
% delta_data: initial thickness data
% heat_data: desired heating data
delta_data = dlmread('thickness_init.csv',',',0,0);
heat_data = dlmread('qj_test.csv',',',0,0);
x = delta_data(:,1);
y = delta_data(:,2);
delta = delta_data(:,3);
q_des = heat_data(:,3);

% Continuous function of desired q
q_des_c = fit([x,y],q_des,'linearinterp');
% Continuous function of delta
delta_c = fit([x,y],delta,'linearinterp');

%% PDE model setup
% Create pde model
model = createpde();

%% PDE model gemeotry
% In this case, a square
geom_descrip = [3;4;0;1;1;0;0;0;1;1;];
dl = decsg(geom_descrip);
geometryFromEdges(model,dl);
figure(1);
pdegplot(model,'EdgeLabels','on');
title('Geometry')
xlabel('x (m)')
ylabel('y (m)')
set(gca, 'FontSize', 14)

%% Apply boundary conditions
applyBoundaryCondition(model,'dirichlet','Edge',1,'r',-V);
applyBoundaryCondition(model,'dirichlet','Edge',3,'r',V);
applyBoundaryCondition(model,'neumann','Edge',[2,4],'g',0);

%% Apply PDE coefficient
conductance_handle = @(location,state) delta_c(location.x,location.y);
specifyCoefficients(model,'m',0,'d',0,'c',conductance_handle,'a',0,'f',0,'face',1);

%% Create mesh
generateMesh(model);

%% Solve the PDE
figure(2);
result = solvepde(model);
voltage = result.NodalSolution;
% Plot solution
pdeplot(model,'XYData',voltage,'ZData',voltage);
title('Initial Voltage');
colorbar
xlabel('x (m)')
ylabel('y (m)')
drawnow

%% Calculate joule heating
[e_x, e_y] = evaluateGradient(result,x,y);
q_cal = delta.*(e_x.^2+e_y.^2)*sigma;
err = abs(q_cal-q_des);
err = mean(err);

%% Begin iteration
i = 0;

while (err > tol)
    if mod(i,4) == 0
        % Print iteration number and error for every major iteration
        fprintf('Iteration number %g \n', i/4+1)
        fprintf('The error is %g \n \n', err)
        
        % Plot for solution movie
        fig = figure(5);
        set(fig, 'Visible', 'off')
        plot(0:0.1:1, delta_c(0.4*ones(1,length(0:0.1:1)), 0:0.1:1), 'LineWidth', 2);
        set(gca, 'FontSize', 14);
        ytickformat('%.2f')
        xlabel('y (m)');
        ylabel('Thickness (m)');
        title('Thickness Solution');
        axis([0 1 0.8E-8 2.2E-8])
        hold on;
        frame = getframe(gcf);
        writeVideo(v,frame);
    end
    
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
    figure(2);
    voltage = result.NodalSolution;
    pdeplot(model,'XYData',voltage,'ZData',voltage);
    
    % Caculate joule heating 
    [e_x, e_y] = evaluateGradient(result,x,y);
    q_cal = delta.*(e_x.^2+e_y.^2)*sigma;
    
    % Calcuate error
    err = abs(q_cal-q_des);
    err = mean(err);
    
    % Update iteration count
    i = i + 1;
end

% Print final total iteration count, error, and standard deviation
fprintf('After %g iterations \n', i);
fprintf('The error is %g \n', mean(q_cal-q_des));
fprintf('The standard deviation is %g \n', std(q_cal-q_des));

% Editing voltage plot
figure(2)
title('Converged Voltage');
colorbar
xlabel('x (m)')
ylabel('y (m)')
drawnow

%% Plot desired heating
figure(3)
plot(q_des_c, [x,y], q_des)
xlabel('x (m)');
ylabel('y (m)');
title('Heat Profile for Linear Thickness');
colorbar;
set(gca, 'FontSize', 14)

%% Plot converged thickness
figure(4)
plot(delta_c,[x,y],delta);
xlabel('x (m)');
ylabel('y (m)');
title('Predicted Thickness Profile');
colorbar;
set(gca, 'FontSize', 14)

%% Plot for solution movie
fig = figure(5);
set(fig, 'Visible', 'on')
plot(0:0.1:1, delta_c(0.4*ones(1,length(0:0.1:1)), 0:0.1:1), 'LineWidth', 2);
set(gca, 'FontSize', 14);
ytickformat('%.2f')
xlabel('y (m)');
ylabel('Thickness (m)');
title('Thickness Solution');
axis([0 1 0.8E-8 2.2E-8])
hold on;
legend('1', '2', '3', '4', '5', '6', 'Location', 'southeast');
frame = getframe(gcf);
writeVideo(v,frame);
close(v);
