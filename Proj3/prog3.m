%% setup
clear;
close all;
clc;

%% data import and interp
tol = .5;
global sigma; 
sigma = 1E6;
delta_data = dlmread('thickness_init.csv',',',1,0);
heat_data = dlmread('heat_desired.csv',',',1,0);
x = delta_data(:,1);
y = delta_data(:,2);
delta = delta_data(:,3);
q_des = heat_data(:,3);


%% create PDE model
model = createpde();

x_bottom = 0.0254*[0 3.622 7.411 11.107 14.724 18.279 21.772 25.214 28.611 32.037 35.463];
y_bottom = 0.0254*[0 0.815 3.101 5.532 8.081 10.715 13.431 16.211 19.046 22.421 25.795];

x_top = 0.0254*[0 2.712 5.378 7.994 10.566 13.09 15.567 18.005 20.404 22.764 25.09];
y_top = 0.0254*[20.902 22.085 23.37 24.754 26.217 27.76 29.38 31.057 32.79 34.576 36.404];

geom_descrip = [2;22;x_bottom'; flip(x_top');y_bottom';flip(y_top')];
dl = decsg(geom_descrip);
geometryFromEdges(model,dl);
figure(1);
pdegplot(model,'EdgeLabels','on');

applyBoundaryCondition(model,'dirichlet','Edge',[1:1:10],'r',-115);
applyBoundaryCondition(model,'dirichlet','Edge',[12:1:21],'r',115);
applyBoundaryCondition(model,'neumann','Edge',[11,22],'g',0);

generateMesh(model);
x_mesh = model.Mesh.Nodes(1,:)';
y_mesh = model.Mesh.Nodes(2,:)';

% continuous function of desired q
q_des_c = fit([x,y],q_des,'linearinterp');

% continuous function of delta
delta_c = fit([x,y],delta,'linearinterp');

conductance_handle = @(location,state) delta_c(location.x,location.y);
specifyCoefficients(model,'m',0,'d',0,'c',conductance_handle,'a',0,'f',0,'face',1);

result = solvepde(model);


voltage = result.NodalSolution;
figure(2);
pdeplot(model,'XYData',voltage,'ZData',voltage);
title('voltage');

delta_mesh = delta_c(x_mesh,y_mesh);
q_mesh = q_des_c(x_mesh,y_mesh);

[e_x, e_y] = evaluateGradient(result,x_mesh,y_mesh);

q_cal = delta_mesh.*(e_x.^2+e_y.^2).^2*sigma;
err = q_cal-q_mesh;
err = mean(err);

%% Iterate
i = 0;
err

while (abs(err) > tol)
    i
    err
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
    
    delta_c = fit([x_mesh,y_mesh],delta_mesh,'linearinterp');
    
    conductance_handle = @(location,state) delta_c(location.x,location.y);
    specifyCoefficients(model,'m',0,'d',0,'c',conductance_handle,'a',0,'f',0,'face',1);
    generateMesh(model);

    result = solvepde(model);
    voltage = result.NodalSolution;
    pdeplot(model,'XYData',voltage,'ZData',voltage);
%     drawnow;
    title('voltage');
    
    [e_x, e_y] = evaluateGradient(result,x_mesh,y_mesh);

    q_cal = delta_mesh.*(e_x.^2+e_y.^2)*sigma;
    err = q_cal-q_mesh;
    err = mean(err);
    i = i + 1;
end






