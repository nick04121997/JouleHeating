%% setup
clear;
close all;
clc;

%% data import and interp
tol = 10;
global sigma; 
sigma = 1E6;
delta_data = dlmread('thickness_init.csv',',',1,0);
heat_data = dlmread('heat_desired.csv',',',1,0);
x = delta_data(:,1);
y = delta_data(:,2);
delta = delta_data(:,3);
q_des = heat_data(:,3);

% continuous function of desired q
q_des_c = fit([x,y],q_des,'linearinterp');
% continuous function of delta
delta_c = fit([x,y],delta,'linearinterp');

%% create PDE model
model = createpde();

geom_descrip = [3;4;0;1;1;0;0;0;1;1;];
dl = decsg(geom_descrip);
geometryFromEdges(model,dl);
pdegplot(model,'EdgeLabels','on');

applyBoundaryCondition(model,'dirichlet','Edge',1,'r',-115);
applyBoundaryCondition(model,'dirichlet','Edge',3,'r',115);
applyBoundaryCondition(model,'neumann','Edge',[2,4],'g',0);

conductance_handle = @(location,state) delta_c(location.x,location.y);
specifyCoefficients(model,'m',0,'d',0,'c',conductance_handle,'a',0,'f',0,'face',1);
generateMesh(model);

result = solvepde(model);
voltage = result.NodalSolution;
pdeplot(model,'XYData',voltage,'ZData',voltage);
title('voltage');

[e_x, e_y] = evaluateGradient(result,x,y);
q_cal = delta.*(e_x.^2+e_y.^2).^2*sigma;
err = q_cal-q_des;
err = err.^2;
err = sum(err);

%% Iterate
i = 0;

while (err > tol)
    i
    err
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
    
    delta_c = fit([x,y],delta,'linearinterp');
    
    conductance_handle = @(location,state) delta_c(location.x,location.y);
    specifyCoefficients(model,'m',0,'d',0,'c',conductance_handle,'a',0,'f',0,'face',1);
    generateMesh(model);

    result = solvepde(model);
    voltage = result.NodalSolution;
    pdeplot(model,'XYData',voltage,'ZData',voltage);
%     drawnow;
    title('voltage');
    
    [e_x, e_y] = evaluateGradient(result,x,y);
    q_cal = delta.*(e_x.^2+e_y.^2)*sigma;
    err = q_cal-q_des;
    err = err.^2;
    err = sum(err);
    i = i + 1;
end






