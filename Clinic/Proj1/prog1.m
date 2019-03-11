%% setup
clear;
close all;
clc;

%% data import and interp
global sigma; 
sigma = 1E6;
thickness_data = dlmread('thickness_init.csv',',',1,0);
heat_data = dlmread('heat_desired.csv',',',1,0);

x = thickness_data(:,1);
y = thickness_data(:,2);
delta_init = thickness_data(:,3);
q_des = heat_data(:,3);

q_des_cont = fit([x,y],q_des,'linearinterp');
delta_init_cont = fit([x,y],delta_init,'linearinterp');

%% create PDE model
model = createpde();

geom_descrip = [3;4;0;1;1;0;0;0;1;1;];
dl = decsg(geom_descrip);
geometryFromEdges(model,dl);
pdegplot(model,'EdgeLabels','on');

applyBoundaryCondition(model,'dirichlet','Edge',1,'r',-115);
applyBoundaryCondition(model,'dirichlet','Edge',3,'r',115);
applyBoundaryCondition(model,'neumann','Edge',[2,4],'g',0);

conductance_handle = @(location,state) delta_init_cont(location.x,location.y);
specifyCoefficients(model,'m',0,'d',0,'c',conductance_handle,'a',0,'f',0,'face',1);

generateMesh(model);

result = solvepde(model);
voltage = result.NodalSolution;
pdeplot(model,'XYData',voltage,'ZData',voltage);
title('voltage');

q_calc(result,.5,.5,delta_init_cont)


