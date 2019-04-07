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

% q bar
q_des_cont = fit([x,y],q_des,'linearinterp');
% delta bar
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

[e_x, e_y] = evaluateGradient(result,x,y);
e_x_cont = fit([x,y],e_x,'cubicinterp');
e_y_cont = fit([x,y],e_y,'cubicinterp');

% q calcuated on that iteration 
q_n = q_calc(result,x',y',delta_init_cont,e_x_cont,e_y_cont);
% q tilde
q_err = q_n - q_des;
q_err_cont = fit([x,y],q_err,'linearinterp');

%% solve the voltage error PDE (v tilde)
f = @(location,state) f_coeff(result,q_err_cont,e_x_cont,e_y_cont,location.x,location.y);
c = @(location,state) c_coeff(result, delta_init_cont, location.x, location.y,e_x_cont,e_y_cont);

model_v = createpde();

geometryFromEdges(model_v,dl);
figure(2);
pdegplot(model,'EdgeLabels','on');

applyBoundaryCondition(model,'dirichlet','Edge',1,'r',0);
applyBoundaryCondition(model,'dirichlet','Edge',3,'r',0);
applyBoundaryCondition(model,'neumann','Edge',[2,4],'g',0);

specifyCoefficients(model_v,'m',0,'d',0,'c',c,'a',0,'f',f,'face',1);

generateMesh(model_v);

result_v = solvepde(model_v);


