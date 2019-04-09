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

x_bottom = 0.0254*[0 3.622 7.411 11.107 14.724 18.279 21.772 25.214 28.611 32.037 35.463];
y_bottom = 0.0254*[0 0.815 3.101 5.532 8.081 10.715 13.431 16.211 19.046 22.421 25.795];

x_top = 0.0254*[0 2.712 5.378 7.994 10.566 13.09 15.567 18.005 20.404 22.764 25.09];
y_top = 0.0254*[20.902 22.085 23.37 24.754 26.217 27.76 29.38 31.057 32.79 34.576 36.404];


x_dat = zeros(11,11);
y_dat = zeros(11,11);

num_slope = y_top - y_bottom;
den_slope = x_top-x_bottom;

slope = num_slope./den_slope;

for i=1:11
    if slope(i) == inf
        slope(i) = 0;
    else
    end
end

b = y_bottom - slope.*x_bottom;

x_diff = x_bottom-x_top;

for i=1:11
    if x_diff(i) ~= 0
        x_dat(:,i) = x_top(i):x_diff(i)/10:x_bottom(i)';
    else
        x_dat(2:end-1,i) = x_top(i)*ones(9,1);
    end
    
    if slope~=0
        y_dat(:,i) = slope(i)*x(:,i) + b(i);
    else
        y_dat(:,i) = fliplr(y_bottom(i):(y_top(i)-y_bottom(i))/10:y_top(i))';
    end

end

geom_descrip = [2;22;x_bottom'; flip(x_top');y_bottom';flip(y_top')];
dl = decsg(geom_descrip);
geometryFromEdges(model,dl);
figure(1);
pdegplot(model,'EdgeLabels','on');

applyBoundaryCondition(model,'dirichlet','Edge',[1:1:10],'r',-115);
applyBoundaryCondition(model,'dirichlet','Edge',[12:1:21],'r',115);
applyBoundaryCondition(model,'neumann','Edge',[11,22],'g',0);

conductance_handle = @(location,state) delta_c(location.x,location.y);
specifyCoefficients(model,'m',0,'d',0,'c',conductance_handle,'a',0,'f',0,'face',1);
generateMesh(model);

result = solvepde(model);
voltage = result.NodalSolution;
figure(2);
pdeplot(model,'XYData',voltage,'ZData',voltage);
title('voltage');

[e_x, e_y] = evaluateGradient(result,x,y);
e_x_matrix = reshape(e_x', [11 11]);
e_y_matrix = reshape(e_y', [11 11]);

for j=1:11
    for k=1:11
        if isnan(e_x_matrix(j,k))
            if j == 1
                e_x_matrix(j,k) = e_x_matrix(j+1,k);
                e_y_matrix(j,k) = e_y_matrix(j+1,k);
            elseif j == 11
                e_x_matrix(j,k) = e_x_matrix(j-1,k);
                e_y_matrix(j,k) = e_y_matrix(j-1,k);
            else
                e_x_matrix(j,k) = .5*(e_x_matrix(j-1,k) + e_x_matrix(j+1,k));
                e_y_matrix(j,k) = .5*(e_y_matrix(j-1,k) + e_y_matrix(j+1,k));
            end
        end
    end
end
e_x = reshape(e_x_matrix, [121 1]);
e_y = reshape(e_y_matrix, [121 1]);

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
    e_x_matrix = reshape(e_x', [11 11]);
    e_y_matrix = reshape(e_y', [11 11]);

    for j=1:11
        for k=1:11
            if isnan(e_x_matrix(j,k))
                if j == 1
                    e_x_matrix(j,k) = e_x_matrix(j+1,k);
                    e_y_matrix(j,k) = e_y_matrix(j+1,k);
                elseif j == 11
                    e_x_matrix(j,k) = e_x_matrix(j-1,k);
                    e_y_matrix(j,k) = e_y_matrix(j-1,k);
                else
                    e_x_matrix(j,k) = .5*(e_x_matrix(j-1,k) + e_x_matrix(j+1,k));
                    e_y_matrix(j,k) = .5*(e_y_matrix(j-1,k) + e_y_matrix(j+1,k));
                end
            end
        end
    end
    e_x = reshape(e_x_matrix, [121 1]);
    e_y = reshape(e_y_matrix, [121 1]);
    q_cal = delta.*(e_x.^2+e_y.^2)*sigma;
    err = q_cal-q_des;
    err = err.^2;
    err = sum(err);
    i = i + 1;
end






