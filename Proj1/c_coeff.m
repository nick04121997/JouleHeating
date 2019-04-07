function c = c_coeff(result,delta,x,y,e_x_cont,e_y_cont)
    e_field_str = e_field_calc(result,e_x_cont,e_y_cont,x,y);
    e_x = e_field_str.x;
    e_y = e_field_str.y;
    c_11 = 1-2*e_x.^2./(e_x.^2+e_y.^2);
    c_12 = -2*e_x.*e_y./(e_x.^2+e_y.^2);
    c_22 = 1-2*e_y.^2./(e_x.^2+e_y.^2);
    nr = numel(x);
    c = zeros(3,nr);
    c(1,:) = delta(x,y)*c_11;
    c(2,:) = delta(x,y)*c_12;
    c(3,:) = delta(x,y)*c_22;
end