function f = f_coeff(result,q_err_cont,e_x_cont,e_y_cont,x,y)
    q_diff_str = q_error_diff(q_err_cont,x,y);
    q_x = q_diff_str.x;
    q_y = q_diff_str.y;
    
    e_field_str = e_field_calc(result,e_x_cont,e_y_cont,x,y);
    e_x = e_field_str.x';
    e_y = e_field_str.y';
    e_xx = e_field_str.xx;
    e_yx = e_field_str.yx;
    e_yy = e_field_str.yy;
    nr = numel(x);
    f = zeros(1,nr);
%     f(1,:) = (-q_x.*e_x-q_y.*e_y)./(e_x.^2+e_y.^2);
    f(1,:) = -(q_x.*e_x+q_err_cont(x,y).*e_xx-q_err_cont(x,y).*e_x./(e_x.^2+e_y.^2).*(2*e_x.*e_xx+2*e_y.*e_yx) + ...
        q_y.*e_y+q_err_cont(x,y).*e_yy-q_err_cont(x,y).*e_y./(e_x.^2+e_y.^2).*(2*e_x.*e_yx+2*e_y.*e_yy))./(e_x.^2+e_y.^2);
end