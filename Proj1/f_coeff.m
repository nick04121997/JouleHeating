function f = f_coeff(result,q_err_cont,x,y)
    q_diff_str = q_error_diff(q_err_cont,x,y);
    q_x = q_diff_str.x;
    q_y = q_diff_str.y;
    e_field_str = e_field_calc(result,x,y);
    e_x = e_field_str.x';
    e_y = e_field_str.y';
    nr = numel(x)
    f = zeros(1,nr);
    f(1,:) = (-q_x.*e_x-q_y.*e_y)./(e_x.^2+e_y.^2);
end