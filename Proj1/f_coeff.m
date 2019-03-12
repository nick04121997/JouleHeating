function f = f_coeff(result,q_err_cont,x,y)
    q_x = q_error_diff(q_err_cont,x,y).x;
    q_y = q_error_diff(q_err_cont,x,y).y;
    e_x = e_field_calc(result,x,y).x;
    e_y = e_field_calc(result,x,y).y;
    f = (-q_x*e_x-q_y*e_y)/(e_x^2+e_y^2);
end