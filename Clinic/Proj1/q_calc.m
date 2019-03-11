function q = q_calc(result,x,y,delta)
    global sigma;
    [e_field_x, e_field_y] = e_field(result,x,y);
    e_field_sq = e_field_x^2 + e_field_y^2;
    q = sigma*e_field_sq*delta(x,y);
    sigma
end
