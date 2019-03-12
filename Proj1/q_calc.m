function q = q_calc(result,x,y,delta)
    global sigma;
    e_field = e_field_calc(result,x,y);
    e_field_sq = e_field.x.^2 + e_field.y.^2;
    q = sigma*e_field_sq.*delta(x,y)';
end
