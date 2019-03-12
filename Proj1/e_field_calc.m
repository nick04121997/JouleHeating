function elec_field = e_field_calc(result,x,y)
    [e_field_x,e_field_y] = evaluateGradient(result,x,y);
end

