function e_field = e_field_calc(result,x,y)
    [e_field_x,e_field_y] = evaluateGradient(result,x,y);
    e_field.x = e_field_x;
    e_field.y = e_field_y;
end

