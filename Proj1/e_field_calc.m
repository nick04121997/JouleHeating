function e_field = e_field_calc(result,e_x_cont,e_y_cont,x,y)
    [e_field_x,e_field_y] = evaluateGradient(result,x,y);
    e_field.x = e_field_x;
    e_field.y = e_field_y;
    
    e_xx = differentiate(e_x_cont,x,y);
    [e_yx, e_yy] = differentiate(e_y_cont,x,y);
    
    e_field.xx = e_xx;
    e_field.yx = e_yx;
    e_field.yy = e_yy;
end
    

