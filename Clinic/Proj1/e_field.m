function [e_field_x,e_field_y] = e_field(result,x,y)
    [e_field_x,e_field_y] = evaluateGradient(result,x,y);
end

