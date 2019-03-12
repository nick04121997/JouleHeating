function q_err_diff = q_error_diff(q_err_cont,x,y)
    [q_err_diff.x, q_err_diff.y] = differentiate(q_err_cont,x,y);
end

