function [electrons] = compute_electron_positions_w_pillbox(electrons, L, W, d_t, box)

    electrons.v_x = electrons.v_x + electrons.a_x*d_t;
    electrons.v_y = electrons.v_y + electrons.a_y*d_t;

    x_temp = electrons.x + electrons.v_x*d_t;
    y_temp = electrons.y + electrons.v_y*d_t;
    
    [electrons] = correct_out_of_region(electrons, x_temp, y_temp, L, W);

    
    inds_in_upper_box = in_upper_box(x_temp, y_temp, L, W, box);
    inds_in_lower_box = in_lower_box(x_temp, y_temp, L, W, box);
    
    electrons = correct_box_side_collision(electrons, box, L, W, inds_in_upper_box, inds_in_lower_box);

    electrons.x = electrons.x + electrons.v_x*d_t;
    electrons.y = electrons.y + electrons.v_y*d_t;

end


function [electrons] = correct_out_of_region(electrons, x_temp, y_temp, L, W)
    find_gt = @(a,b) (a > b);
    find_lt = @(a,b) (a < b);
    
    out_of_boundary_y = find_lt(y_temp,0) | find_gt(y_temp,W);
    
    neg_out_of_boundary_x = find_lt(x_temp,0);
    pos_out_of_boundary_x = find_gt(x_temp,L);
    
    electrons.v_y(out_of_boundary_y)  = -electrons.v_y(out_of_boundary_y);
    
    electrons.x(neg_out_of_boundary_x) = electrons.x(neg_out_of_boundary_x) + L;
    electrons.x(pos_out_of_boundary_x) = electrons.x(pos_out_of_boundary_x) - L;

end

function [inds_in_upper_box] = in_upper_box(x_temp, y_temp, L, W, box)
    x = x_temp;
    y = y_temp;
    
    inds_in_upper_box = (x >= (L - box.L)/2) & (x <= (L + box.L)/2) & (y > (W + box.gap)/2);
end

function [inds_in_lower_box] = in_lower_box(x_temp, y_temp, L, W, box)
    x = x_temp;
    y = y_temp;

    inds_in_lower_box = (x >= (L - box.L)/2) & (x <= (L + box.L)/2) & (y < (W - box.gap)/2);
end

function [electrons] = correct_box_side_collision(electrons, box, L, W, in_upper_box, in_lower_box)
    x = electrons.x;
    y = electrons.y;

    in_box = (in_upper_box | in_lower_box);
    
    find_rh = x < (L - box.L)/2;
    find_lh = x > (L + box.gap)/2;
    find_m = (x > (L - box.L)/2) & (x < (L + box.L)/2) & (y > (W - box.gap)/2)  & (y < (W + box.gap)/2);

    electrons.v_x(in_box & (find_rh | find_lh)) = -electrons.v_x(in_box & (find_rh | find_lh));
    electrons.v_y(in_box & find_m) = -electrons.v_y(in_box & find_m);
end


