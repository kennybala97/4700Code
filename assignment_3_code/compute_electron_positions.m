function [electrons] = compute_electron_positions(electrons, L, W, d_t)

    [electrons] = correct_out_of_bounds_particles_rethermalized(electrons, L, W, d_t);
    
    electrons.v_x = electrons.v_x + electrons.a_x*d_t;
    electrons.v_y = electrons.v_y + electrons.a_y*d_t;
    
    electrons.x = electrons.x + electrons.v_x*d_t;
    electrons.y = electrons.y + electrons.v_y*d_t;
end

function [electrons] = correct_out_of_bounds_particles(electrons, L, W, d_t)
    find_gt = @(a,b) (a > b);
    find_lt = @(a,b) (a < b);

    x_temp = electrons.x + (electrons.v_x + electrons.a_x*d_t)*d_t;
    y_temp = electrons.y + (electrons.v_y + electrons.a_y*d_t)*d_t;
    
    out_of_boundary_y = find_lt(y_temp,0) | find_gt(y_temp,W);
    
    neg_out_of_boundary_x = find_lt(x_temp,0);
    pos_out_of_boundary_x = find_gt(x_temp,L);
    
    electrons.v_y(out_of_boundary_y)  = -electrons.v_y(out_of_boundary_y);
    
    electrons.x(neg_out_of_boundary_x) = electrons.x(neg_out_of_boundary_x) + L;
    electrons.x(pos_out_of_boundary_x) = electrons.x(pos_out_of_boundary_x) - L;
end

function [electrons] = correct_out_of_bounds_particles_rethermalized(electrons, L, W, d_t)
    find_gt = @(a,b) (a > b);
    find_lt = @(a,b) (a < b);

    x_temp = electrons.x + (electrons.v_x + electrons.a_x*d_t)*d_t;
    y_temp = electrons.y + (electrons.v_y + electrons.a_y*d_t)*d_t;
    
    out_of_boundary_y = find_lt(y_temp,0) | find_gt(y_temp,W);
    
    neg_out_of_boundary_x = find_lt(x_temp,0);
    pos_out_of_boundary_x = find_gt(x_temp,L);
    
    electrons.v_y(out_of_boundary_y)  = -electrons.v_y(out_of_boundary_y);
    
    electrons.x(neg_out_of_boundary_x) = electrons.x(neg_out_of_boundary_x) + L;
    electrons.x(pos_out_of_boundary_x) = electrons.x(pos_out_of_boundary_x) - L;
    
    [v_x_new, v_y_new, v_mag_new] = compute_maxwell_boltzmann_velocities(electrons);

    electrons.v_x(pos_out_of_boundary_x) = abs(v_x_new(pos_out_of_boundary_x));
    electrons.v_y(pos_out_of_boundary_x) = v_y_new(pos_out_of_boundary_x);
    electrons.v_mag(pos_out_of_boundary_x) = v_mag_new(pos_out_of_boundary_x);
end

