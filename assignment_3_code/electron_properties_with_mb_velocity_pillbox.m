function [electron_properties] = electron_properties_with_mb_velocity_pillbox(T, L, W, N, box)
    electron_properties = Electron_Properties_MB(N,T,L,W);
    
    electron_properties.x = uniform_random_gen(0,L,N);
    electron_properties.y = uniform_random_gen(0,W,N);
    
    [inds_in_upper_box] = in_upper_box(electron_properties, L, W, box);
    [inds_in_lower_box] = in_lower_box(electron_properties, L, W, box);
    
    while sum(inds_in_upper_box + inds_in_lower_box) > 0
        num_remaining = sum(inds_in_upper_box | inds_in_lower_box);
        
        electron_properties.x(inds_in_upper_box | inds_in_lower_box) = uniform_random_gen(0,L,num_remaining);
        electron_properties.y(inds_in_upper_box | inds_in_lower_box) = uniform_random_gen(0,W,num_remaining);
        
        [inds_in_upper_box] = in_upper_box(electron_properties, L, W, box);
        [inds_in_lower_box] = in_lower_box(electron_properties, L, W, box);
    end
    
    electron_properties = get_maxwell_boltzmann_velocity_properties(electron_properties);
end

function electron_properties = get_maxwell_boltzmann_velocity_properties(electron_properties)
    [v_x,v_y,v_mag] = compute_maxwell_boltzmann_velocities(electron_properties);
    electron_properties.v_x = v_x;
    electron_properties.v_y = v_y;
    electron_properties.v_mag = v_mag;

    electron_properties.temperature = compute_electron_temperature(electron_properties);
end

function [inds_in_upper_box] = in_upper_box(electrons, L, W, box)
    x = electrons.x;
    y = electrons.y;
    
    inds_in_upper_box = (x >= (L - box.L)/2) & (x <= (L + box.L)/2) & (y > (W + box.gap)/2);
end

function [inds_in_lower_box] = in_lower_box(electrons, L, W, box)
    x = electrons.x;
    y = electrons.y;

    inds_in_lower_box = x >= (L - box.L)/2 & x <= (L + box.L)/2 & y < (W - box.gap)/2;
end
