function [electron_properties] = electron_properties_with_mb_velocity(T, L, W, N)
    electron_properties = Electron_Properties_MB(N,T,L,W);
    
    electron_properties.x = uniform_random_gen(0,L,N);
    electron_properties.y = uniform_random_gen(0,W,N);
    
    electron_properties = get_maxwell_boltzmann_velocity_properties(electron_properties);
end

function electron_properties = get_maxwell_boltzmann_velocity_properties(electron_properties)
    [v_x,v_y,v_mag] = compute_maxwell_boltzmann_velocities(electron_properties);
    electron_properties.v_x = v_x;
    electron_properties.v_y = v_y;
    electron_properties.v_mag = v_mag;

    electron_properties.temperature = compute_electron_temperature(electron_properties);
end