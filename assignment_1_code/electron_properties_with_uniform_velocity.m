function [electron_properties] = electron_properties_with_uniform_velocity(T, L, W, N)
    electron_properties = Electron_Properties(N,T,L,W);
    
    electron_properties.x = uniform_random_gen(0,L,N);
    electron_properties.y = uniform_random_gen(0,W,N);
    
    electron_properties = get_uniform_velocity_properties(electron_properties);
end

function electron_properties = get_uniform_velocity_properties(electron_properties)
    electron_properties.theta = uniform_random_gen(-pi,pi,electron_properties.N);

    electron_properties.v_x = electron_properties.v_th.*cos(electron_properties.theta);
    electron_properties.v_y = electron_properties.v_th.*sin(electron_properties.theta);

    electron_properties.temperature = compute_electron_temperature(electron_properties);
end