function [v_x,v_y,v_mag] = compute_maxwell_boltzmann_velocities(electron_properties)
    variance = electron_properties.k*electron_properties.T/electron_properties.m_eff;

    v_x = randn(electron_properties.N,1)*sqrt(variance);
    v_y = randn(electron_properties.N,1)*sqrt(variance);

    v_mag = sqrt(electron_properties.v_x.^2 + electron_properties.v_y.^2);
end

