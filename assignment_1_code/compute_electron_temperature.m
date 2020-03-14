function [temperature] = compute_electron_temperature(electron_properties)
    velocity_magnitude = sqrt(electron_properties.v_x.^2 + electron_properties.v_y.^2);
    temperature = (velocity_magnitude.^2*electron_properties.m_eff)/(3*electron_properties.k);
end