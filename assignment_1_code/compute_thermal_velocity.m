function [v_th] = compute_thermal_velocity(electrons)
    v_th = sqrt( (2*electrons.k*electrons.T)/electrons.m_eff );
end