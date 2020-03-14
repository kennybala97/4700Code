classdef Electron_Properties
    
    properties(Constant)
        m_eff = 9.10938215e-31 * 0.26;
        k = 1.38064852e-23;
    end
    
    properties
        N;
        T;
        
        x_len;
        y_len;
        
        v_th;
    
        x;
        y;

        theta;

        v_x;
        v_y;

        temperature;
    end
    
    methods
        
        function electron_properties = Electron_Properties(N, T, x_len, y_len)
            electron_properties.N = N;
            electron_properties.T = T;
            electron_properties.x_len = x_len;
            electron_properties.y_len = y_len;
            
            electron_properties.v_th = compute_thermal_velocity(electron_properties);
        end
    
        function [v_th] = compute_thermal_velocity(electron_properties)
            v_th = sqrt( ( 2 * electron_properties.k * electron_properties.T ) / electron_properties.m_eff );
        end
        
    end
    
end

