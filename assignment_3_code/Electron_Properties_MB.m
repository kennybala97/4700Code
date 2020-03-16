classdef Electron_Properties_MB < Electron_Properties
    properties
        a;
        v_mag;
    end
    
    methods
        function mb_electron_properties = Electron_Properties_MB(N, T, x_len, y_len)
            mb_electron_properties = mb_electron_properties@Electron_Properties(N, T, x_len, y_len);
            mb_electron_properties.a = sqrt(mb_electron_properties.k*mb_electron_properties.T/mb_electron_properties.m_eff);
        end
        
    end
end

