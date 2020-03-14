clear all;
close all;
clc;

N = 10;
T = 300;
L = 200e-9;
W = 100e-9;

tau_mn = 0.2e-12;

electron_properties = electron_properties_with_mb_velocity(T,L,W,N);

d_t = (W/100)/electron_properties.v_th;
n = 1;
n_final = 1000;
p_scat = 1-exp(-d_t/tau_mn);

f = figure;

F(n_final) = struct('cdata',[],'colormap',[]);

while n < n_final
    [electron_properties] = compute_electron_positions(electron_properties, L, W, d_t);

    electron_properties.temperature = compute_electron_temperature(electron_properties);
    
    scattering_electron_indices = p_scat > rand(N,1);
    
    [v_x_new, v_y_new, v_mag_new] = compute_maxwell_boltzmann_velocities(electron_properties);
    
    electron_properties.v_x(scattering_electron_indices) = v_x_new(scattering_electron_indices);
    electron_properties.v_y(scattering_electron_indices) = v_y_new(scattering_electron_indices);
    electron_properties.v_mag(scattering_electron_indices) = v_mag_new(scattering_electron_indices);
    
    scatter(electron_properties.x, electron_properties.y,10,'.k');hold on;
%     quiver(electron_properties.x, electron_properties.y, electron_properties.v_x, electron_properties.v_y,0.1,'r');

%     hold off;

    xlim([0 L]);
    ylim([0 W]);   
    pause(0.01);
    drawnow
    F(n) = getframe(f);
    
    n = n+1;
end

