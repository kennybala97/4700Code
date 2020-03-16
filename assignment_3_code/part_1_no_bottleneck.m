function [] = part_1_no_bottleneck()
    clear all;
    close all;
    clc;
    
    plot_positions()
    plot_current()
end

function [] = plot_positions()

    figure;
    [final_x, final_y, ~, ~] = simulate(0.1, 0, 5);
    
    for i=1:1:5
       scatter(final_x(i,:)*1e9,final_y(i,:)*1e9,'.');hold on;
    end
    title('Position Over Time for Five Different Particles, $V_o=0.1V$', 'interpreter', 'latex');
    xlabel('X(nm)', 'interpreter', 'latex');
    ylabel('Y(nm)', 'interpreter', 'latex');
    grid on;

end

function [] = plot_current()
    figure;

    marks = ["r.","b.","g."];
    voltages = 0:0.2:0.4;
    
    for i = 1:3
        
        [~, ~, I, d_t] = simulate(voltages(i), 0, 500);
        
        scatter((1:1:1000)*d_t*1e9,I,marks(i));
        hold on;
    end
    
    title('Average Drift Current $v_{dx}$ for Different Voltages', 'interpreter', 'latex');
    xlabel('Time(ns)', 'interpreter', 'latex');
    ylabel('Current(A)', 'interpreter', 'latex');
    legend("$V_o = 0V$","$V_o = 0.2V$","$V_o = 0.4V$","interpreter", "latex");
    grid on;

end

function [final_x, final_y, I, d_t] = simulate(V_INIT, V_DIR, N)

    conc = 10^15;

    T = 300;
    L = 200e-9;
    W = 100e-9;

    tau_mn = 0.2e-11;

    box.L = 40e-9;
    box.gap = 20e-9;


    electron_properties = electron_properties_with_mb_velocity(T,L,W,N);

    e_field = V_INIT/L;
    force_e_field = electron_properties.q*e_field;

    a_e_field = force_e_field/electron_properties.m_eff;

    electron_properties.a_x = a_e_field*cosd(V_DIR);
    electron_properties.a_y = a_e_field*sind(V_DIR);

    d_t = (W/1000)/electron_properties.v_th;
    n = 1;
    n_final = 1000;
    p_scat = 1-exp(-d_t/tau_mn);

    % f = figure;

    [xq,yq] = meshgrid(0:1e-9:L, 0:1e-9:W);

    v_dx_av = zeros(1,n_final);

    final_x = zeros(N,n_final);
    final_y = zeros(N,n_final);

    while n < n_final
        [electron_properties] = compute_electron_positions(electron_properties, L, W, d_t);

        electron_properties.temperature = compute_electron_temperature(electron_properties);

        scattering_electron_indices = p_scat > rand(N,1);

        [v_x_new, v_y_new, v_mag_new] = compute_maxwell_boltzmann_velocities(electron_properties);

        v_dx_av(n) = mean(electron_properties.v_x);

        electron_properties.v_x(scattering_electron_indices) = v_x_new(scattering_electron_indices);
        electron_properties.v_y(scattering_electron_indices) = v_y_new(scattering_electron_indices);
        electron_properties.v_mag(scattering_electron_indices) = v_mag_new(scattering_electron_indices);

        final_x(:,n) = electron_properties.x;
        final_y(:,n) = electron_properties.y;

    %     scatter(electron_properties.x, electron_properties.y,10,'.k');hold on;
    %     quiver(electron_properties.x, electron_properties.y, electron_properties.v_x, electron_properties.v_y,0.1,'r');
    %     hold off;

    %     xlim([0 L]);
    %     ylim([0 W]); 
    % 
    %     pause(0.01);
    %     drawnow
    %     frame = getframe(f);
    % 
    %     im{n} = frame2im(frame);


        n = n+1;
    end
    
    I = abs(electron_properties.q)*conc*v_dx_av;
end
