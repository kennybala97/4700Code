function [] = part_1_no_bottleneck()
    clear all;
    close all;
    clc;
    
    simulate_part_d(0.1, 0, 10000)
%     plot_positions()
%     plot_current()
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

    voltages = 0:0.3:1.2;
    
    for i = 1:length(voltages)
        
        [~, ~, I, d_t] = simulate(voltages(i), 0, 500);
        
        scatter((1:1:1000)*d_t*1e12,I*100,'.');
        hold on;
        
        legend_arr(i) = sprintf("$V_o = %5.2f$",voltages(i));
    end
    
    title('Average Drift Current $I_{d}$ for Different Voltages', 'interpreter', 'latex');
    xlabel('Time(ps)', 'interpreter', 'latex');
    ylabel('Current(mA)', 'interpreter', 'latex');
    legend(legend_arr,"interpreter", "latex");
    grid on;

end

function [final_x, final_y, I, d_t] = simulate(V_INIT, V_DIR, N)

    conc = (10^15)*(100^2);

    T = 300;
    L = 200e-9;
    W = 100e-9;

    tau_mn = 0.2e-12;

    box.L = 40e-9;
    box.gap = 20e-9;


    electron_properties = electron_properties_with_mb_velocity(T,L,W,N);

    e_field = V_INIT/L;
    force_e_field = electron_properties.q*e_field;

    a_e_field = force_e_field/electron_properties.m_eff;

    electron_properties.a_x = a_e_field*cosd(V_DIR);
    electron_properties.a_y = a_e_field*sind(V_DIR);

    d_t = (W/10)/electron_properties.v_th;
    n = 1;
    n_final = 50;
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

        n = n+1;        
    end
    
    I = abs(electron_properties.q)*conc*v_dx_av*W;
end


function [] = simulate_part_d(V_INIT, V_DIR, N)

    T = 300;
    L = 200e-9;
    W = 100e-9;

    tau_mn = 0.2e-14;

    electron_properties = electron_properties_with_mb_velocity(T,L,W,N);

    e_field = V_INIT/L;
    force_e_field = electron_properties.q*e_field;

    a_e_field = force_e_field/electron_properties.m_eff;

    electron_properties.a_x = a_e_field*cosd(V_DIR);
    electron_properties.a_y = a_e_field*sind(V_DIR);

    d_t = (W/1000)/electron_properties.v_th;
    n = 1;
    n_final = 100;
    p_scat = 1-exp(-d_t/tau_mn);

    f = figure('units','normalized','outerposition',[0 0 1 1]);

    [xq,yq] = meshgrid(0:1e-9:L, 0:1e-9:W);

    while n < n_final
        [electron_properties] = compute_electron_positions(electron_properties, L, W, d_t);

        electron_properties.temperature = compute_electron_temperature(electron_properties);

        scattering_electron_indices = p_scat > rand(N,1);

        [v_x_new, v_y_new, v_mag_new] = compute_maxwell_boltzmann_velocities(electron_properties);

        v_dx_av(n) = mean(electron_properties.v_x);

        electron_properties.v_x(scattering_electron_indices) = v_x_new(scattering_electron_indices);
        electron_properties.v_y(scattering_electron_indices) = v_y_new(scattering_electron_indices);
        electron_properties.v_mag(scattering_electron_indices) = v_mag_new(scattering_electron_indices);

        subplot(211);
        scatter(electron_properties.x, electron_properties.y,10,'.k');hold on;
        
        quiver(electron_properties.x, electron_properties.y, electron_properties.v_x, electron_properties.v_y,0.1,'r');
        hold off;

        xlim([0 L]);
        ylim([0 W]); 
        title('Position and Velocity for Part 1 Simulation, $V_o=0.1V$', 'interpreter', 'latex');
        xlabel('X(m)', 'interpreter', 'latex');
        ylabel('Y(m)', 'interpreter', 'latex');
        
        subplot(212);

        heatmap = griddata(electron_properties.x, electron_properties.y, electron_properties.temperature,xq,yq);


        pcolor(xq, yq, heatmap);
        title('Temperature Map(K) For Part 1 Simulation, $V_o=0.1V$', 'interpreter', 'latex');
        xlabel('X(m)', 'interpreter', 'latex');
        ylabel('Y(m)', 'interpreter', 'latex');
        
        view(2)

        colormap jet;

        colorbar;
        caxis([0,1000]);

        
%         subplot(133);
%         compass(mean(electron_properties.v_x),mean(electron_properties.v_y));
%         title('Compass Plot of Velocity Direction, $V_o=100V$', 'interpreter', 'latex');
        pause(0.01);
%         drawnow
%         frame = getframe(f);
% 
%         im{n} = frame2im(frame);
        
        n = n+1;
        
    end
        
%     filename = 'part_1_100V.gif'; % Specify the output file name
% 
%     n_images = length(im);
% 
%     for idx = 1:n_images
%         [A,map] = rgb2ind(im{idx},256);
%         if idx == 1
%             imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0.01);
%         else
%             imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.01);
%         end
%     end
end

