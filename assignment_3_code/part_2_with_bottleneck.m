function [bottlenecks,I_av] = part_2_with_bottleneck()
    clear all;
    close all;
    clc;
    
%     [final_x, final_y, ~, ~] = simulate(10000,1);
    
%     plot_positions()
    [bottlenecks,I_av] = plot_current();
    
end

function [] = plot_positions()

    figure;
    [final_x, final_y, ~, ~] = simulate(1000,100);
    
    for i=1:10:1000
       scatter(final_x(i,:)*1e9,final_y(i,:)*1e9,'.');hold on;
    end
    title('Part 2 : Position Over Time for Different Particles, $V_o=100V$', 'interpreter', 'latex');
    xlabel('X(nm)', 'interpreter', 'latex');
    ylabel('Y(nm)', 'interpreter', 'latex');
    grid on;

end

function [bottlenecks,I_av] = plot_current()
    figure;

    bottlenecks = 0.1:0.1:10;
    
    for i = 1:length(bottlenecks)
        
        [~, ~, I, d_t] = simulate(1000, bottlenecks(i));
        I_av(i) = mean(I);
    end
    
    scatter(bottlenecks,I_av,'.');
    
    title('Average Drift Current $I_{d}$ for Different Bottlenecks', 'interpreter', 'latex');
    xlabel('Bottleneck Width(nm)', 'interpreter', 'latex');
    ylabel('Current(mA)', 'interpreter', 'latex');
    grid on;

end

function [final_x, final_y, I, d_t] = simulate(N, V0)

    addpath("./laplace_solver");
    
    conc = (10^15)*(100^2);

    T = 300;
    L = 200e-9;
    W = 100e-9;

    tau_mn = 0.2e-12;

    box.L = 40e-9;
    box.gap = 20e-9;
%     box.gap = bottleneck;
    
    nx = 100;
    ny = 50;
    
    delta = L/nx;
    
    electron_properties = electron_properties_with_mb_velocity_pillbox(T,L,W,N,box);
    
    [V, E_x, E_y, J_x, J_y, ~, ~] = part_2_solver(nx,ny,delta,0.01,1,box.L*1e9,box.gap*1e9,V0);
    
    [X, Y] = ndgrid((0:1:(nx - 1))*delta, (0:1:(ny - 1))*delta);
        
    electron_properties = calculate_a_e_field(electron_properties, X, Y, E_x, E_y);

    d_t = (W/1000)/electron_properties.v_th;
    n = 1;
    n_final = 200;
    p_scat = 1-exp(-d_t/tau_mn);

    % f = figure;

    v_dx_av = zeros(1,n_final);

    final_x = zeros(N,n_final);
    final_y = zeros(N,n_final);

%     f = figure;
    
    [xq,yq] = meshgrid(0:1e-9:L, 0:1e-9:W);
    
    while n < n_final
        electron_properties = calculate_a_e_field(electron_properties, X, Y, E_x, E_y);

        [electron_properties] = compute_electron_positions_w_pillbox(electron_properties, L, W, d_t, box);
        electron_properties.temperature = compute_electron_temperature(electron_properties);

        scattering_electron_indices = p_scat > rand(N,1);

        [v_x_new, v_y_new, v_mag_new] = compute_maxwell_boltzmann_velocities(electron_properties);

        v_dx_av(n) = mean(electron_properties.v_x);
%     
%         subplot(211);
%         scatter(electron_properties.x, electron_properties.y,10,'.k');hold on;
%         
%         quiver(electron_properties.x, electron_properties.y, electron_properties.v_x, electron_properties.v_y,0.1,'r');
%         hold off;
% 
%         xlim([0 L]);
%         ylim([0 W]); 
%         title('Position and Velocity for Part 3 Simulation, $V_o=1V$', 'interpreter', 'latex');
%         xlabel('X(m)', 'interpreter', 'latex');
%         ylabel('Y(m)', 'interpreter', 'latex');
%         
%         subplot(212);
% 
%         heatmap = griddata(electron_properties.x, electron_properties.y, electron_properties.temperature,xq,yq);
%         
%         heatmap(xq > (L - box.L)/2 & xq < (L + box.L)/2 & yq > (W + box.gap)/2) = 0;
%         heatmap(xq > (L - box.L)/2 & xq < (L + box.L)/2 & yq < (W - box.gap)/2) = 0;
%         
%         pcolor(xq, yq, heatmap);
%         title('Temperature Map(K) For Part 3 Simulation, $V_o=1V$', 'interpreter', 'latex');
%         xlabel('X(m)', 'interpreter', 'latex');
%         ylabel('Y(m)', 'interpreter', 'latex');
%         
%         view(2)
% 
%         colormap jet;
% 
%         colorbar;
%         caxis([0,1000]);
% 
%         pause(0.01);
%         drawnow
%         frame = getframe(f);
% 
%         im{n} = frame2im(frame);
        
        electron_properties.v_x(scattering_electron_indices) = v_x_new(scattering_electron_indices);
        electron_properties.v_y(scattering_electron_indices) = v_y_new(scattering_electron_indices);
        electron_properties.v_mag(scattering_electron_indices) = v_mag_new(scattering_electron_indices);

        final_x(:,n) = electron_properties.x;
        final_y(:,n) = electron_properties.y;
        
%         pause(0.01);
        
        n = n+1;
    end
    
    I = abs(electron_properties.q)*conc*v_dx_av*W;        
        
%     filename = 'part_2_1V.gif'; % Specify the output file name
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


function [electron_properties] = calculate_a_e_field(electron_properties, X, Y, E_x, E_y)

    f_E_x = griddedInterpolant(X,Y,E_x);
    f_E_y = griddedInterpolant(X,Y,E_y);

    E_x_q = f_E_x(electron_properties.x,electron_properties.y);
    E_y_q = f_E_y(electron_properties.x,electron_properties.y);
    
    force_e_field_x = electron_properties.q*E_x_q;
    force_e_field_y = electron_properties.q*E_y_q;
    
    electron_properties.a_y = force_e_field_x/electron_properties.m_eff;
    electron_properties.a_x = force_e_field_y/electron_properties.m_eff;
end