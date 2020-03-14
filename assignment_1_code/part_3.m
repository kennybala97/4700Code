clear all;
close all;
clc;

N = 100000;
T = 300;
L = 200e-9;
W = 100e-9;

tau_mn = 0.2e-10;

box.L = 40e-9;
box.gap = 20e-9;

electron_properties = electron_properties_with_mb_velocity_pillbox(T,L,W,N, box);

d_t = (W/100)/electron_properties.v_th;
n = 1;
n_final = 50;
p_scat = 1-exp(-d_t*0/tau_mn);

f = figure;

[xq,yq] = meshgrid(0:1e-9:L, 0:1e-9:W);

while n < n_final
    [electron_properties] = compute_electron_positions_w_pillbox(electron_properties, L, W, d_t, box);

    electron_properties.temperature = compute_electron_temperature(electron_properties);

    scattering_electron_indices = p_scat > rand(N,1);
    
    [v_x_new, v_y_new, v_mag_new] = compute_maxwell_boltzmann_velocities(electron_properties);
    
    subplot(2,1,1);
    
    electron_properties.v_x(scattering_electron_indices) = v_x_new(scattering_electron_indices);
    electron_properties.v_y(scattering_electron_indices) = v_y_new(scattering_electron_indices);
    electron_properties.v_mag(scattering_electron_indices) = v_mag_new(scattering_electron_indices);
    
    scatter(electron_properties.x, electron_properties.y,10,'.k');hold on;
    quiver(electron_properties.x, electron_properties.y, electron_properties.v_x, electron_properties.v_y,0.1,'r');
    hold off;

    xlim([0 L]);
    ylim([0 W]); 
    
    subplot(2,1,2);
    
    heatmap = griddata(electron_properties.x, electron_properties.y, electron_properties.temperature,xq,yq);

    heatmap(xq > (L - box.L)/2 & xq < (L + box.L)/2 & yq > (W + box.gap)/2) = 0;
    heatmap(xq > (L - box.L)/2 & xq < (L + box.L)/2 & yq < (W - box.gap)/2) = 0;
    

    pcolor(heatmap);
    view(2)
    
    colormap jet;
    
    colorbar;

  
    pause(0.01);
    drawnow
    frame = getframe(f);
    
    im{n} = frame2im(frame);
    
    n = n+1;
end

filename = 'part_3.gif'; % Specify the output file name

n_images = length(im);

for idx = 1:n_images
    [A,map] = rgb2ind(im{idx},256);
    if idx == 1
        imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0.01);
    else
        imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.01);
    end
end