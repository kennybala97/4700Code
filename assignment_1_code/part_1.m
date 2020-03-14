clear all;
close all;
clc;

N = 1000;
T = 300;
L = 200e-9;
W = 100e-9;

electrons = electron_properties_with_uniform_velocity(T,L,W,N);

d_t = (W/100)/electrons.v_th;
n = 1;
n_final = 1000;

f = figure;

F(n_final) = struct('cdata',[],'colormap',[]);

while n < n_final
    [electrons] = compute_electron_positions(electrons, L, W, d_t);
    electrons.temperature = compute_electron_temperature(electrons);
        
    scatter(electrons.x, electrons.y,10,'.r');hold on;
    quiver(electrons.x, electrons.y, electrons.v_x, electrons.v_y,0.1,'b');
    
    hold off;

    xlim([-L L]);
    ylim([-W W]);   
    
    drawnow
    F(n) = getframe(f);
    
    n = n+1;
end

fig = figure;
movie(fig,F,1);
