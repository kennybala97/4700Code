clear all;
close all;
clc;

% This code generates electrons in a given distribution and then
% uses a spring model for force between electrons

n_t = 1000;
n = 0;

d_t = 1e-13;

q = -1.60217662e-19;
m_o = 9.10938215e-31;
r_o = 1e-9;

L = 10e-9;
W = 10e-9;

N = 500;

E = -1000;

% Here I define a variety of useful anonymous functions

compute_dist = @(d_x, d_y) sqrt( (d_x).^2 + (d_y).^2 );

compute_angle = @(x,y) atan2(y,x);

% Anonymous function to generate "N" uniformly distributed random numbers
% between a specified interval of "(a,b)"

uniform_random_position_gen = @(a,b,N,n) a + (b-a).*rand(N,n);

% Anonymous function to compute force by using the spring potential

force_spring_potential = @(r, r_o, k) (exp(k*(r)))*k*(r - r_o);

% Generate the electron positions.

electron_positions_x = zeros(N,1);
electron_positions_y = uniform_random_position_gen(-W/2,W/2,N,1);

electron_velocities_x = zeros(N,1);
electron_velocities_y = zeros(N,1);

velocity_to_plot = zeros(n_t,1);

% Compute the forces between each electron. This is done by using a spring
% potential model

ax_positions = subplot(1,2,1);
scatter_plot = scatter(ax_positions,electron_velocities_x, electron_velocities_y,'r.');
xlim([-L L]);
ylim([-W W]);

ax_velocity = subplot(1,2,2);
velocity_plot = plot(ax_velocity,velocity_to_plot);
hold(ax_velocity,'on');

while n < n_t
    for i=1:1:N
        x = electron_positions_x(i);
        y = electron_positions_y(i);

        v_x = electron_velocities_x(i);
        v_y = electron_velocities_y(i);            
        
        scattering_probability = rand();
        
        if scattering_probability > 0.95
           v_x = 0; 
        end
        
        F_x = q*E;

        a_x = F_x/m_o;

        v_x = v_x + a_x*d_t;

        x = x + v_x*d_t;
        
        if x > L
            x = x - 2*L;
        elseif x < -L
            x = x + 2*L;
        end

        if y > W
            y = y - 2*W;
        elseif x < -W
            y = y + 2*W;
        end
        
        electron_positions_x(i) = x;
        electron_positions_y(i) = y;
        electron_velocities_x(i) = v_x;
        electron_velocities_y(i) = v_y;
    end

    velocity_to_plot(n+1) = electron_velocities_x(1);
    
    set(scatter_plot,...
        'xData',electron_positions_x,...
        'yData',electron_positions_y);
    
    set(velocity_plot,'yData',velocity_to_plot);
    
    pause(0.1);
    
    n = n+1;
end



