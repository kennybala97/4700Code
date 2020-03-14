clear all;
close all;
clc;

% This code generates electrons in a given distribution and then
% uses a spring model for force between electrons

n_t = 1000;
n = 0;

d_t = 1e-13;

m_o = 9.10938215e-31;
r_o = 1e-9;

L = 10e-9;
W = 10e-9;

N = 500;

k = -0.0000000002;

% Here I define a variety of useful anonymous functions

compute_dist = @(d_x, d_y) sqrt( (d_x).^2 + (d_y).^2 );

compute_angle = @(x,y) atan2(y,x);

% Anonymous function to generate "N" uniformly distributed random numbers
% between a specified interval of "(a,b)"

uniform_random_position_gen = @(a,b,N) a + (b-a).*rand(N,1);

% Anonymous function to compute force by using the spring potential

force_spring_potential = @(r, r_o, k) (exp(k*(r)))*k*(r - r_o);

% Generate the electron positions.

electron_positions_x = uniform_random_position_gen(-L/2,L/2,N);
electron_positions_y = uniform_random_position_gen(-W/2,W/2,N);

electron_velocities_x = zeros(N,1);
electron_velocities_y = zeros(N,1);

% Compute the forces between each electron. This is done by using a spring
% potential model

f = figure;

while n < n_t
    for i=1:1:N
        x = electron_positions_x(i);
        y = electron_positions_y(i);

        v_x = electron_velocities_x(i);
        v_y = electron_velocities_y(i);

        for j=1:1:N
            if j == i
                continue
            end
            
            d_x = x - electron_positions_x(j);
            d_y = y - electron_positions_y(j);

            r = compute_dist(d_x, d_y);
            angle = compute_angle(d_x,d_y);

            F = force_spring_potential(r,r_o,k);

            F_x = F*cos(angle);
            F_y = F*sin(angle);

            a_x = F_x/m_o;
            a_y = F_y/m_o;

            v_x = v_x + a_x*d_t;
            v_y = v_y + a_y*d_t;
        end
        x = x + v_x*d_t;
        y = y + v_y*d_t;
        
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

    scatter(electron_positions_x, electron_positions_y, '.');
    ylim([-W W]);
    xlim([-L L]);
    drawnow;
    frame = getframe(f);
    im = frame2im(frame);
    
    [imind,cm] = rgb2ind(im,256); 
    % Write to the GIF File 
    if n == 0 
      imwrite(imind,cm,'md_spring.gif','gif', 'Loopcount',inf); 
    else 
      imwrite(imind,cm,'md_spring.gif','gif','WriteMode','append'); 
    end 
    
    n = n+1;
end



