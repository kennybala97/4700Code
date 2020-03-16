clear all;
close all;
clc;

nx = 150;
ny = 100;
delta = 1;

a_cond = 1e-2;
b_cond = 1;

L = nx*delta;
W = ny*delta;

V0 = 100;

box_w = 50;
box_gap = 50;

x_c = nx/2;
y_c = ny/2;

[X,Y] = meshgrid(1:1:nx,1:1:ny);
X = X';
Y = Y';

c_map = ones(nx,ny)*b_cond;

box_top_idx = (X > x_c - box_w/2) & (X < x_c + box_w/2) & (Y > y_c + box_gap/2);
box_bottom_idx = (X > x_c - box_w/2) & (X < x_c + box_w/2) & (Y < y_c - box_gap/2);

c_map(box_top_idx | box_bottom_idx) = a_cond;

[V_numeric,E_x,E_y,J_x,J_y] = laplace_solver_2d_part_2(nx,ny,delta,c_map);
 
figure(1);
subplot(2,1,1);
pcolor(c_map);
xlabel('x(cm)');
ylabel('y(cm)');    
h = colorbar;
set(get(h,'title'),'string','\sigma(S)');
title('Assignment 2 Part 2 A - Conductivity \sigma(x,y)');
colormap jet;

subplot(2,1,2);
pcolor(V_numeric);
xlabel('x(cm)');
ylabel('y(cm)');    
h = colorbar;
set(get(h,'title'),'string','Voltage(V)');
title('Assignment 2 Part 2 A - Voltage');
colormap jet;
print(figure(1),'part_2_a_sigma_and_V','-dpdf','-fillpage');

figure(2);
subplot(2,1,1);
pcolor(E_x);
xlabel('x(cm)');
ylabel('y(cm)');    
h = colorbar;
set(get(h,'title'),'string','E_x(V/m)');
title('Assignment 2 Part 2 A - E_x');
colormap jet;

subplot(2,1,2);
pcolor(E_y);
xlabel('x(cm)');
ylabel('y(cm)');    
h = colorbar;
set(get(h,'title'),'string','E_y(V/m)');
title('Assignment 2 Part 2 A - E_y');
colormap jet;
print(figure(2),'part_2_a_E_y_E_x','-dpdf','-fillpage');

figure(3);
subplot(2,1,1);
quiver(Y,X, E_x,E_y);
xlabel('x(cm)');
ylabel('y(cm)');    
title('Assignment 2 Part 2 A - Electric Field Vector');

subplot(2,1,2);
quiver(Y,X,J_x,J_y);
xlabel('x(cm)');
ylabel('y(cm)');    
title('Assignment 2 Part 2 A - Current Density');

print(figure(3),'part_2_a_E_J_quiver','-dpdf','-fillpage');
