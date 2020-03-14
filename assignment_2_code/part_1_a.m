clear all;
close all;
clc;

nx = 100;
ny = 100;

delta = 1e-2;

V = 1;

initial_mesh = nan(nx,ny);

initial_mesh(1,:) = V;
initial_mesh(end,:) = 0;

[V] = laplace_solver_1d(nx,ny,delta,initial_mesh);

figure(1);
surf(V);
title('Assignment 2 Part 1A - Solution of Laplace''s Equation for the case of \partial{V}/\partial{y} = 0');
xlabel('x(cm)');
ylabel('y(cm)');
zlabel('Voltage(V)');
view(+45,-45);
colormap jet;
print(figure(1),'part_1_a_assignment','-dpdf','-fillpage')
