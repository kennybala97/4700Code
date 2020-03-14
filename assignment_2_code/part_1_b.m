clear all;
close all;
clc;

% filename = 'part_1_b_analytic_vs_numerical_x_150_x_100.gif';

nx = 150;
ny = 100;
delta = 1e-2;

L = nx*delta;
W = ny*delta;

V0 = 1;

initial_mesh = nan(nx,ny);

initial_mesh(1,:) = 0;
initial_mesh(end,:) = 0;

initial_mesh(1,:) = V0;
initial_mesh(end,:) = V0;

[V_numeric] = laplace_solver_2d(nx,ny,delta,initial_mesh);
% surf(V_numeric);
% colormap jet;
% xlabel('x(cm)');
% ylabel('y(cm)');
% title('Assignment 2 Part 1 B - Numerical Solution');
% h = colorbar;
% set(get(h,'title'),'string','Voltage(V)');
% print(figure(1),'part_1_b_assignment_numerical','-depsc')


[X,Y] = meshgrid(1:1:nx,1:1:ny);
X = (X')*delta;
Y = (Y')*delta;

analytic_soln = zeros(nx,ny);

a = W;
b = L/2;
x = X - L/2;
y = Y;

f = figure('units','normalized','outerposition',[0 0 1 1]);

colormap(f,jet);

for n = 1:2:200
    analytic_soln = analytic_soln + (1/n).*(cosh((n.*pi.*x)./a))./(cosh((n.*pi.*b)./a)).*sin((n.*pi.*y)./a);
    
%     subplot(3,1,1);
%     pcolor(analytic_soln');
%     xlabel('x(cm)');
%     ylabel('y(cm)');
%     title('Assignment 2 Part 1 B - Analytic Solution');
%     colorbar;
% 
%     subplot(3,1,2);    
%     pcolor(V_numeric');
%     xlabel('x(cm)');
%     ylabel('y(cm)');
%     colorbar;
%     title('Assignment 2 Part 1 B - Numerical Solution');
%     
%     subplot(3,1,3);
%     pcolor(analytic_soln' - V_numeric');
%     xlabel('x(cm)');
%     ylabel('y(cm)');    
%     colorbar;
%     title('Assignment 2 Part 1 B - Analytic vs. Numerical Solution Difference');
% 
%     sgtitle(sprintf('L = %i cm  W = %i cm n = %i',nx,ny,n));
%     
%     pause(0.1);
%     
%     drawnow
%     
%   % Capture the plot as an image 
% 
%   frame = getframe(f); 
% 
%   im = frame2im(frame); 
%   
%   [imind,cm] = rgb2ind(im,256); 
% 
%   % Write to the GIF File 
% 
%   if n == 1 
%       imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
%   else 
%       imwrite(imind,cm,filename,'gif','WriteMode','append'); 
%   end
end

subplot(3,1,1);
pcolor(analytic_soln');
xlabel('x(cm)');
ylabel('y(cm)');
title('Assignment 2 Part 1 B - Analytic Solution');
colorbar;

subplot(3,1,2);    
pcolor(V_numeric');
xlabel('x(cm)');
ylabel('y(cm)');
colorbar;
title('Assignment 2 Part 1 B - Numerical Solution');

subplot(3,1,3);
pcolor(analytic_soln' - V_numeric');
xlabel('x(cm)');
ylabel('y(cm)');    
colorbar;
title('Assignment 2 Part 1 B - Analytic vs. Numerical Solution Difference');

sgtitle(sprintf('L = %i cm  W = %i cm n = %i',nx,ny,n));
print(f,'part_1_b_assignment_analytical_numerical','-dpdf','-fillpage')

