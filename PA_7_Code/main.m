clear all;
close all;
clc;

R1 = 1;
R2 = 2;
R3 = 10;
R4 = 0.1;
RO = 1000;
C = 0.25;
L = 0.2;
alpha = 100;
omega = linspace(0,100,100);
V_in = linspace(-10,10,100);

C = [0, 0, 0, 0, 0, 0;...
     C,  -C, 0, 0, 0, 0;...
     0,  0, 0, 0, 0, 0;...
     0,  0, 0, 0, 0, 0;...
     0,  0, 0, 0, 0, -L;...
     0,  0, 0, 0, 0, 0];

G = [1, 0, 0, 0, 0, 0;...
     1/R1, -1/R1 - 1/R2, 0, 0, 0, -1;...
      0,  0, -1/R3, 0, 0, 1;...
      0,  0, 0, 1/R4, -(1/R4 + 1/RO), 0;...
      0,  1, -1, 0, 0, 0;...
      0,  0, 0, 1, 0, -alpha];
 
 for n = 1:100
    F = [V_in(n);0;0;0;0;0];
    V1(n,:) = G\F;    
 end
 
 figure;
 plot(V_in)
 hold on;
 plot(V1(:,5));
 
  for n = 1:100
    F = [10;0;0;0;0;0];
    V2(n,:) = (G + 1j*omega(n)*C)\F;    
 end
 
 figure;
 plot(omega,abs(V2(:,5)));
 
   for n = 1:100
    F = [10;0;0;0;0;0];
    V2(n,:) = (G + 1j*omega(n)*C)\F;    
 end
 
 figure;
 plot(omega,abs(V2(:,5)));