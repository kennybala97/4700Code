clear all;
close all;
clc;

R1 = 1;
R2 = 2;
R3 = 97;
R4 = 0.1;
RO = 1000;
C = 0.25;
L = 10.2;
alpha = 100;
omega = linspace(0,100,1000);
V_in = linspace(-10,10,1000);

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
 
 for n = 1:1000
    F = [V_in(n);0;0;0;0;0];
    V1(n,:) = G\F;    
 end
 
figure;
plot(V_in, V1(:,5));
title('DC Sweep of Circuit', 'interpreter', 'latex');
ylabel('Output Voltage $V_O$', 'interpreter', 'latex');
xlabel('Input Voltage $V_{in}$', 'interpreter', 'latex');
 
  for n = 1:1000
    F = [10;0;0;0;0;0];
    V2(n,:) = (G + 1j*omega(n)*C)\F;    
  end
  
figure;
subplot(121);
plot(omega, real(V2(:,5)));

title('AC Sweep of Circuit', 'interpreter', 'latex');
ylabel('Output Voltage $V_O$', 'interpreter', 'latex');
xlabel('Frequency $\omega$', 'interpreter', 'latex');

subplot(122);
loglog(omega, 10*log10(real(V2(:,5))./10));

title('Circuit Gain', 'interpreter', 'latex');
ylabel('Gain $\frac{V_O}{V_i}$', 'interpreter', 'latex');
xlabel('Frequency $\omega$', 'interpreter', 'latex');
 
 
  for n = 1:1000
    C_rand = 0.05*randn(6,6).*C;
    F = [10;0;0;0;0;0];
    V3(n,:) = (G + 1j*pi*C_rand)\F;    
  end
  
figure;

histogram(real(V3(:,5))./10);

title('Circuit Gain with Random Perturbation', 'interpreter', 'latex');
xlabel('Gain $\frac{V_O}{V_i}$', 'interpreter', 'latex');
 
 