clear all;
close all;
clc;

components.R1 = 1;
components.R2 = 2;
components.R3 = 97;
components.R4 = 0.1;
components.RO = 1000;
components.C = 0.25;
components.L = 0.2;
components.alpha = 100;

n = 1000;
dt = 1/n;
t = 0:dt:(dt*(n-1));

V_in_step = zeros(1,n);
V_in_step(t > 0.03) = 1;

[V] = transient_simulator(components,V_in_step,dt,n);

figure;
subplot(121);
plot(t,V_in_step);
hold on;
plot(t,V(5,:));
title('Transient Response To A Step', 'interpreter', 'latex');
ylabel('Voltage', 'interpreter', 'latex');
xlabel('Time(s)', 'interpreter', 'latex');
legend('$V_{in}$','$V_O$', 'interpreter', 'latex');

input_spect = fftshift(fft(V_in_step));
output_spect = fftshift(fft(V(5,:)));

subplot(122);
semilogy(10*log10(abs(input_spect(500:end))));
hold on;
semilogy(10*log10(abs(output_spect(500:end))));

title('Welch PSD Estimate', 'interpreter', 'latex');
ylabel('Power', 'interpreter', 'latex');
xlabel('Frequency', 'interpreter', 'latex');
legend('$V_{in}$','$V_O$', 'interpreter', 'latex');


V_in_sine = sin(2*pi*t./0.03);

[V] = transient_simulator(components,V_in_sine,dt,n);

figure;
subplot(121);
plot(t,V_in_sine);
hold on;
plot(t,V(5,:));
title('Transient Response To A Step', 'interpreter', 'latex');
ylabel('Voltage', 'interpreter', 'latex');
xlabel('Time(s)', 'interpreter', 'latex');
legend('$V_{in}$','$V_O$', 'interpreter', 'latex');

subplot(122);
[pxx,f] = pwelch(V_in_sine,128,0,0:0.1:10,1/dt,'onesided');
semilogy(f,10*log10(pxx));
hold on;
[pxx,f] = pwelch(V(5,:),128,0,0:0.1:10,1/dt,'onesided');
semilogy(f,10*log10(pxx));
title('Welch PSD Estimate', 'interpreter', 'latex');
ylabel('Power', 'interpreter', 'latex');
xlabel('Frequency', 'interpreter', 'latex');
legend('$V_{in}$','$V_O$', 'interpreter', 'latex');




