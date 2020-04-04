clear all;
close all;
clc;

[bottlenecks,I_av] = part_2_with_bottleneck();

eqn = @(V)