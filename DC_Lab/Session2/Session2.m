%% Exp 2
% Teacher : MS_Jafari
% Author: [SeyedAli] - [SeyedHosseini]
% E-mail: [alishosseini79@aut.ac.ir] 
%Student-Number : [9723042]
% University: Amirkabir University of Technology
%%
clc;
close all;
clear;
%% Initialization
Rs = 1e4;

SPS = 8;
alpha = 0.2;
Span = 10;
%% part 1
M = 16;
N = 1e3; 
D1 = randi([0,M-1 ],1,N);
clc;
% D = D - 1 ;
D = pskmod(D1,M);

H = [0,1,1,1,0,1,0,0,1,1,1,0,1,1,0,0,1,1,0,0,0,...
    0,1,1,1,0,0,1,1,1,0,1,1,0,0,0,1,1,0,1,0,1,0,0,1....
    ,0,0,0,0,1,1,0,0,1,1,0,0,0,0,0,1,0,1,0,0,1,1,1,....
    0,0,0,0,1,1,0,1,1,0,1,1,1,1,1,1,1,1,0,0,1,0,1....
    ,0,0,0,0,0,0,0,1,0,0,0,1,1,0,1,1,1,0,1,1,1,0,1....
    ,0,1,0,1,0,0,1,1,1,0,1,0,1,0];

hh = numel(H);
%% Plotting
clc;
figure(1)
plot(abs(xcorr((-1).^H)),"r")
grid on;
ylabel("Mag")
xlabel("time")
title("AutoCorr Header")
axis([0 260 -5 140])

%% Pilot
Ds = D/(std(D));
clc;
P = (1 + i)/(sqrt(2));
pp=20;
pilot = repmat(P,1,pp);
header = (-1).^H;

%% Frame
clc;
f1 = [pilot Ds];
f11 = repmat(f1,1,2);
f22 = [header Ds];
frame = [f22 f11]; %frame 
L = length(frame);
%% Scatter Plot
clc;
% figure(2)
scatterplot(frame)
grid on;
ylabel("imag")
xlabel("real")
title("Constl")
% axis([0 260 -5 140])

