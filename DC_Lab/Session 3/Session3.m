%% Exp 3
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
D2 = qammod(D1,M);
stdD = std(D2);
D = D2/stdD;

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

 %% Coefficients
% frameSize = size(frame);
%  numFrames = 1; nSamples = numFrames*frameSize;
%  DampingFactor = 1.21; NormalizedLoopBandwidth = 1e-4;
%      %% Calculate range estimates
%      NormalizedPullInRange = min(1, 2*pi*sqrt(2)*DampingFactor*...
%      NormalizedLoopBandwidth);
%      MaxFrequencyLockDelay = (4*NormalizedPullInRange^2)/...
%      (NormalizedLoopBandwidth)^3;
%      MaxPhaseLockDelay = 1.3/(NormalizedLoopBandwidth);
%      %% Impairments
%      frequencyOffsetHz = Rs*(NormalizedPullInRange);
%      PhaseRecoveryLoopBandwidth = NormalizedLoopBandwidth*SPS;
%      PhaseRecoveryGain = SPS;
%      PhaseErrorDetectorGain = log2(M); DigitalSynthesizerGain = -1;
%      theta = PhaseRecoveryLoopBandwidth/...
%      ((DampingFactor + 0.25/DampingFactor)*SPS);
%     delta = 1 + 2*DampingFactor*theta + theta*theta;
%      % G1
%      KP = (4*DampingFactor*theta/delta)/...
%      (PhaseErrorDetectorGain*PhaseRecoveryGain);
%      % G3
%      KI = (4/SPS*theta*theta/delta)/...
%      (PhaseErrorDetectorGain*PhaseRecoveryGain);