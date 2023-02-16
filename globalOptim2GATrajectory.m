clear, clc, close all;

options = optimoptions('ga','Display','iter');

launchDateRange = [juliandate(2031,2,1),juliandate(2031,3,1)];
sequence = [3,4,5,7];
uranus_vinf = 0;

lb = [0, 0.1, 0.1, 0.1];
ub = [1, 1, 1, 1];
nvars = 4;

[x,fval] = ga(@(input)twoGA_Trajectory(input,launchDateRange,sequence,2),nvars,[],[],...
        [],[],lb,ub,...
        @(input)nonLCons(input,launchDateRange,sequence,uranus_vinf),options);

plotTraj(x,launchDateRange,sequence)
printOrbitInfo(x,launchDateRange,sequence)