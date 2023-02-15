clear, clc, close all;

options = optimoptions('fmincon', 'Algorithm', 'sqp',...
    'MaxFunctionEvaluations',1e4,'Display','iter');

launchDate = juliandate(2028,1,1) + yr2day(5.3);
sequence = [3,4,5,7];
uranus_vinf = 6;

TOFleg1 = 0.7*get_HohmanTOF(sequence(1),sequence(2));
TOFleg2 = 0.6*get_HohmanTOF(sequence(2),sequence(3));
TOFleg3 = 0.3*get_HohmanTOF(sequence(3),sequence(4));

IC = [TOFleg1,TOFleg2,TOFleg3,3.1e3,1.1e4];
lb = [0.25*get_HohmanTOF(sequence(1),sequence(2)),...
    0.6*get_HohmanTOF(sequence(2),sequence(3)),...
    0.25*get_HohmanTOF(sequence(3),sequence(4)),3e3,1e4];
ub = [get_HohmanTOF(sequence(1),sequence(2)),...
    get_HohmanTOF(sequence(2),sequence(3)),...
    get_HohmanTOF(sequence(3),sequence(4)),1e5,1e7];

nump = 100;

[x,fval] = fmincon(@(input)twoGA_Trajectory(input,launchDate,sequence,2),IC,[],[],...
        [],[],lb,ub,...
        @(input)nonLCons(input,launchDate,sequence,uranus_vinf),options);

TOFleg1 = x(1);
TOFleg2 = x(2);
TOFleg3 = x(3);

figure;
plotPlanetaryOrbits(launchDate,sequence)
[r1,vp1] = extractEphem(launchDate,sequence(1),true);
r2 = extractEphem(launchDate+TOFleg1,sequence(2),false);
r3 = extractEphem(launchDate+TOFleg1+TOFleg2,sequence(3),false);
r4 = extractEphem(launchDate+TOFleg1+TOFleg2+TOFleg3,sequence(4),false);

mu = 1.32712440018E11;
[v1,v2] = lambert(r1,r2,TOFleg1,0,mu);

plotLambertArc(r1,r2,TOFleg1,mu)
plotLambertArc(r2,r3,TOFleg2,mu)
plotLambertArc(r3,r4,TOFleg3,mu)