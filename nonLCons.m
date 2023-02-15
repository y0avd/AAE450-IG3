% TOFleg1 = input(1); %days
% TOFleg2 = input(2); %days
% TOFleg3 = input(3); %days
% flybyr1 = input(4); %km
% flybyr2 = input(5); %km

function [c, ceq] = nonLCons(input,launchDate, sequence, uranus_vinf)
% This function contains all non linear constraints, namely matching dV
% for gravity assists and vinf at uranus

%% Initializations
mu = 1.32712440018E11;

dVmax = 15;
TOFmax = 18*365.2422;

TOFleg1 = input(1); %days
TOFleg2 = input(2); %days
TOFleg3 = input(3); %days
flybyr1 = input(4); %km
flybyr2 = input(5); %km

%% Calculating Trajectory
[r1, vp1] = extractEphem(launchDate,sequence(1),true);
[r2, vp2] = extractEphem(launchDate+TOFleg1,sequence(2),true);
[r3, vp3] = extractEphem(launchDate+TOFleg1+TOFleg2,sequence(3),true);
[r4, vp4] = extractEphem(launchDate+TOFleg1+TOFleg2+TOFleg3,sequence(4),true);

[v1,v2i,~,~] = lambert(r1,r2,TOFleg1,0,mu);
[v2f,v3i,~,~] = lambert(r2,r3,TOFleg2,0,mu);
[v3f,v4,~,~] = lambert(r3,r4,TOFleg3,0,mu);

%% Matching Uranus Vinf
ceq(1) = norm(v4 - vp4') - uranus_vinf;

%% Matching dV for GAs
dvGA1_req = v2f - v2i;
deltaGA1_req = abs(acos(dot(v2i,dvGA1_req)/(norm(v2i)*norm(dvGA1_req))));
dvGA2_req = v3f - v3i;
deltaGA2_req = abs(acos(dot(v3i,dvGA2_req)/(norm(v3i)*norm(dvGA2_req))));

[dvGA1, deltaGA1] = get_dvGA(vp2,v2i,flybyr1,sequence(2));
[dvGA2, deltaGA2] = get_dvGA(vp3,v3i,flybyr2,sequence(3));

% matching dV magnitude
ceq(2) = norm(dvGA1_req) - dvGA1;
ceq(3) = norm(dvGA2_req) - dvGA2;

% matching dV angle
ceq(4) = deltaGA1_req - deltaGA1;
ceq(5) = deltaGA2_req - deltaGA2;

%% Max dV and TOF
totaldV = norm(v1 - vp1);
totalTOF = TOFleg1 + TOFleg2 + TOFleg3;

c(1) = totaldV - dVmax;
c(2) = totalTOF - TOFmax;

end

