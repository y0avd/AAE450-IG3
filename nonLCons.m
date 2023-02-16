% K_TOFleg1 = input(1); %0 to 1
% K_TOFleg2 = input(2); %0 to 1
% K_TOFleg3 = input(3); %0 to 1

function [c, ceq] = nonLCons(input,launchDateRange, sequence, uranus_vinf)
% This function contains all non linear constraints, namely matching dV
% for gravity assists and vinf at uranus

%% Initializations
mu = 1.32712440018E11;

%dVmax = 40;
%TOFmax = 18*365.2422;

K_launchDate = input(1);
K_TOFleg1 = input(2);
K_TOFleg2 = input(3);
K_TOFleg3 = input(4);

launchDate = launchDateRange(1) +...
    K_launchDate*(launchDateRange(2) - launchDateRange(1)); %days
TOFleg1 = K_TOFleg1*get_HohmanTOF(sequence(1),sequence(2)); %days
TOFleg2 = K_TOFleg2*get_HohmanTOF(sequence(2),sequence(3)); %days
TOFleg3 = K_TOFleg3*get_HohmanTOF(sequence(3),sequence(4)); %days


%% Calculating Trajectory
[r1, vp1] = extractEphem(launchDate,sequence(1),true);
[r2, vp2] = extractEphem(launchDate+TOFleg1,sequence(2),true);
[r3, vp3] = extractEphem(launchDate+TOFleg1+TOFleg2,sequence(3),true);
[r4, vp4] = extractEphem(launchDate+TOFleg1+TOFleg2+TOFleg3,sequence(4),true);

[v1,v2i,~,~] = lambert(r1,r2,TOFleg1,0,mu);
[v2f,v3i,~,~] = lambert(r2,r3,TOFleg2,0,mu);
[v3f,v4,~,~] = lambert(r3,r4,TOFleg3,0,mu);

%% Matching Uranus Vinf
%ceq(1) = (norm(v4 - vp4) - uranus_vinf)/uranus_vinf;

%% Matching dV for GAs
dv_reqGA1 = v2f - v2i;
v_relGA1 = v2i - vp2;
alphaGA1 = pi - angleBetween(v_relGA1,dv_reqGA1);
deltaGA1 = pi - 2*alphaGA1;

[dvGA1, passbyr1] = get_dvGA(v_relGA1,deltaGA1,sequence(2));

dv_reqGA2 = v3f - v3i;
v_relGA2 = v3i - vp3;
alphaGA2 = pi - angleBetween(v_relGA2,dv_reqGA2);
deltaGA2 = pi - 2*alphaGA2;

[dvGA2, passbyr2] = get_dvGA(v_relGA2,deltaGA2,sequence(3));

% matching dV magnitude
K = 1/10; % constraint scaling
ceq(1) = K*(norm(dv_reqGA1) - dvGA1);
ceq(2) = K*(norm(dv_reqGA2) - dvGA2);

%% Max dV and TOF
%totaldV = norm(v1 - vp1);
%totalTOF = TOFleg1 + TOFleg2 + TOFleg3;

%c(1) = (totaldV - dVmax)/dVmax;
%c(2) = (totalTOF - TOFmax)/TOFmax;

%% Verifying s/c does not intersect planet
c(1) = (getPlanetRpLim(sequence(2)) - passbyr1)/getPlanetRpLim(sequence(2));
c(2) = (getPlanetRpLim(sequence(3)) - passbyr2)/getPlanetRpLim(sequence(3));

%% Verifying delta of gravity assist is reasonable
K = -1/pi; %scaling
c(3) = K*deltaGA1;
c(4) = K*deltaGA2;

end

