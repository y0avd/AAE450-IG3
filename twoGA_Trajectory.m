% TOFleg1 = input(1); %days
% TOFleg2 = input(2); %days
% TOFleg3 = input(3); %days
% flybyr1 = input(4); %km
% flybyr2 = input(5); %km

% min_var = 1 (TOF) 2 (dV)

function output = twoGA_Trajectory(input,launchDate,sequence,min_var)
%% Initializations
mu = 1.32712440018E11;

TOFleg1 = input(1); %days
TOFleg2 = input(2); %days
TOFleg3 = input(3); %days

%% Calulate Trajectory
[r1,vp1] = extractEphem(launchDate,sequence(1),true);
r2 = extractEphem(launchDate+TOFleg1,sequence(2),false);
[v1,~,~,~] = lambert(r1,r2,TOFleg1,0,mu);

if min_var == 1
    output = TOFleg1 + TOFleg2 + TOFleg3;
elseif min_var == 2
    % Calculate deltaV from LEO
    output = norm(v1 - vp1);

    % DOES NOT INCLUDE dV at Uranus (either with or without aerocapture)

    % POTENTIALLY ADD OBERTH dV
else
    disp('ERROR: Invalid min_var')
end