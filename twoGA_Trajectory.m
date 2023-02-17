% K_launchDate = input(1) %0 to 1
% K_TOFleg1 = input(2); %0 to 1
% K_TOFleg2 = input(3); %0 to 1
% K_TOFleg3 = input(4); %0 to 1

% min_var = 1 (TOF) 2 (dV)

function output = twoGA_Trajectory(input,launchDateRange,sequence,min_var)
%% Initializations
mu = 1.32712440018E11;

K_launchDate = input(1);
K_TOFleg1 = input(2);
K_TOFleg2 = input(3);
K_TOFleg3 = input(4);

launchDate = launchDateRange(1) +...
    K_launchDate*(launchDateRange(2) - launchDateRange(1)); %days
TOFleg1 = K_TOFleg1*get_HohmanTOF(sequence(1),sequence(2)); %days
TOFleg2 = K_TOFleg2*get_HohmanTOF(sequence(2),sequence(3)); %days
TOFleg3 = K_TOFleg3*get_HohmanTOF(sequence(3),sequence(4)); %days

%% Calulate Trajectory
[r1,vp1] = extractEphem(launchDate,sequence(1),true);
r2 = extractEphem(launchDate+TOFleg1,sequence(2),false);
r3 = extractEphem(launchDate+TOFleg1+TOFleg2,sequence(3),false);
[r4,vp4] = extractEphem(launchDate+TOFleg1+TOFleg2+TOFleg3,sequence(4),true);
[v1,~,~,~] = lambert(r1,r2,TOFleg1,0,mu);
[~,v4,~,~] = lambert(r3,r4,TOFleg3,0,mu);

if min_var == 1
    output = TOFleg1 + TOFleg2 + TOFleg3;
elseif min_var == 2
    % Calculate deltaV from LEO
    earth_vinf = norm(v1 - vp1);
    output = earth_vinf;
%     uranus_vinf = norm(v4 - vp4);
%     output = earth_vinf + uranus_vinf;
    % DOES NOT INCLUDE dV at Uranus (either with or without aerocapture)
    % POTENTIALLY ADD OBERTH dV
else
    disp('ERROR: Invalid min_var')
end
