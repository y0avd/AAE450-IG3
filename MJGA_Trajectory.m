function output = twoGA_Trajectory(input,launchDate,Sequence,uranus_vinf,min_var)
%% Initializations
mu = 1.32712440018E11;

TOFleg1 = input(1); %days
TOFleg2 = input(2); %days
TOFleg3 = input(3); %days
flybyr1 = input(4); %km
flybyr2 = input(5); %km

%% Calulate Trajectory
[r1, vp1] = extractEphem(launchDate,3);
[r2, vp2] = extractEphem(launchDate + TOFleg1,Sequence(1));
[r3, vp3] = extractEphem(launchDate + TOFleg1 + TOFleg2,Sequence(2));
[r4, vp4] = extractEphem(launchDate + TOFleg1 + TOFleg2 + TOFleg3,7);

[v1,v2i,~,~] = lambert(r1',r2',TOFleg1,0,mu);
[v2f,v3i,~,~] = lambert(r2',r3',TOFleg2,0,mu);
[v3f,v4,~,~] = lambert(r3',r4',TOFleg3,0,mu);

%% Calculate dV required
dvreq1 = v2f - v2i;
dvreq2 = v3f - v3i;

%% Calculate dV from gravity assist
dvG1 = get_dvGA(r2,vp2,v2i,flybyr1,Sequence(1));
dvG2 = get_dvGA(r3,vp3,v3i,flybyr2,Sequence(2));

if min_var == "TOF"
    output = totalTOF;
else
    % Calculate deltaV
    output = totaldV;
end
end