function plotTraj(x,launchDateRange,sequence)
    mu = 1.32712440018E11;
    
    K_launchDate = x(1);
    K_TOFleg1 = x(2);
    K_TOFleg2 = x(3);
    K_TOFleg3 = x(4);

    launchDate = launchDateRange(1) +...
        K_launchDate*(launchDateRange(2) - launchDateRange(1)); %days
    TOFleg1 = K_TOFleg1*get_HohmanTOF(sequence(1),sequence(2)); %days
    TOFleg2 = K_TOFleg2*get_HohmanTOF(sequence(2),sequence(3)); %days
    TOFleg3 = K_TOFleg3*get_HohmanTOF(sequence(3),sequence(4)); %days

    %% Calculating Trajectory
    [r1, ~] = extractEphem(launchDate,sequence(1),true);
    [r2, vp2] = extractEphem(launchDate+TOFleg1,sequence(2),true);
    [r3, vp3] = extractEphem(launchDate+TOFleg1+TOFleg2,sequence(3),true);
    [r4, ~] = extractEphem(launchDate+TOFleg1+TOFleg2+TOFleg3,sequence(4),true);
    
    [~,v2i,~,~] = lambert(r1,r2,TOFleg1,0,mu);
    [v2f,v3i,~,~] = lambert(r2,r3,TOFleg2,0,mu);
    [v3f,~,~,~] = lambert(r3,r4,TOFleg3,0,mu);

    %% Calculating GA information
    dv_reqGA1 = v2f - v2i;
    v_relGA1 = v2i - vp2;
    alphaGA1 = pi - angleBetween(v_relGA1,dv_reqGA1);
    deltaGA1 = pi - 2*alphaGA1;
    
    dv_reqGA2 = v3f - v3i;
    v_relGA2 = v3i - vp3;
    alphaGA2 = pi - angleBetween(v_relGA2,dv_reqGA2);
    deltaGA2 = pi - 2*alphaGA2;

    % Plotting
    figure;
    title("Interplanetary Trajectory")
    hold on; grid on; axis equal
    plotPlanetaryOrbits(launchDate,sequence)
    
    plotPlanetV(launchDate,sequence(1))
    plotPlanetV(launchDate+TOFleg1,sequence(2))
    plotPlanetV(launchDate+TOFleg1+TOFleg2,sequence(3))
    plotPlanetV(launchDate+TOFleg1+TOFleg2+TOFleg3,sequence(4))

    plotLambertArc(r1,r2,TOFleg1,mu)
    plotLambertArc(r2,r3,TOFleg2,mu)
    plotLambertArc(r3,r4,TOFleg3,mu)

    figure;
    title("GA1 Passby Trajectory")
    hold on; grid on; axis equal
    plotGA(v_relGA1,deltaGA1,sequence(2))

    figure;
    title("GA2 Passby Trajectory")
    hold on; grid on; axis equal
    plotGA(v_relGA2,deltaGA2,sequence(3))
end

