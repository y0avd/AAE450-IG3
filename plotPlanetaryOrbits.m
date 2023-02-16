function plotPlanetaryOrbits(launchDate, sequence)
    nump = 1e3;

    r_sun = 696000;
    [X,Y,Z] = sphere;
    surf(X*r_sun,Y*r_sun,Z*r_sun)

    for i = 1:length(sequence)
        T = getPlanetT(sequence(i))/(3600*24);

        tspan = linspace(0,T,nump);

        for k = 1:nump
            r = extractEphem(launchDate+tspan(k),sequence(i),false);
            x(k) = r(1); %#ok<AGROW> 
            y(k) = r(2); %#ok<AGROW> 
            z(k) = r(3); %#ok<AGROW> 
        end

        plot3(x,y,z,'k')
    end
end

