function [dvGA, delta] = get_dvGA(vp,vi,flybyr,planetID)
    mu = getPlanetMu(planetID);

    v_rel_in = vi - vp;
    v_inf = norm(v_rel_in);
    specE = v_inf^2/2;
    a = -mu/(2*specE);
    e = 1 - flybyr/a;
    delta = 2*asin(1/e); %rads

    dvGA = sqrt(2*v_inf^2 - 2*v_inf^2*cos(delta));
end

