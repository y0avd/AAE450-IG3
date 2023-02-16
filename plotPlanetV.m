function plotPlanetV(t,planetID)
    % velocity vector scaling
    K = 1e7;

    [r,vp] = extractEphem(t,planetID,true);

    plot3([r(1),r(1) + K*vp(1)],...
        [r(2),r(2) + K*vp(2)],...
        [r(3),r(3) + K*vp(3)],'b')
end

