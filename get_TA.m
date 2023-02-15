function TA = get_TA(r1,r2)
    % finds transfer angle between two 3D position vectors
    
    cosTA = dot(r1,r2)/(norm(r1)*norm(r2));
    sinTA = cross(r1,r2)/(norm(r1)*norm(r2));
    TA = rad2deg(atan2(norm(sinTA),cosTA));

    if sign(sinTA(3)) < 0
        % Transfer angle is greater than 180
        TA = 360 - TA;
    end
end

