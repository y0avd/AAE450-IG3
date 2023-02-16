[r1,ve] = extractEphem(2463464.5,3,true);

[r2,vu] = extractEphem(2.467160751064000e+06,7,true);

[v1,v2] = lambert(r1,r2,10.12*365.2422,0,1.32712440018E11);

% [v1b,v2b] = LAMBERTBATTIN(r1'.*1000,r2'.*1000,'pro',10.12*365.2422*86400);
% v1b = v1b./1000;
% v2b = v2b./1000;

vinf = norm(v1-ve)
a = -(398600.4)/vinf^2;
vp = sqrt(398600.4*((2/6578)+(1/-a)))
dVe = vp - sqrt(398600.4*((2/6578)-(1/21182)));

vinf = norm(v2-vu);
a = -(5.793939e6)/vinf^2;
vp = sqrt(5.793939e6*((2/33500)+(1/-a)));
dVu = vp - sqrt(5.793939e6*((2/33500)-(1/7.8468e+05)));

dVtot = dVe+dVu

