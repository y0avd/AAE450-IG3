function delV = calcDelV(launch,TOF)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

mu = 1.32712440018E11;

TOF = TOF*365.25;

[r1,ve] = extractEphem(launch,3,true);

[r2,vu] = extractEphem(launch+TOF,7,true);

[v1,v2] = lambert(r1,r2,TOF,0,mu);

vinf = norm(v1-ve);
a = -(398600.4)/vinf^2;
vp = sqrt(398600.4*((2/6578)+(1/-a)));
dVe = vp - sqrt(398600.4*((2/6578)-(1/21182)));

vinf = norm(v2-vu);
a = -(5.793939e6)/vinf^2;
vp = sqrt(5.793939e6*((2/33500)+(1/-a)));
dVu = vp - sqrt(5.793939e6*((2/33500)-(1/7.8468e+05)));

delV = dVe+dVu;

end