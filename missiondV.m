function [dVtot,dVe,dVu] = missiondV(vInfe,vInfu)

%% CONSTANTS
muE = 3.986004418e5;
muU = 5.793939e6;

%% ORBIT PARAMETERS
rpE = 6578;
raE = 42164;

rpU = 33500;
raU = 1535860;

%% CALCULATIONS
aEe = 0.5*(rpE+raE);
aHe = -(muE)/vInfe^2;
vpHe = sqrt(muE*((2/rpE)-(1/aHe)));
vpEe = sqrt(muE*((2/rpE)-(1/aEe)));
dVe = vpHe - vpEe;

aEu = 0.5*(rpU+raU);
aHu = -(muU)/vInfu^2;
vpHu = sqrt(muU*((2/rpU)-(1/aHu)));
vpEu = sqrt(muU*((2/rpU)-(1/aEu)));
dVu = vpHu - vpEu;

dVtot = dVe+dVu;
end