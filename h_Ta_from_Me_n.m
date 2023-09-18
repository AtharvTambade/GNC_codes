function [h, Ta] = h_Ta_from_Me_n(e,M,n)
deg = pi/180;
M=120.173*deg;
e= 0.0074922;
n=14.77061887;
u=398600;
E=kepler_E(e, M);
n0=n*2*pi/(24*60*60);

Ta=2*atan(((1+e)/(1-e))^(1/2)*tan(E/2))
h=(u^2/n0)^(1/3)*sqrt(1-e^2)
a=(h^2/u)*(1/(1-e^2))