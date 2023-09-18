function dfdt = rates(~,f)
% 
%...Constants;
mu = 398600; %Gravitational parameter (km^3/s^2)
RE = 6378; %Earth’s radius (km)
wE = [ 0 0 7.2921159e-5]'; %Earth’s angular velocity (rad/s)
%...Satellite data:
CD = 2.2; %Drag codfficient
m = 100; %Mass (kg)
A = pi/4*(1^2) ; %Frontal area (m^2)
%
% This function calculates the spacecraft acceleration from its
% position and velocity at time t.
% ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
x    = f(1);
y    = f(2);
z    = f(3);
vx   = f(4);
vy   = f(5);
vz   = f(6);
R = [x y z]'; %Position vector (km/s)
r = norm(R); %Distance from earth’s center (km)
alt = r - RE; %Altitude (km)
rho = atmosphere(alt); %Air density from US Standard Model (kf/m^3)
V = [vx vy vz]'; %Velocity vector (km/s)
Vrel = V - cross(wE,R); %Velocity relative to the atmosphere (km/s)
vrel = norm(Vrel); %Speed relative to the atmosphere (km/s)
uv = Vrel/vrel; %Relative velocity unit vector
ap = -CD*A/m*rho*... %Acceleration due to drag (m/s^2)
(1000*vrel)^2/2*uv; %(converting units of vrel from km/s to m/s)
a0 = -mu*R/r^3; %Gravitational ecceleration (km/s^2)
a = a0 + ap/1000; %Total acceleration (km/s^2)
dfdt = [V a]'; %Velocity and the acceleraion returned to ode45
end %rates