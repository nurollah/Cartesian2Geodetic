% In the name of Allah
function [p]=Geo2Cart(G,a,e)
% This code is used to convert geodetic coordinates to cartesian
% coordinates
% aboute inputs:
% G: is a vector that contains [phi, Lamda, Height]
% a: is the semi-major axis of the biaxial ellipsoid
% e: is eccentricity of the biaxial ellipsoid
format long;
phi=G(:,1);
la=G(:,2);
h=G(:,3);
N=a*ones(size(h,1),1)./sqrt(1-e^2*((sind(phi)).^2));% Eq (2).
p(:,1)=(N+h).*cosd(phi).*cosd(la);% X coordinate. Eq (1)
p(:,2)=(N+h).*cosd(phi).*sind(la);% Y coordinate. Eq (1)
p(:,3)=(N*(1-e^2)+h).*sind(phi);% Z coordinate. Eq (1)

