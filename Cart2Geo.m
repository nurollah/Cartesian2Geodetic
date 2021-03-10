% In th name of Allah
function [G]=Cart2Geo(p,a,e)
% This code is used to convert  cartesian coordinates to geodetic coordinates 
% If you have used this function, please also cite the following paper:
% Tatar, N., Farzaneh, S. (2018). Transforming Geocentric Cartesian 
% Coordinates to Geodetic Coordinates by a New Initial Value Calculation Paradigm.
% Journal of the Earth and Space Physics, 44(4), 19-28. 
% doi: 10.22059/jesphys.2018.246251.1006946
% data :  2018
% School of Surveying and Geomatics Engineering, College of Engineering,
% University of Tehran, Iran 
% this code implemented by Nurollah Tatar; Email: n.tatar@ut.ac.ir
% ---------------------------Note ----------------------------------------
% This code is allowed to use only for research purpose, and we
% don't provide any warranty. 
% p: is a vector that contains [X, Y, Z]
% a: is the semi-major axis of the biaxial ellipsoid
% e: is eccentricity of the biaxial ellipsoid
format long

b = a*sqrt(1-e^2);
e2 = sqrt(1-e^2);
a1 = 1/a;
b1 = 1/b;
a2 = a^2;
b2 = b^2;

x = p(:,1);
y = p(:,2);
z = p(:,3);
%
%
PG = sqrt(x.^2+y.^2); % Eq(5).

lamda = 2*atand(y./(x+PG)); % Eq(4).
lamda=360*(lamda<0)+lamda;

k = sqrt((PG*a1).^2+(z*b1).^2);   % Eq(9).
dd = (k-1).*(PG.^2+z.^2);
t0 = e2*(k.^2*a2+dd).*z./((k.^2*b2+dd).*PG+0.000001); % Eq(18).
C = ones(size(z,1),1)./sqrt(e2^2+t0.^2);  % Eq(20).
h = (e2*PG+z.*t0-b*sqrt(1+t0.^2)).*C; % Eq(19).
N = sqrt((a2-e^2*(PG-e2*h.*C).^2))/e2;    % Eq(25).
phi = atand((N+h).*z./((N*e2^2+h).*PG+0.000001)); % Eq(26).

G(:,1:3) = [phi,lamda,h];
