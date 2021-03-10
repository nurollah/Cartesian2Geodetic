% In the name of Allah
% -------------------------- Demo Code ------------------------------------
% If you have used this function, please also cite the following paper:
% paper : Tatar, N., Farzaneh, S. (2018). Transforming Geocentric Cartesian 
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
%==========================================================================
clc
clear 
format long;

a = 6378137;% the semi-major axis of WGS84
b = 6356752.314;% the semi-minor axis of WGS84
e = sqrt((a^2-b^2)/a^2);%is eccentricity of the biaxial ellipsoid

% load data
[filename,pathname,filterindex] = uigetfile({'*.txt','TEXT(*.txt)';'*.*','all File(*.*)'},'Open');
p1 = load([pathname filename]);

p = Geo2Cart(p1,a,e);
G = Cart2Geo(p,a,e);
% dphi, dlamda
d_coordinate = max(abs(G(:,1:2)-p1(:,1:2)))*pi/180
dh = max(abs(G(:,3)-p1(:,3)))