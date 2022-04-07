function [N,zN] = buoyfreq(z,T,dz)

% Calculate the buoyancy frequency given a temperature profile.
% input z and output zN are positive downwards
% dz is spacing for interpolation
% N is in rad/s
% DJW 2/12/08
% DJW 12/22/21 to make clearer

g = 9.81;

% Default dz = 0.1 m unless otherwise specified
if nargin < 3
    dz = 0.1;
end

% Interpolate onto regular grid of dz
zN = 0:dz:max(z);
Ti = interp1(z,T,zN);

% GLEON: Uses DJW function density based on Chen and Millero
rho = density(Ti);

% drhodz calculated as a center difference
drhodz = nan*ones(length(rho),1);
for i = 2:length(rho)-1
    drhodz(i) = (rho(i+1)-rho(i-1))/(zN(i+1)-zN(i-1));
end

% Because z is positive downwards, drhodz should always be positive
ind = drhodz<0;
drhodz(ind) = 0;

N = sqrt(g/mean(rho,'omitnan')*drhodz);


