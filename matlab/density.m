function rho = density(T)

% Calculate the density of water for a given temperature from Chen and
% Millero for pure water
% DJW 2/15/08

rho = 0.9998395+6.7914*10^(-5)*T-9.0894*10^(-6)*T.^2+1.0171*10^(-7)*T.^3-1.2846*10^(-9)*T.^4+1.1592*10^(-11)*T.^5-5.0125*10^(-14)*T.^6;
rho = rho*1000;