% FIND_NPROF_T1_T2_MODES.M
% Uses the temperature profile (z,T) and computes profile of buoyancy
% frequency (zN,N) to solve for the V1 and V2 seiche modal structure. The
% thermocline depth (zt) is defined as the depth of maximum displacement
% during a V1 seiche and is the peak in the V1 modal structure. The
% boundaries of the metalimnion (zm1, zm2) are the depths of maximum
% displacement during a V2 seiche and are identified from the peaks in the
% V2 modal structure. From this and the basin length L, periods of the V1H1
% (T1)and V2H1 (T2) seiches are computed. This is based on Munnich et al
% 1992 and solves for the modes using many layers as opposed to a two or
% three layer approximation.
%
% Depth input is positive downward as that is how most people's data are.
% 
% DJW, JAS, MSK 6/10/09
% DJW has updated many times since and not logged the changes
% DJW 11/9/21 Cleaning up for GLEON project

function [N,zN,T1,T2,zt,zm1,zm2] = find_V1V2_modes(z,T,L)

% GLEON: Uses DJW function buoyfreq, could be replaced by either calling
% function from Lake Analyser in this function or changing the input to a
% buoyancy frequency profile instead of a temperature profile. The profile
% output does need to be on a regular spaced vertical grid.
[N,zN] = buoyfreq(z,T);

% For solving modes, there needs to be a full profile of N from surface to
% bottom. This extends the profile to the top and the bottom.
indx = find(~isnan(N));
N(1:indx(1)) = N(indx(1));
N(indx(end):end) = N(indx(end));

% Normalize N^2/mean(N)^2
Nmean = mean(N);
gamma2 = (N/Nmean).^2;

% Normalize z by max depth
H = max(zN);
D = zN/H;

% Solve the eigenvalue problem
nseg = length(N)+1;

A = zeros*ones(nseg-1,nseg-1);
A(1,1) = -2;
A(1,2) = 1;
A(nseg-1,nseg-2) = 1;
A(nseg-1,nseg-1) = -2;
for i = 2:nseg-2
    A(i,i-1) = 1;
    A(i,i) = -2;
    A(i,i+1) = 1;
end

delta = 1/nseg;

B = zeros*ones(nseg-1,nseg-1);
for i = 1:nseg-1
    B(i,i) = gamma2(i);
end
B = -delta^2*B;

[phi,lambda] = eig(A,B); % phi is the amplitude

% Find two smallest eigenvalues, which correspond to V1 and V2
dl = diag(lambda);
indxinf = dl == -Inf;
dl(indxinf) = nan;

lmin = min(dl);
index(1) = find(dl == lmin);
lambda1 = lmin; %1st vertical mode eigenvalue
dl(index(1)) = nan; %replace smallest with nan

lmin = min(dl);
index(2) = find(dl == lmin);
lambda2 = lmin; %2nd vertical mode eigenvalue

% wave speeds from eigenvalues
c1 = sqrt(1/lambda1)*Nmean*H;
c2 = sqrt(1/lambda2)*Nmean*H;

% period of the v1h1 mode
T1 = 2*L/c1;
% period of the v2h1 mode
T2 = 2*L/c2;

% phi contains information on the modal structure. Finding the minimum two
% eigenvalues above, we have found which columns contain the modal
% structure for those two modes. The peak in the V1 profile corresponds
% with the thermocline and the peaks in the V2 profile correspond to the
% boundaries of the metalimnion from an internal seiche energetics
% standpoint as there are the depths of the the maximum seiche amplitude. 
phi1 = abs(phi(:,index(1)));
indx_zT = phi1==max(phi1);
zt = zN(indx_zT); 

phi2 = (phi(:,index(2)));
indx_zm1 = phi2==max(phi2);
temp1 = zN(indx_zm1);
indx_zm2 = phi2==min(phi2);
temp2 = zN(indx_zm2);
zm1 = min([temp1 temp2]);
zm2 = max([temp1 temp2]);
