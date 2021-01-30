clear all;
close all;
clc;

pkg load symbolic

% Const
c = 2.99792458E8;
%tag true x,y,z position:
display('Input tag coordinates');
xt = 0.32
yt = 4.65
zt = 1.23
% Anchors x,y,z positions:
W=[
0 0 1.35    % anchor 0
6.05 0.1 1.72   % anchor 1
4.3 3.40 1.6  % anchor 2
0.40 3.80 1.62     % anchor 3
];

% number of anchors
N = length(W);

% compute true distances between a tag and each anchor
for i=1:N
    di(i) = sqrt((xt-W(i,1))^2 + (yt-W(i,2))^2 + (zt-W(i,3))^2);

end

% compute true time differences of arrival due to anchor 0 with added noise (measurement error)
for i=1:N
    dti0(i) = (di(i)-di(1))/c + rand*1E-12;
end

display('time differences');
dti0
%%------------------------------------------ ETAP 2
% Now let's calculate the tag coordinates based on the given time differences of arrival

% Inputs:
    % anchors coordinates W
    % Time difference of arrival of the tag's blink broadcast due to anchor 0

% create equations to solve
% 
% di^2 = (x-xi)^2 + (y-yi)^2 + (z-zi)^2
% 
% di^2 - d0^2 = (x-xi)^2 + (y-yi)^2 + (z-zi)^2 - (x-x0)^2 - (y-y0)^2 - (z-z0)^2 = ... = (di-d0)(di+d0)
% 
% 
% 
% Coefficients matrix
% 
% A = [x0-x1 y0-y1 z0-z1
%  x0-x2 y0-y2 z0-z2
%  x0-x3 y0-y3 z0-z3]
% 
% [X] = [x
%        y
%        z]
%
%  [A][X] = [B] ... let's find [X]
% 
A = zeros(3,3);
for i=1:N-1
    A(i,1) = W(1,1) - W(i+1,1);
    A(i,2) = W(1,2) - W(i+1,2);
    A(i,3) = W(1,3) - W(i+1,3);
end

% Compute inverse of the matrix A
Ainverse = invertMatrix3rdDegree(A)

% Compute distance differences due to anchor 0
for i=1:N
    ddi(i) = c*dti0(i);
end

% init unknowns
% d0 - distance between the tag and anchor 0 (reference anchor)
% x,y,z - tag's coordinates
syms d0 x y z;

% [B]
% k coefficients
display('true distances');
di
% dti0*1E12
display('distance differences');
ddi
for i=1:N-1
    k(i) = 0.5*(ddi(i+1)^2 + (W(1,1)^2 + W(1,2)^2 + W(1,3)^2) - (W(i+1,1)^2 + W(i+1,2)^2 + W(i+1,3)^2));
end
k

% m,n coefficients
for i=1:N-1
    m(i) = 0;
    n(i) = 0;
    for j=1:N-1
        m(i) = m(i) + Ainverse(i,j)*k(j);
        n(i) = n(i) + Ainverse(i,j)*ddi(j+1);
    end
end
m
n
% Compute distance d0
% To do this an alfa, beta and gamma coefficients of the quadratic equation have to be computed

% alfa coefficient
alfa = -1;
for i=1:N-1
    alfa = alfa + n(i)^2;
end
alfa

% beta coefficient
beta = 0;
for i=1:N-1
    beta = beta + 2*(m(i)-W(1,i))*n(i);
end
beta

% gamma coefficient
gamma = 0;
for i=1:N-1
    gamma = gamma + (m(i)-W(1,i))^2;
end
gamma

% solve quadratic equation: alfa*d0^2 +beta*d0 + gamma = 0;
delta = beta^2 - 4*alfa*gamma;
d01 = (-beta + sqrt(delta))/(2*alfa);
d02 = (-beta - sqrt(delta))/(2*alfa);

% choose solution > 0
d0 = d02

% Compute coordinates of the tag:

% [X] = [A]^-1 [B]
display('Computed tag coordinates');
x = d0*n(1) + m(1)
y = d0*n(2) + m(2)
z = d0*n(3) + m(3)

%% OPTIONAL: calculate the 3D angles at which the tag is visible from each anchor

% Compute the distance from the tag to each anchor
for i=1:N
    di(i) = d0 + ddi(i);
end

% fi and gamma angles at which tag is visible from each anchor
for i=1:N
    fi(i) = atan((y-W(i,2))/(x-W(i,1)));
    gamma(i) = atan((z-W(i,3))/sqrt((x-W(i,1))^2 + (y-W(i,2))^2));
end
%% OPTIONAL

% Displacements from real location
display('Displacements from real location');
dx = abs(xt-x)
dy = abs(yt-y)
dz = abs(zt-z)

% Absolute displacement from real location
display('Absolute displacement from real location');
dp = sqrt(dx^2+dy^2+dz^2)


disp('press return to continue') 
pause() 