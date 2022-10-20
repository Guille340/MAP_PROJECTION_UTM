
% Shared Variables
tol = 0.01/3600; % tolerance for convergence [deg]
ellipsoidName = 'INTER'; % reference ellipsoid (International 1924)

%% TEST 01: Convergence from Geographic
% Test Inputs
lat_Dms = [73  0  0 ;  30  0  0 ;    72  4 32.110 ];
lon_Dms = [45  0  0 ; 102  0  0 ; -[113 54 43.321]];

% Expected Outputs
C_Dms = [0 0 0.00; -[1 30 3.76]; -[2 46 15.31]]; 

% Test
nTests = size(lat_Dms,1);
for m = 1:nTests
    lat = lat_Dms(m,1) + lat_Dms(m,2)/60 + lat_Dms(m,3)/3600;
    lon = lon_Dms(m,1) + lon_Dms(m,2)/60 + lon_Dms(m,3)/3600;
    C0 = utmConvergence(lat,lon,ellipsoidName);
    C = C_Dms(m,1) + C_Dms(m,2)/60 + C_Dms(m,3)/3600;
    assert(abs(C0 - C) <= tol)    
end

%% TEST 02: Convergence from Grid
% Test Inputs
zone = [     '47R' ;      '31P' ;      '43X' ;      '31F'];
N =    [3322824.08 ; 1000000.00 ; 9000000.00 ; 4000329.42];
E =    [ 789411.59 ;  200000.00 ;  500000.00 ;  307758.89];

% Expected Outputs
C_Dms   = [  1 30  3.96 ;-[0 25 43.950];  0  0  0.000; -[2 22 58.830]];

% Test
nTests = length(N);
for m = 1:nTests
    C0 = utmConvergence(zone(m,:),N(m),E(m),ellipsoidName);
    C = C_Dms(m,1) + C_Dms(m,2)/60 + C_Dms(m,3)/3600;
    assert(abs(C0 - C) <= tol)    
end
