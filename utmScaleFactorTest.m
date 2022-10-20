
% Shared Variables
tol = 1e-8; % tolerance for scale factor
ellipsoidName = 'INTER'; % reference ellipsoid (International 1924)

%% TEST 01: Scale Factor from Geographic
% Test Inputs
lat_Dms = [73  0  0 ;  30  0  0 ;    72  4 32.110 ];
lon_Dms = [45  0  0 ; 102  0  0 ; -[113 54 43.321]];

% Expected Outputs (from Geographic)
k = [0.99960000 ; 1.00063354 ; 0.99972228];

% Test
nTests = size(lat_Dms,1);
for m = 1:nTests
    lat = lat_Dms(m,1) + lat_Dms(m,2)/60 + lat_Dms(m,3)/3600;
    lon = lon_Dms(m,1) + lon_Dms(m,2)/60 + lon_Dms(m,3)/3600;
    k0 = utmScaleFactor(lat,lon,ellipsoidName);
    assert(abs(k0 - k(m)) <= tol)    
end

%% TEST 02: Scale Factor from Grid
% Test Inputs
zone = [     '47R' ;      '31P' ;      '43X' ;      '31F'];
N =    [3322824.08 ; 1000000.00 ; 9000000.00 ; 4000329.42];
E =    [ 789411.59 ;  200000.00 ;  500000.00 ;  307758.89];

% Expected Outputs
k = [1.00063346 ; 1.00071386 ; 0.99960000 ; 1.00005345];

% Test
nTests = length(N);
for m = 1:nTests
    k0 = utmScaleFactor(zone(m,:),N(m),E(m),ellipsoidName);
    assert(abs(k0 - k(m)) <= tol)    
end



