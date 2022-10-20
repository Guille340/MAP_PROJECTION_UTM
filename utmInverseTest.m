%% TEST 01: Grid to Geographic Conversion
ellipsoidName = 'INTER'; % reference ellipsoid (International 1924)

% Tolerance
tol_geo = 0.01; % tolerance for latitude and longitude [deg]
tol_C = 0.01/3600; % tolerance for convergence [deg]
tol_k = 1e-8; % tolerance for scale factor

% Test Inputs
zone = [     '47R' ;      '31P' ;      '43X' ;      '31F'];
N    = [3322824.08 ; 1000000.00 ; 9000000.00 ; 4000329.42];
E    = [ 789411.59 ;  200000.00 ;  500000.00 ;  307758.89];

% Expected Outputs
lat_Dms = [ 30  0  6.489 ;  9  2 10.706 ; 81 3 30.487 ; -[54  6 28.992]];
lon_Dms = [101 59 59.805 ;  0 16 17.099 ; 75 0  0.000 ;    0  3 33.695 ];
C_Dms   = [  1 30  3.96  ;-[0 25 43.950];  0  0  0.000; -[ 2 22 58.830]];
k       = [   1.00063346 ;   1.00071386 ;  0.99960000 ;     1.00005345 ];

% Test
nTests = size(lat_Dms,1);
for m = 1:nTests
    [lat0,lon0,C0,k0] = utmInverse(zone(m,:),N(m),E(m),ellipsoidName);
    lat = lat_Dms(m,1) + lat_Dms(m,2)/60 + lat_Dms(m,3)/3600;
    lon = lon_Dms(m,1) + lon_Dms(m,2)/60 + lon_Dms(m,3)/3600;
    C = C_Dms(m,1) + C_Dms(m,2)/60 + C_Dms(m,3)/3600;
    assert(abs(lat0 - lat) <= tol_geo)
    assert(abs(lon0 - lon) <= tol_geo)
    assert(abs(C0 - C) <= tol_C)   
    assert(abs(k0 - k(m)) <= tol_k)
end

