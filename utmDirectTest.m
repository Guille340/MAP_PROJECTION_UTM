%% TEST 01: Geographic to Grid Conversion
ellipsoidName = 'INTER'; % reference ellipsoid (International 1924)

% Tolerance
tol_utm = 0.01; % tolerance for Northing and Easting [m]
tol_C = 0.01/3600; % tolerance for convergence [deg]
tol_k = 1e-8; % tolerance for scale factor

% Test Inputs
lat_Dms = [73  0  0 ;  30  0  0 ;    72  4 32.110 ];
lon_Dms = [45  0  0 ; 102  0  0 ; -[113 54 43.321]];

% Expected Outputs
zone = [     '38X' ;      '48R' ;      '12X'];
N =    [8100702.90 ; 3322624.35 ; 8000000.01];
E =    [ 500000.00 ;  210577.93 ;  400000.00];
C_Dms = [0 0 0.00; -[1 30 3.76]; -[2 46 15.31]]; 
k = [0.99960000 ; 1.00063354 ; 0.99972228];

% Test
nTests = size(lat_Dms,1);
for m = 1:nTests
    lat = lat_Dms(m,1) + lat_Dms(m,2)/60 + lat_Dms(m,3)/3600;
    lon = lon_Dms(m,1) + lon_Dms(m,2)/60 + lon_Dms(m,3)/3600;
    [zone0,N0,E0,C0,k0] = utmDirect(lat,lon,ellipsoidName);
    C = C_Dms(m,1) + C_Dms(m,2)/60 + C_Dms(m,3)/3600;
    assert(strcmp(zone0,zone(m,:)))
    assert(abs(N0 - N(m)) <= tol_utm)
    assert(abs(E0 - E(m)) <= tol_utm)
    assert(abs(C0 - C) <= tol_C)   
    assert(abs(k0 - k(m)) <= tol_k)
end

