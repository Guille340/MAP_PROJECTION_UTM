function [zone,N,E,C,k] = utmDirect(lat,lon,ellipsoidName)

%  [zone,N,E,C,k] = universalTransverseMercator(lat,lon,ellipsoidName)
%
%  DESCRIPTION: Calculates the projected coordinates (Northing and Easting)
%  of a point from its geographic coordinates (latitude and longitude) and
%  reference ellipsoid using the Universal Transverse Mercator Projection 
%  (TMP). The function applies the standard formulas from the Defense Mapping 
%  Agency (DMA 8358.2, 1989). The grid origin (lat0,lon0), false Easting E0, 
%  false northing F0 and scale factor at the central meridian k0 are parameters
%  needed for calculation of the projected coordinates. These parameters are 
%  known and adopt the following values for the UTM projection: 
%
%     lat0 = 0°N
%     lon0 = central meridian (CM)
%     N0 = 10,000,000 m for South hemisphere and 0 m for North hemisphere
%     E0 = 500,000 m
%     k0 = 0.9996
%
%  INPUT VARIABLES
%  - lat: latitude of the point [deg]
%  - lon: longitude of the point [deg]
%  - ellipsoidName: name of reference ellipsoid. Check refEllip.m help for 
%    string options.
%
%  OUTPUT VARIABLES
%  - zone: three-character of the form 'XXY', with XX representing the
%    longitude zone (01-60) and Y representing the latitude zone (C-X, 
%    excluding zone I and O). For details on the limits of each zone refer
%    to the document DMA 8358.1 (1990).
%  - N: Northing [m]
%  - E: Easting [m]
%  - C: convergence of grid north to true north for selected point [deg]
%  - k: scale factor for selected point
%
%  INTERNALLY CALLED FUNCTIONS
%  - utmZone
%  - refEllip
%  - meridianarc
%
%  CONSIDERATIONS & LIMITATIONS
%  - Computations are accurate to the nearest 0.01 m.
%  - Latitude must be within 80°S and 84°N, limits of the UTM projection.
%  - The Northing is referred to the equator according to the UTM standard,
%    not to the latitudinal zone.
%
%  REFERENCES
%  - https://en.wikipedia.org/wiki/Universal_Transverse_Mercator_coordinate_system
%  - http://www.dmap.co.uk/utmworld.htm
%  - Hager,J. W., Fry,L. L., Jacks,S. S., Hill,D. R. (1990) Datums, Ellipsoids, 
%    Grids, and Grid Reference Systems, Defense Mapping Agency, Technical Manual 
%    TM8358.1.
%  - Hager,J. W., Behensky,J. F., and Drew,B. W. (1989) The Universal Grids: 
%    Universal Transverse Mercator (UTM) and Universal Polar Stereographic 
%    (UPS), Defense Mapping Agency, Technical Manual TM8358.2.

%  VERSION HISTORY
%  ===============
%  VERSION 1.0.0, 09 Jan 2020
%  Guillermo Jimenez Arranz
%  email: gjarranz@gmail.com
% ______________________________

latRadians = lat * pi/180;

% Error Control
if lat > 84 && lat < -80
    error('Latitude must be within the limits 84 and - 80 (84°N to 80°S)')
end

if lon < -180 && lon > 180
    error('Longitude must be within the limits -180 and 180 (180W to 180E)')
end

ellipsoidList = {'WGS84','WGS72','WGS66','WGS60','GRS80','NAD83','GDA94',...
    'AIR30','MdAIR','AusNS','AGD66','AGD84','INTER','IAU65','IAU68','GRS67',...
    'MdGRS','SAD69','CLK80','CLK66','NAD27','KRASO','ATS77','EVRST','BESSL'};
if ~ismember(ellipsoidName,ellipsoidList)
    error('Reference ellipsoid not recognised (see refEllip.m for valid names)')
end

% UTM Projection Parameters
zone = utmZone(lat,lon);
lonZoneWidth = 6; % width of UTM longitude zones [deg]
lonLowerLimit = -180; % longitude of UTM left boundary (zone 1) [deg]
lonZone = str2double(zone(1:2)); % longitude zone number
lat0 = 0; % Latitude of grid origin [deg]
lon0 = lonLowerLimit + lonZone*lonZoneWidth - lonZoneWidth/2; % CM lon [deg]
       ... NOTE: this same formula is applied to zones 32V, 31X and 37X.
lambda = (lon - lon0); % longitude difference [deg]
lambdaRadians = lambda * pi/180; % longitude difference [rad]
k0 = 0.9996; % Scale factor of Central Meridian (CM) for UTM projection
E0 = 5e5; % False Easting
N0 = 1e7 * (lat < 0); % False Northing

% Ellipsoid Parameters
[a,b,f] = refEllip(ellipsoidName); % parameters of the ref. ellipsoid
ec = sqrt(f*(2-f)); % eccentricity
ec2 = ec/(1-f); % second eccentricity
v = a/sqrt(1 - (ec*sin(latRadians))^2); % prime vertical curvature radius [m]                
M = meridianarc(lat,[a b],'UTM'); % meridian arc at point [m]
M0 = meridianarc(lat0,[a b],'UTM'); % meridian arc at grid origin [m]

% Pre-Stored Constants
sinLat = sin(latRadians);
cosLat = cos(latRadians);
cosLatPow2 = cosLat * cosLat;
cosLatPow3 = cosLatPow2 * cosLat;
cosLatPow4 = cosLatPow3 * cosLat;
cosLatPow5 = cosLatPow4 * cosLat;
cosLatPow6 = cosLatPow5 * cosLat;
cosLatPow7 = cosLatPow6 * cosLat;
cosLatPow8 = cosLatPow7 * cosLat;
tanLat = tan(latRadians);
tanLatPow2 = tanLat^2;
tanLatPow4 = tanLatPow2^2;
tanLatPow6 = tanLatPow2^3;
ec2Pow2 = ec2^2;
ec2Pow4 = ec2Pow2^2;
ec2Pow6 = ec2Pow2^3;
ec2Pow8 = ec2Pow4^2;

% Direct Formula (Constants)
T1 = M*k0;
T2 = v*k0*sinLat*cosLat/2;
T3 = v*k0*sinLat*cosLatPow3/24*(5 - tanLatPow2 + 9*ec2Pow2*cosLatPow2 ...
    + 4*ec2Pow4*cosLatPow4);
T4 = v*k0*sinLat*cosLatPow3/720*(61 - 58*tanLatPow2 + tanLatPow4 ...
    + 270*ec2Pow2*cosLatPow2 - 330*tanLatPow2*ec2Pow2*cosLatPow2 ...
    + 445*ec2Pow4*cosLatPow4 + 324*ec2Pow6*cosLatPow6 ...
    - 680*tanLatPow2*ec2Pow4*cosLatPow4 + 88*ec2Pow8*cosLatPow8 ...
    - 600*tanLatPow2*ec2Pow6*cosLatPow6 ...
    - 192*tanLatPow2*ec2Pow8*cosLatPow8);
T5 = v*k0*sinLat*cosLatPow7/40320*(1385 - 3111*tanLatPow2 + 543*tanLatPow4 ...
    - tanLatPow6);
T6 = v*k0*cosLat;
T7 = v*k0*cosLatPow3/6*(1 - tanLatPow2 + ec2Pow2*cosLatPow2);
T8 = v*k0*cosLatPow5/120*(5 - 18*tanLatPow2 + tanLatPow4 ...
    + 14*ec2Pow2*cosLatPow2 - 58*tanLatPow2*ec2Pow2*cosLatPow2 ...
    + 13*ec2Pow4*cosLatPow4 + 4*ec2Pow6*cosLatPow6 ...
    - 64*tanLatPow2*ec2Pow4*cosLatPow4 ...
    - 24*tanLatPow2*ec2Pow6*cosLatPow6);
T9 = v*k0*cosLatPow7/5040*(61 - 479*tanLatPow2 + 179*tanLatPow4 - tanLatPow6);

% Direct Formula
N = T1 + lambdaRadians^2*T2 + lambdaRadians^4*T3 + lambdaRadians^6*T4 ...
    + lambdaRadians^8*T5 + N0 - k0*M0; % corrected Northing (+)
E = lambdaRadians*T6 + lambdaRadians^3*T7 + lambdaRadians^5*T8 ...
    + lambdaRadians^7*T9 + E0; % corrected Easting (+)

% Convergence and Scale Factor
C = utmConvergence(lat,lon,ellipsoidName); % convergence at point [deg]
k = utmScaleFactor(lat,lon,ellipsoidName); % scale factor at point

