function [lat,lon,C,k] = utmInverse(zone,N,E,ellipsoidName)

%  [lat,lon,C,k] = utmInverse(zone,N,E,ellipsoidName)
%
%  DESCRIPTION: Calculates the geographic coordinates (latitude and longitude) 
%  of a point from its projected coordinates (Easting and Northing), grid zone 
%  and reference ellipsoid using the Universal Transverse Mercator Projection 
%  (UTM). The function applies the standard formulas from the Defense Mapping 
%  Agency (DMA 8358.2, 1989). The grid origin (lat0,lon0), false Easting E0, 
%  false northing F0 and scale factor at the central meridian k0 are parameters
%  needed for calculation of the geographic coordinaes. These parameters are 
%  known and adopt the following values for the UTM projection: 
%
%     lat0 = 0°N
%     lon0 = central meridian (CM)
%     N0 = 10,000,000 m for South hemisphere and 0 m for North hemisphere
%     E0 = 500,000 m
%     k0 = 0.9996
%
%  INPUT VARIABLES
%  - zone: three-character of the form 'XXY', with XX representing the
%    longitude zone (01-60) and Y representing the latitude zone (C-X, 
%    excluding zone I and O). For details on the limits of each zone refer
%    to the document DMA 8358.1 (1990).
%  - N: Northing [m]
%  - E: Easting [m]
%  - ellipsoidName: name of reference ellipsoid. Check refEllip.m help for 
%    string options.
%
%  OUTPUT VARIABLES
%  - lat: latitude of the point [deg]
%  - lon: longitude of the point [deg]
%
%  INTERNALLY CALLED FUNCTIONS
%  - refEllip
%  - meridianarc
%  - meridianarcInverse
%
%  CONSIDERATIONS & LIMITATIONS
%  - Computations appear to be accurate to within 1 arcsec.
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

% UTM Zone Extraction
lonZone = str2double(zone(1:2));
latZone = zone(3); % 
northHemisphere = latZone >= 'N'; % true for North hemisphere (zones N-X)

% UTM Projection Parameters (General)
k0 = 0.9996; % Scale factor of Central Meridian (CM) for UTM projection
E0 = 5e5; % False Easting
N0 = 1e7 * (~northHemisphere); % False Northing

% Grid Origin Calculation
lonZoneWidth = 6; % width of UTM longitude zones [deg]
lonLowerLimit = -180; % longitude of UTM left boundary (zone 1) [deg]
lat0 = 0; % Latitude of grid origin [deg]
lon0 = lonLowerLimit + lonZone*lonZoneWidth - lonZoneWidth/2; % CM lon [deg]
lon0Radians = lon0 * pi/180; % Longitude of grid origin (i.e. CM lon) [rad]

% Ellipsoid Parameters
[a,b,f] = refEllip(ellipsoidName); % parameters of the ref. ellipsoid
ec = sqrt(f*(2-f)); % eccentricity
ec2 = ec/(1-f); % second eccentricity
M0 = meridianarc(lat0,[a b],'UTM'); % meridian arc at grid origin [m]

% Easting & Northing
y = N - N0 + k0*M0; % real Northing (+/-, relative to equator)
x = E - E0; % real Easting (+/-)

% Ellipsoid Parameters (Cont.)
Mf = y/k0; % distance from the equator to the footpoint [m]
latf = meridianarcInverse(Mf,[a b],'UTM'); % footpoint latitude [deg]     
latfRadians = latf*pi/180; % footpoint latitude [rad]
vf = a/sqrt(1 - (ec*sin(latfRadians))^2); % prime vertical curvature radius [m]                
pf = vf^3*(1-ec^2)/a^2; % meridian curvature radius [m]

% Pre-Stored Constants
cosLatf = cos(latfRadians);
cosLatfPow2 = cosLatf * cosLatf;
cosLatfPow4 = cosLatfPow2 * cosLatfPow2;
cosLatfPow6 = cosLatfPow4 * cosLatfPow2;
cosLatfPow8 = cosLatfPow6 * cosLatfPow2;
tanLatf = tan(latfRadians);
tanLatfPow2 = tanLatf * tanLatf;
tanLatfPow4 = tanLatfPow2 * tanLatfPow2;
tanLatfPow6 = tanLatfPow4 * tanLatfPow2;
ec2Pow2 = ec2 * ec2;
ec2Pow4 = ec2Pow2 * ec2Pow2;
ec2Pow6 = ec2Pow4 * ec2Pow2;
ec2Pow8 = ec2Pow6 * ec2Pow2;

% Direct Formula (Constants)
T10 = tanLatf/(2*pf*vf*k0^2);
T11 = tanLatf/(24*pf*vf^3*k0^4) * (5 + 3*tanLatfPow2 + ec2Pow2*cosLatfPow2 ...
    - 4*ec2Pow4*cosLatfPow4 - 9*tanLatfPow2*ec2Pow2*cosLatfPow2);
T12 = tanLatf/(720*pf*vf^5*k0^6) * (61 + 90*tanLatfPow2 ...
    + 46*ec2Pow2*cosLatfPow2 + 45*tanLatfPow4 ...
    - 252*tanLatfPow2*ec2Pow2*cosLatfPow2 - 3*ec2Pow4*cosLatfPow4 ...
    + 100*ec2Pow6*cosLatfPow6 - 66*tanLatfPow2*ec2Pow4*cosLatfPow4 ...
    - 90*tanLatfPow4*ec2Pow2*cosLatfPow2 + 88*ec2Pow8*cosLatfPow8 ...
    + 225*tanLatfPow4*ec2Pow4*cosLatfPow4 ...
    + 84*tanLatfPow2*ec2Pow6*cosLatfPow6 ...
    - 192*tanLatfPow2*ec2Pow8*cosLatfPow8);
T13 = tanLatf/(40320*pf*vf^7*k0^8) * (1385 + 3633*tanLatfPow2 ...
    + 4095*tanLatfPow4 + 1575*tanLatfPow6);
T14 = 1/(vf*k0*cosLatf);
T15 = 1/(6*vf^3*k0^3*cosLatf) * (1 + 2*tanLatfPow2 + ec2Pow2*cosLatfPow2);
T16 = 1/(120*vf^5*k0^5*cosLatf) * (5 + 6*ec2Pow2*cosLatfPow2 ...
    + 28*tanLatfPow2 - 3*ec2Pow4*cosLatfPow4 ...
    + 8*tanLatfPow2*ec2Pow2*cosLatfPow2 + 24*tanLatfPow4 ...
    - 4*ec2Pow6*cosLatfPow6 + 4*tanLatfPow2*ec2Pow4*cosLatfPow4 ...
    + 24*tanLatfPow2*ec2Pow6*cosLatfPow6);
T17 = 1/(5040*vf^7*k0^7*cosLatf) * (61 + 662*tanLatfPow2 + 1320*tanLatfPow4 ...
    + 720*tanLatfPow6);

% Inverse Formula
latRadians = latfRadians - x^2*T10 + x^4*T11 - x^6*T12 + x^8*T13;
lonRadians = x*T14 - x^3*T15 + x^5*T16 - x^7*T17 + lon0Radians;
lat = latRadians * 180/pi; % latitude of point [deg]
lon = lonRadians * 180/pi; % longitude of point [deg]

% Convergence and Scale Factor
C = utmConvergence(zone,N,E,ellipsoidName); % convergence at point [deg]
k = utmScaleFactor(zone,N,E,ellipsoidName); % scale factor at point
