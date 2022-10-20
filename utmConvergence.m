function C = utmConvergence(varargin)

%  C = utmConvergence(varargin)
%
%  DESCRIPTION: Calculates the covergence at the working point, i.e. the angle
%  in degrees between the true north and the projected north at that location.
%  The point can be given in geographic coordinates (latitude,longitude) or
%  as grid coordinates (northing,easting). In the last case the zone must
%  also be specified. The function uses the equations from the Defence Mapping 
%  Agency (1989).
%
%  INPUT VARIABLES (OPTION 1: Convergence from Geographic Coordinates)
%  - lat (varargin{1}): latitude of the point [deg]
%  - lon (varargin{2}): longitude of the point [deg]
%  - ellipsoidName (varargin{3}): name of reference ellipsoid. Check refEllip.m 
%    help for string options.
%
%  INPUT VARIABLES (OPTION 2: Convergence from Grid Coordinates)
%  - zone (varargin{1}): three-character of the form 'XXY', with XX 
%    representing the longitude zone (01-60) and Y representing the latitude 
%    zone (C-X, excluding zone I and O). For details on the limits of each zone 
%    refer to the document DMA 8358.1 (1990).
%  - N (varargin{2}): Northing [m]
%  - E (varargin{3}): Easting [m]
%  - ellipsoidName (varargin{4}): name of reference ellipsoid. Check refEllip.m 
%    help for string options.
%
%  OUTPUT VARIABLES
%  - C: convergence at selected location [deg]
%
%  INTERNALLY CALLED FUNCTIONS
%  - None
%
%  CONSIDERATIONS & LIMITATIONS
%  - The convergence is 0 deg at the Central Meridian used for the grid, where
%    the projection causes no distortion.
%  - A simple way to calculate the convergence from a point given in geographic 
%    coordinates is to use the following equation: atan(tan(lon-lon0)*sin(lat), 
%    with lon0 being the longitude of the grid origin or Central Meridian.
%
%  REFERENCES
%  - https://en.wikipedia.org/wiki/Universal_Transverse_Mercator_coordinate_system
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

switch nargin
    case 3
        C = convergenceFromGeographic(varargin{1},varargin{2},varargin{3});
    case 4
        C = convergenceFromGrid(varargin{1},varargin{2},varargin{3},...
            varargin{4});
end
end

function C = convergenceFromGeographic(lat,lon,ellipsoidName)

latRadians = lat * pi/180;

% Error Control
if ~isnumeric(lat) || ~isnumeric(lon)
    error('Latitude and longitude must be numeric values')
end
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
lon0 = lonLowerLimit + lonZone*lonZoneWidth - lonZoneWidth/2; % CM lon [deg]
       ... NOTE: this same formula is applied to zones 32V, 31X and 37X.
lambda = (lon - lon0); % longitude difference [deg]
lambdaRadians = lambda * pi/180; % longitude difference [rad]

% Ellipsoid Parameters
[~,~,f] = refEllip(ellipsoidName); % parameters of the ref. ellipsoid
ec = sqrt(f*(2-f)); % eccentricity
ec2 = ec/(1-f); % second eccentricity

% Pre-Stored Constants
sinLat = sin(latRadians);
cosLat = cos(latRadians);
cosLatPow2 = cosLat * cosLat;
cosLatPow4 = cosLatPow2 * cosLatPow2;
cosLatPow6 = cosLatPow4 * cosLatPow2;
cosLatPow8 = cosLatPow6 * cosLatPow2;
tanLat = tan(latRadians);
tanLatPow2 = tanLat * tanLat;
tanLatPow4 = tanLatPow2 * tanLatPow2;
ec2Pow2 = ec2 * ec2;
ec2Pow4 = ec2Pow2 * ec2Pow2;
ec2Pow6 = ec2Pow4 * ec2Pow2;
ec2Pow8 = ec2Pow6 * ec2Pow2;

% Direct Formula (Constants)
T18 = sinLat;
T19 = sinLat*cosLatPow2/3 * (1 + 3*ec2Pow2*cosLatPow2 ...
    + 2*ec2Pow4*cosLatPow4);
T20 = sinLat*cosLatPow4/15 * (2 - tanLatPow2 + 15*ec2Pow2*cosLatPow2 ...
    + 35*ec2Pow4*cosLatPow4 - 15*tanLatPow2*ec2Pow2*cosLatPow2 ...
    + 33*ec2Pow6*cosLatPow6 - 50*tanLatPow2*ec2Pow4*cosLatPow4 ...
    + 11*ec2Pow8*cosLatPow8 - 60*tanLatPow2*ec2Pow6*cosLatPow6 ...
    - 24*tanLatPow2*ec2Pow8*cosLatPow8);
T21 = sinLat*cosLatPow6/315 * (17 - 26*tanLatPow2 + 2*tanLatPow4);

% Inverse Formula
CRadians = lambdaRadians*T18 + lambdaRadians^3*T19 + lambdaRadians^5*T20 ...
    + lambdaRadians^7*T21; % convergence [rad]
C = CRadians * 180/pi; % convergence [deg]
end

function C = convergenceFromGrid(zone,N,E,ellipsoidName)

% Error Control
if ~ischar(zone) || length(zone)~=3
    error('Non-valid format for the specified geographic zone')
end
if ~isnumeric(N) || ~isnumeric(E)
    error('Northing and easting must be numeric values')
end
ellipsoidList = {'WGS84','WGS72','WGS66','WGS60','GRS80','NAD83','GDA94',...
    'AIR30','MdAIR','AusNS','AGD66','AGD84','INTER','IAU65','IAU68','GRS67',...
    'MdGRS','SAD69','CLK80','CLK66','NAD27','KRASO','ATS77','EVRST''BESSL'};
if ~ismember(ellipsoidName,ellipsoidList)
    error('Reference ellipsoid not recognised (see refEllip.m for valid names)')
end

% UTM Zone
latZone = zone(3); % latitude zone (C-X, except I and O)
northHemisphere = latZone >= 'N'; % true for North hemisphere (zones N-X)

% UTM Projection Parameters (General)
k0 = 0.9996; % Scale factor of Central Meridian (CM) for UTM projection
E0 = 5e5; % False Easting
N0 = 1e7 * (~northHemisphere); % False Northing
lat0 = 0; % Latitude of grid origin [deg]

% Ellipsoid Parameters
[a,b,f] = refEllip(ellipsoidName); % parameters of the ref. ellipsoid
ec = sqrt(f*(2-f)); % eccentricity
ec2 = ec/(1-f); % second eccentricity
M0 = meridianarc(lat0,[a b],'UTM'); % meridian arc at grid origin [m]

% Easting & Northing
y = N - N0 + k0*M0; % real Northing (+/-, relative to equator)
x = E - E0; % real Easting (+/-)

% Ellipsoid Parameters (Cont.)
Mf = abs(y/k0); % distance from the equator to the footpoint [m]
latf = meridianarcInverse(Mf,[a b],'UTM'); % footpoint latitude [deg]     
latfRadians = latf*pi/180; % footpoint latitude [rad]
vf = a/sqrt(1 - (ec*sin(latfRadians))^2); % prime vertical curvature radius [m]                

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

% Inverse Formula (Constants)
T22 = tanLatf/(vf*k0);
T23 = tanLatf/(3*vf^3*k0^3) * (1 + tanLatfPow2 - ec2Pow2*cosLatfPow2 ...
    - 2*ec2Pow4*cosLatfPow4);
T24 = tanLatf/(15*vf^5*k0^5) * (2 + 5*tanLatfPow2 + 2*ec2Pow2*cosLatfPow2 ...
    + 3*tanLatfPow4 + tanLatfPow2*ec2Pow2*cosLatfPow2 ...
    + 9*ec2Pow4*cosLatfPow4 + 20*ec2Pow6*cosLatfPow6 ...
    - 7*tanLatfPow2*ec2Pow4*cosLatfPow4 ...
    - 27*tanLatfPow2*ec2Pow6*cosLatfPow6 + 11*ec2Pow8*cosLatfPow8 ...
    - 24*tanLatfPow2*ec2Pow8*cosLatfPow8);
T25 = tanLatf/(315*vf^7*k0^7) * (17 + 77*tanLatfPow2 + 105*tanLatfPow4 ...
    + 45*tanLatfPow6);

% Inverse Formula
CRadians = x*T22 - x^3*T23 + x^5*T24 - x^7*T25; % convergence [rad]
C = CRadians * 180/pi; % convergence [deg]
end