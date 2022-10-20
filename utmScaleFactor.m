function k = utmScaleFactor(varargin)

%  k = utmScaleFactor(varargin)
%
%  DESCRIPTION: Calculates the scale factor at the working point. The scale
%  factor is a parameter that represents the distortion in distance of a point
%  in the projected grid. For the Central Meridian (CM) of the zone where the
%  specified location falls the scale factor for UTM is k0 = 0.9996e. If k < 1, 
%  the projected distance is smaller than the real (geographic) distance; for 
%  k > 1, the projected distance is larger than the real distance. No distortion
%  in distance occurs for k = 1, which for the UTM projection coincides with 
%  eastings 320 km and 680 km. Between those easting limits distance distortion 
%  is minimised (k0 < k < 1). The function uses the equations from the Defence 
%  Mapping Agency (1989).
%
%  INPUT VARIABLES (OPTION 1: Scale Factor from Geographic Coordinates)
%  - lat (varargin{1}): latitude of the point [deg]
%  - lon (varargin{2}): longitude of the point [deg]
%  - ellipsoidName (varargin{3}): name of reference ellipsoid. Check refEllip.m 
%    help for string options.
%
%  INPUT VARIABLES (OPTION 2: Scale Factor from Grid Coordinates)
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
%  - k: scale factor at selected location [deg]
%
%  INTERNALLY CALLED FUNCTIONS
%  - None
%
%  CONSIDERATIONS & LIMITATIONS
%  - The scale factor is k0 at the Central Meridian and 1 at E = 320,000 m and
%    E = 680,000 m.
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
        k = scaleFactorFromGeographic(varargin{1},varargin{2},varargin{3});
    case 4
        k = scaleFactorFromGrid(varargin{1},varargin{2},varargin{3},...
            varargin{4});
end

end

function k = scaleFactorFromGeographic(lat,lon,ellipsoidName)

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
    'MdGRS','SAD69','CLK80','CLK66','NAD27','KRASO','ATS77','EVRST''BESSL'};
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
k0 = 0.9996; % Scale factor of Central Meridian (CM) for UTM projection

% Ellipsoid Parameters
[~,~,f] = refEllip(ellipsoidName); % parameters of the ref. ellipsoid
ec = sqrt(f*(2-f)); % eccentricity
ec2 = ec/(1-f); % second eccentricity

% Pre-Stored Constants
cosLat = cos(latRadians);
cosLatPow2 = cosLat * cosLat;
cosLatPow4 = cosLatPow2 * cosLatPow2;
cosLatPow6 = cosLatPow4 * cosLatPow2;
tanLat = tan(latRadians);
tanLatPow2 = tanLat * tanLat;
tanLatPow4 = tanLatPow2 * tanLatPow2;
ec2Pow2 = ec2 * ec2;
ec2Pow4 = ec2Pow2 * ec2Pow2;
ec2Pow6 = ec2Pow4 * ec2Pow2;

% Direct Formula (Constants)
T26 = cosLatPow2/2 * (1 + ec2Pow2*cosLatPow2);
T27 = cosLatPow4/24 * (5 - 4*tanLatPow2 + 14*ec2Pow2*cosLatPow2 ...
    + 13*ec2Pow4*cosLatPow4 - 28*tanLatPow2*ec2Pow2*cosLatPow2 ...
    + 4*ec2Pow6*cosLatPow6 - 48*tanLatPow2*ec2Pow4*cosLatPow4 ...
    - 24*tanLatPow2*ec2Pow6*cosLatPow6);
T28 = cosLatPow6/720 * (61 - 148*tanLatPow2 + 16*tanLatPow4);

% Direct Formula
k = k0*(1 + lambdaRadians^2*T26 + lambdaRadians^4*T27 + lambdaRadians^6*T28);
end

function k = scaleFactorFromGrid(zone,N,E,ellipsoidName)

% Error Control
if ~ischar(zone) || length(zone)~=3
    error('Non-valid format for the specified geographic zone')
end
if ~isnumeric(N) || ~isnumeric(E)
    error('Northing and easting must be numeric values')
end
ellipsoidList = {'WGS84','WGS72','WGS66','WGS60','GRS80','NAD83','GDA94',...
    'AIR30','MdAIR','AusNS','AGD66','AGD84','INTER','IAU65','IAU68','GRS67',...
    'MdGRS','SAD69','CLK80','CLK66','NAD27','KRASO','ATS77','EVRST','BESSL'};
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
tanLatf = tan(latfRadians);
tanLatfPow2 = tanLatf * tanLatf;
ec2Pow2 = ec2 * ec2;
ec2Pow4 = ec2Pow2 * ec2Pow2;
ec2Pow6 = ec2Pow4 * ec2Pow2;

% Inverse Formula (Constants)
T29 = 1/(2*vf^2*k0^2) * (1 + ec2Pow2*cosLatfPow2);
T30 = 1/(24*vf^4*k0^4) * (1 + 6*ec2Pow2*cosLatfPow2 ...
    + 9*ec2Pow4*cosLatfPow4 + 4*ec2Pow6*cosLatfPow6 ...
    - 24*tanLatfPow2*ec2Pow4*cosLatfPow4 ...
    - 24*tanLatfPow2*ec2Pow6*cosLatfPow6);
T31 = 1/(720*vf^6*k0^6);

% Inverse Formula
k = k0*(1 + x^2*T29 + x^4*T30 + x^6*T31); % scale factor
end
