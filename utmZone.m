function zone = utmZone(lat,lon)

%  zone = utmZone(lat,lon)
%
%  DESCRIPTION: returns the grid zone of a given geographic point as a three-
%  character string of the form 'XXY', where XX is the longitude zone and 'Y'
%  the latitude zone. There are 60 longitudinal zones numbered 1 to 60, 
%  starting at 180°W. The longitudinal zones are 6° wide, except for a few 
%  exceptions around Norway and Svalbard. The projected grid is also split in
%  20 latitude zones named C to X, excluding letters I and O. All latitudinal
%  zones are 8° south-north, apart from zone X which is 12° south-north.
%
%  INPUT VARIABLES
%  - lat: latitude of the point [deg]
%  - lon: longitude of the point [deg]
%
%  OUTPUT VARIABLES
%  - zone: three-character string of the form 'XXY', with XX representing the
%    longitudinal zone (01-60) and Y representing the latitudinal zone (C-X). 
%    For details on zone limits refer to the document DMA 8358.1 (1990).
%
%  INTERNALLY CALLED FUNCTIONS
%  - None
%
%  REFERENCES
%  - https://en.wikipedia.org/wiki/Universal_Transverse_Mercator_coordinate_system
%  - http://www.dmap.co.uk/utmworld.htm
%  - Hager,J. W., Fry,L. L., Jacks,S. S., Hill,D. R. (1990) Datums, Ellipsoids, 
%    Grids, and Grid Reference Systems, Defense Mapping Agency, Technical Manual 
%    TM8358.1.

%  VERSION HISTORY
%  ===============
%  VERSION 1.0.0, 09 Jan 2020
%  Guillermo Jimenez Arranz
%  email: gjarranz@gmail.com
% ______________________________

% UTM Longitude Zone Calculation
lonZoneWidth = 6; % width of UTM longitude zones [deg]
lonLowerLimit = -180; % longitude of UTM left boundary (zone 1) [deg]
lonZone = floor((lon - lonLowerLimit)/lonZoneWidth) + 1;

% UTM Latitude Zone Calculation
latZoneWidth = 8; % width of UTM latitude zones between 80S and 72N [deg]
latLowerLimit = -80; % latitude of UTM bottom boundary (zone C) [deg]
latZoneNames = char([67:72, 74:78, 80:88]); % name of latitude zones (C-X)
if lat < 72
    ilatZone = floor((lat - latLowerLimit)/latZoneWidth) + 1;
    latZone = latZoneNames(ilatZone);
else
    latZone = latZoneNames(end); % zone X (72N to 84N, 12 deg width)
end

% Correct UTM Longitude Zone for Latitude Zones V and X
if strcmp(latZone,'V')
    if lon >= 3 && lon < 12
        lonZone = 32;
    end
elseif strcmp(latZone,'X')
    if lon >= 0 && lon < 9
        lonZone = 31;
    elseif lon >= 9 && lon < 21
        lonZone = 33;
    elseif lon >= 21 && lon < 33
        lonZone = 35;
    elseif lon >= 33 && lon < 42
        lonZone = 37;
    end
end

% Zone Name
zone = sprintf('%0.2d%s',lonZone,latZone);