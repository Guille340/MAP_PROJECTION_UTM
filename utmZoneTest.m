%% TEST 01: Zone Format
% Test Inputs
nTests = 5;
lat = rand(1,nTests)*164 -  80;
lon = rand(1,nTests)*360 - 180;
for m = 1:nTests
    zone0 = utmZone(lat(m),lon(m));
    assert(ischar(zone0) && isequal(size(zone0),[1 3]))
end

%% TEST 02: Left Border
% Test Inputs
lat = [ -80;    0;  -29;   11];
lon = [-180; -180; -177; -177];

% Expected Outputs
zone = ['01C';'01N';'01J';'01P'];

% Test
nTests = length(lat);
for m = 1:nTests
    zone0 = utmZone(lat(m),lon(m));
    assert(strcmp(zone0,zone(m,:)))
end

%% TEST 03: Origin of UTM Grid
% Test Inputs
lat = 0;
lon = 0;

% Expected Outputs
zone = ['31N';'30M'];

% Test
zone0 = utmZone(lat,lon);
assert(strcmp(zone0,zone(1,:)) || strcmp(zone0,zone(2,:)))

%% Test 04: Cells of Special Width and Height
% Test Inputs
lat = [78; 78; 78; 78; 78; 78; 60];
lon = [7 ; 10; 19; 22; 31; 34;  4];

% Expected Outputs
zone = ['31X';'33X';'33X';'35X';'35X';'37X';'32V'];

% Test
nTests = length(lat);
for m = 1:nTests
    zone0 = utmZone(lat(m),lon(m));
    assert(strcmp(zone0,zone(m,:)))  
end
