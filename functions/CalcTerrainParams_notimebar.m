function [slope,aspect,eastness,northness] = CalcTerrainParams_notimebar(DEM,cellsize)
%% 
% This function calculates slopes of a DEM using the same methodology as
% ArcGIS. Border cells are set to NaN 
% Author: Christopher Nuth 

%% Initiate individual grids for matrix calculation 

[nrows ncols] = size(DEM);
    slope = nan(size(DEM));
    aspect = nan(size(DEM));
    eastness = nan(size(DEM));
    northness = nan(size(DEM));


%% Calculate Slope


dx = ( (DEM(1:end-2,3:end) + 2.*DEM(2:end-1,3:end) + DEM(3:end,3:end)) - (DEM(1:end-2,1:end-2) + 2.*DEM(2:end-1,1:end-2) + DEM(3:end,1:end-2)) ) / (8 .* cellsize);
dy = ( (DEM(1:end-2,1:end-2) + 2.*DEM(1:end-2,2:end-1) + DEM(1:end-2,3:end)) - (DEM(3:end,1:end-2) + 2.*DEM(3:end,2:end-1) + DEM(3:end,3:end)) ) / (8 .* cellsize);

slope(2:nrows-1,2:ncols-1) = atand(sqrt( dx.^2 + dy.^2 ));

aspect(2:nrows-1,2:ncols-1) = 180 -  atand(dy./dx) + 90.*(dx./abs(dx));

eastness(2:nrows-1,2:ncols-1) = sind(-atand(dx));

northness(2:nrows-1,2:ncols-1) = sind(-atand(dy));

% end
