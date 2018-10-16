function [curvature] = Curvature(DEM,cellsize,r, varargin)
%%  
% This function calculates curvature of a DEM over radius r, which is
% rounded to the grid spaceing of the dem
%  Border cells are set to NaN 
% Author: Louis Sass

% inputs: 
% Z matrix of DEM Z values. This should be the first output from geotiffread or similar. Note the DEM must have square cells
% cellsize from the second output of geotiffread, xextent=yextent=cellsize
% r = radius of teh curvature calculation
% varargin = the method, described below:
%'area-averaged' computes the average angle of prominence for every pixel
%within a given radius. 
%'xy' computes the average of the maximum angle of prominence in both x and both y
%directions within the given radius. 
%For example, if the x's represent the terrain surface and O is the pixel
%of interest:

%           x                       x
%         x   x                   x   x
%xxxxxxxxx      xx       xxxxOxxxx      xx       xxxxxxx
%                  x    x                  x   x
%                     x                      x

% area-averaged yields 0 value of curvature 
% xy yields a negative (concave) value


%%

pnames = {   'method'  };
dflts =  {'area-averaged'};
[method] ...
    = internal.stats.parseArgs(pnames, dflts, varargin{:});

methodNames = {'area-averaged','xy'};
method = internal.stats.getParamVal(method,methodNames,'''method''');

%% Initiate individual grids for matrix calculation 

[nrows ncols] = size(DEM);
cx = DEM;
cy = DEM;

r = round(r/cellsize); %round to the nearest odd number..
%% Calculate Curvature

for i = 1:r
    ii = 2*i + 1;
    kk = 2*i;    
    dxp = zeros(nrows-kk,ncols-kk);
    nxp = zeros(nrows-kk,ncols-kk);
    dxm = zeros(nrows-kk,ncols-kk);
    nxm = zeros(nrows-kk,ncols-kk);
    dyp = zeros(nrows-kk,ncols-kk);
    nyp = zeros(nrows-kk,ncols-kk);
    dym = zeros(nrows-kk,ncols-kk);
    nym = zeros(nrows-kk,ncols-kk);
    switch method
            case 'area-averaged' %adds the cells from each quadrant 
            case 'xy'%sum of the maximum curvature in the x and y directions
                dxp = dxp-100; %to allow negative values in the max loop
                dxm = dxm-100;
                dyp = dyp-100;
                dym = dym-100;            
    end

    for jj = 1:ii-1
        switch method
            case 'area-averaged' %adds the cells from each quadrant 
                dxp = dxp + DEM(jj:(end-ii+jj),ii:end);
                nxp = nxp + ~isnan(DEM(jj:(end-ii+jj),ii:end));
                dxm = dxm + DEM(jj:(end-ii+jj),1:(end-ii+1));
                nxm = nxm + ~isnan(DEM(jj:(end-ii+jj),1:(end-ii+1)));
                dyp = dyp + DEM(1:(end-ii+1), jj:(end-ii+jj));
                nyp = nyp + ~isnan(DEM(1:(end-ii+1), jj:(end-ii+jj)));
                dym = dym + DEM(ii:end, jj:(end-ii+jj));
                nym = nym + ~isnan(DEM(ii:end, jj:(end-ii+jj)));
            case 'xy'
        end
    end
    
    switch method
            case 'area-averaged' %divides the sum to find the average value in each quadrant
                n = ii-1;
                cx(i+1:nrows-i, i+1:ncols-i) = (cx(i+1:nrows-i, i+1:ncols-i) + dxp./(nxp.*2) + dxm./(nxm.*2)) ./2;
                cy(i+1:nrows-i, i+1:ncols-i) = (cy(i+1:nrows-i, i+1:ncols-i) + dyp./(nyp.*2) + dym./(nym.*2)) ./2;
            case 'xy' %sum of the maximum curvature in the x and y directions
                dxp = max(dxp,(-DEM(i+1:nrows-i, i+1:ncols-i)+DEM(i+1:nrows-i,ii:end))./i); 
                dxm = max(dxm,(-DEM(i+1:nrows-i, i+1:ncols-i)+DEM(i+1:nrows-i,1:(end-kk)))./i);
                dyp = max(dyp,(-DEM(i+1:nrows-i, i+1:ncols-i)+DEM(1:(end-kk), i+1:ncols-i))./i);
                dym = max(dym,(-DEM(i+1:nrows-i, i+1:ncols-i)+DEM(ii:end, i+1:ncols-i))./i);
    end

end

switch method
    case 'area-averaged' 
        curvature = -atan((DEM - (cx+cy)./2) ./ cellsize );
    case 'xy'     
        cx(i+1:nrows-i, i+1:ncols-i) = +dxp./2 + dxm./2;
        cy(i+1:nrows-i, i+1:ncols-i) = +dyp./2 + dym./2;
        curvature = -atan(((cx+cy)./2) ./ cellsize ); % - sign so that + curvature is a convex slope and - curvature is a concave slope
end

% end
