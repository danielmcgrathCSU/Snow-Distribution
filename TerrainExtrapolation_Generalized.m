%% Terrain Analysis
% This script analyzes the surface geometry and creates a SWE model based
% on terrain parameters. Extrapolates SWE based on a stepwise linear regression and regression tree.

%inputs: 
% dem in geotiff format with equal X,Y,Z units, like a UTM or Alberts
% coordinate system
% extent geotiff format and snapped to dem
% 
% Louis Sass - written 2016.06.01
% Dan McGrath - edited 

clear 
tic
addpath data/  
addpath functions/
dbstop if error

tic
%% set evaluation parameters here
radius = 500; % set radius over which to evaluate curvature, in the same units as the dem cellsize
distance = [10,600]; % set distance range for shelter term Sx, where minimum must be >= cellsize
near_distance = [10,20]; % set distance ranges for slopebreak term Sb
far_distance = [300,2000];
direction = 30; % set direction for prevailing wind direction for shelter term Sx 

dem = ('2015.08.13.WolvDEMreg_10m.tif');%DEM goes here
mask = ('2015mask_10m.tif');%Extent goes here
pick = ('Wolv_2017_TA.csv');% Radar Pick file goes here
%this pick file contains 10 columns: 1 - trace, 2 - long, 3 - lat, 4 -
%elevation, 5 - twtt, 6 - thickness, 7 - SWE, 8 - profile line, 9 - x
%(UTM), 10 - y (UTM)

findOptimumCurvature = ('no'); % set to yes to find best curvature values, then turn off and re-run using those values above
findOptimumSx = ('no'); % set to yes to find best Sx values, then turn off and re-run using those values above
plots = ('off'); %turn on extra plots for troubleshooting
isHelicopter = ('no'); %set to 'no' for ground data, 'yes' for mixed or helicopter data
maxslope = 90; %set the slope mask here to ditch steep terrain data

names = {'Z','Sb','Sx','curvature','northness','eastness','slope'}; %Possible names: 'Z', 'slope', 'aspect', 'northness', 'eastness', 'curvature','Sx','Sb'
labels = {'elevation [1000 m]' 'slopebreak' 'shelter' 'curvature' 'northness' 'eastness' 'slope [{\circ}]'};

%% import and format the data

[DEM.Z, DEM.ref] = geotiffread(dem); 
if DEM.ref.CellExtentInWorldX == DEM.ref.CellExtentInWorldY
    DEM.cellsize = DEM.ref.CellExtentInWorldX;
else
    error('DEM should be on a square grid')
end
DEM.Z(DEM.Z==DEM.Z(1,1))=NaN; %DEM.Z = the DEM value map
DEM.info = geotiffinfo(dem);
DEM.R = DEM.info.RefMatrix;
DEM.x = repmat((DEM.ref.XWorldLimits(1) + DEM.ref.CellExtentInWorldX/2:DEM.ref.CellExtentInWorldX:DEM.ref.XWorldLimits(2) - DEM.ref.CellExtentInWorldX/2),DEM.ref.RasterSize(1),1);
DEM.y = repmat(flipud((DEM.ref.YWorldLimits(1) + DEM.ref.CellExtentInWorldY/2:DEM.ref.CellExtentInWorldY:DEM.ref.YWorldLimits(2) - DEM.ref.CellExtentInWorldY/2)'),1,DEM.ref.RasterSize(2));

[GL.on, GL.ref] = geotiffread(mask); %GL.on = glacier extent mask
if GL.ref.CellExtentInWorldX == GL.ref.CellExtentInWorldY
    GL.cellsize = GL.ref.CellExtentInWorldX;
else
    error('Extent should be on a square grid')
end
GL.info = geotiffinfo(mask);
GL.R = GL.info.RefMatrix;

%importing GPR pick data
SWE.data = readtable(pick); 

%get x and y locations for each GPR point observation
[SWE.data.x SWE.data.y] = projfwd(DEM.info,SWE.data.lat, SWE.data.long);
swemax = max(SWE.data.SWE); %Max pick value in GPR

SWE.data=SWE.data(~any(ismissing(SWE.data),2),:);


%% check inputs

if GL.cellsize ~= DEM.cellsize 
    error('inputs should have the same cellsize')
end  
        
%% evaluate terrain

%calculating terrain parameters
[DEM.slope,DEM.aspect,DEM.eastness,DEM.northness] = CalcTerrainParams_notimebar(DEM.Z,DEM.cellsize);
DEM.shade = 2.*DEM.eastness + DEM.northness;

%%
% calculating curvature');
[DEM.curvature] = Curvature(DEM.Z,DEM.cellsize,radius, 'method', 'xy'); %method 'xy' is generally much more useful.

% calculating Sx
temp_sx = Sx(DEM.Z, DEM.cellsize, distance(1),distance(2), direction);
[DEM.Sx] = temp_sx;

% calculate Sb
temp_sb = Sb(DEM.Z, DEM.cellsize, near_distance,far_distance, direction);
temp_sb(temp_sb<=0) = 0; % +slopebreaks are well correlated with drifts, -slopebreaks are not necessarily wind scoured, the correlation here is better if you restrict Sb to + values 
[DEM.Sb] = temp_sb;

%% terrain parameter figures for the whole DEM
do=1;
if do==1
    
for n = 1:length(names)
figure ();
colormap(jet)
if strcmp(names{n},'curvature')==1
    clims = [-0.5 0.5];
    imagesc(DEM.(names{n}), 'alphadata', ~isnan(DEM.Z),clims);hold on
else
    imagesc(DEM.(names{n}), 'alphadata', ~isnan(DEM.Z));hold on
end
axis ij;
axis image;   
colorbar;
set(gca,'Fontsize',24);
xlabel('east');
ylabel('north');
%text(100, 100, labels{n});
title(labels{n});
end
end

%%
%Subsequent sections of the code iterate the MVR and regression tree analysis many times through different training and test datasets.  
%Certain calculations are necessary on the entire dataset and that is what
%follows. 

%make figure for that will contain the 100 subplots showing the various
%transects that were used in the calculation
figure, hold on, grid on, box on

%Compile DEM into vectors and ditch NaN locations
Xtemp = [DEM.x(:),DEM.y(:),DEM.Z(:)];
index = ~isnan(Xtemp(:,3));

%now full dataset
Xi = [SWE.data.x,SWE.data.y];
X = double(Xtemp(index,1:2));
T = delaunayn(X);

gprIndicies = dsearchn(X,T,Xi);
gprIndicies = gprIndicies(isnan(gprIndicies)==0);

%entire dataset - for the purposes of standardizing values
aggregate = table;
temp = accumarray(gprIndicies,SWE.data.elev,[],@median);
aggregate.Z = temp(temp~=0);
temp = accumarray(gprIndicies,gprIndicies,[],@median);
aggregate.index = temp(temp~=0);
temp = accumarray(gprIndicies,SWE.data.SWE,[],@median);
aggregate.SWE = temp(temp~=0);

temp = accumarray(gprIndicies,SWE.data.long,[],@median);
aggregate.long = temp(temp~=0);
temp = accumarray(gprIndicies,SWE.data.lat,[],@median);
aggregate.lat = temp(temp~=0);
temp = accumarray(gprIndicies,SWE.data.x,[],@median);
aggregate.x = temp(temp~=0);
temp = accumarray(gprIndicies,SWE.data.y,[],@median);
aggregate.y = temp(temp~=0);
temp = accumarray(gprIndicies,SWE.data.elev,[],@median);
aggregate.elev= temp(temp~=0);

for n = 1:length(names)
   %'calculating terrain parameters at GPR data locations']);
    allvalues = DEM.(names{n});
    vectorvalues = allvalues(:);
    nonanvalues = vectorvalues(index,1);
    aggregate.(names{n}) = nonanvalues(aggregate.index);
end

aggregate=aggregate(~any(ismissing(aggregate),2),:);

ind_entire = find(aggregate.slope<=maxslope);
SWE_entire_c = ((aggregate.x(ind_entire) - DEM.ref.XWorldLimits(1))./DEM.cellsize) + 1/2; %interp
SWE_entire_r = ((DEM.ref.YWorldLimits(2) - aggregate.y(ind_entire))./DEM.cellsize) + 1/2;

clear xySWE
xySWE = aggregate.SWE ./ cosd(aggregate.slope); %corrects for slope
[NormalxySWE SWEmu_full SWEsigma_full] = zscore(xySWE); % standardized
TerrainTable = table(NormalxySWE(ind_entire));
TerrainTable.Properties.VariableNames{'Var1'} = 'SWE';

for n= 1:length(names)
    [add mu sigma] = zscore(aggregate.(names{n})); %standardized
    TerrainTable.(names{n}) = add(ind_entire);
    muV.(names{n}) = mu;
    sigmaV.(names{n})= sigma;
end

%%
%These values/calculations get used in the main "for loop" but don't need to be calculated each
%time
%masking out the glacier data

[m,n] = size(GL.on);
GL.aligned = zeros(m,n);
GL.aligned(GL.on==0) = 1;
addx = (GL.ref.XWorldLimits(1) - DEM.ref.XWorldLimits(1))/DEM.ref.CellExtentInWorldX;
add = zeros(m,addx).*255;
GL.aligned = [add,GL.aligned];
addx = (DEM.ref.XWorldLimits(2) - GL.ref.XWorldLimits(2))/DEM.ref.CellExtentInWorldX;
add = zeros(m,addx).*255;
GL.aligned = [GL.aligned,add];
[m,n] = size(GL.aligned);
addy = (DEM.ref.YWorldLimits(2) - GL.ref.YWorldLimits(2))/DEM.ref.CellExtentInWorldY;
add = zeros(addy,n).*255;
GL.aligned = [add;GL.aligned];
addy = (GL.ref.YWorldLimits(1) - DEM.ref.YWorldLimits(1))/DEM.ref.CellExtentInWorldY;
add = zeros(addy,n).*255;
GL.aligned = [GL.aligned;add];
if size(GL.aligned)~=size(DEM.Z)
    error('the mask is not aligned with the data')
end

GL.area = nansum(nansum(GL.aligned));

%%
P=[];
itr=5; % # of iterations (can be 5, 10, 50, 100 with how it is currently written)

%choose random starting point
pp=randi([1 height(SWE.data)],itr,1);

%get "divisor" for training:test data split
%div= 3.33; % this is for a 70:30 split
div=4; % this is for a 75:25 split
%div=5; %this is for a 80:20 split

for i=1:itr

%find random number within the SWE.data table length
clear p SWE_test SWE_train
p=pp(i,1);

if  p>(height(SWE.data)-(height(SWE.data)/div))
    p=p:(height(SWE.data));
    p=p';
    temp=round((height(SWE.data)/div)-length(p));
    p=[p; (1:max(temp))'];
    p=sort(p);

SWE_test=SWE.data(p,:);   % extract rows in the p vector
SWE_train=SWE.data;
SWE_train(p,:) = [];

else
p=p:p+(height(SWE.data)/div);
p=p';

%partition training table by this random permutation
SWE_test=SWE.data(p,:);   % extract rows in the p vector
SWE_train=SWE.data;
SWE_train(p,:) = [];

end
%% 
clear Xtemp index Xi X T gprIndicies_test gprIndicies_train
Xtemp = [DEM.x(:),DEM.y(:),DEM.Z(:)];
index = ~isnan(Xtemp(:,3));

if strcmp(isHelicopter,'yes')

X = double(Xtemp(index,:));
T = delaunayn(X);

Xi = [SWE.data.x,SWE.data.y,(SWE.data.elev + 12)];% if GPR output changes to ellipsoid heights then remove the +12
gprIndicies = dsearchn(X,T,Xi); 
else
% 'extracting XY positions
X = double(Xtemp(index,1:2));
T = delaunayn(X);

%test dataset
Xi = [SWE_test.x,SWE_test.y];
gprIndicies_test = dsearchn(X,T,Xi);
clear Xi

%now training dataset
Xi = [SWE_train.x,SWE_train.y];
gprIndicies_train = dsearchn(X,T,Xi);
clear Xi
end

gprIndicies_test = gprIndicies_test(isnan(gprIndicies_test)==0);
gprIndicies_train = gprIndicies_train(isnan(gprIndicies_train)==0);

%% aggregate the SWE observations
%training dataset
aggregate_train = table;
temp = accumarray(gprIndicies_train,SWE_train.elev,[],@median);
aggregate_train.Z = temp(temp~=0);
temp = accumarray(gprIndicies_train,gprIndicies_train,[],@median);
aggregate_train.index = temp(temp~=0);
temp = accumarray(gprIndicies_train,SWE_train.SWE,[],@median);
aggregate_train.SWE = temp(temp~=0);

%ag_SWE_train=table;
temp = accumarray(gprIndicies_train,SWE_train.long,[],@median);
aggregate_train.long = temp(temp~=0);
temp = accumarray(gprIndicies_train,SWE_train.lat,[],@median);
aggregate_train.lat = temp(temp~=0);
temp = accumarray(gprIndicies_train,SWE_train.x,[],@median);
aggregate_train.x = temp(temp~=0);
temp = accumarray(gprIndicies_train,SWE_train.y,[],@median);
aggregate_train.y = temp(temp~=0);
temp = accumarray(gprIndicies_train,SWE_train.elev,[],@median);
aggregate_train.elev= temp(temp~=0);

%test dataset
aggregate_test = table;
temp = accumarray(gprIndicies_test,SWE_test.elev,[],@median);
aggregate_test.Z = temp(temp~=0);
temp = accumarray(gprIndicies_test,gprIndicies_test,[],@median);
aggregate_test.index = temp(temp~=0);
temp = accumarray(gprIndicies_test,SWE_test.SWE,[],@median);
aggregate_test.SWE = temp(temp~=0);

temp = accumarray(gprIndicies_test,SWE_test.long,[],@median);
aggregate_test.long = temp(temp~=0);
temp = accumarray(gprIndicies_test,SWE_test.lat,[],@median);
aggregate_test.lat = temp(temp~=0);
temp = accumarray(gprIndicies_test,SWE_test.x,[],@median);
aggregate_test.x = temp(temp~=0);
temp = accumarray(gprIndicies_test,SWE_test.y,[],@median);
aggregate_test.y = temp(temp~=0);
temp = accumarray(gprIndicies_test,SWE_test.elev,[],@median);
aggregate_test.elev= temp(temp~=0);

%% extract terrain parameters associated with each GPR point

for n = 1:length(names)
   %'calculating terrain parameters at GPR data locations']);
    allvalues = DEM.(names{n});
    vectorvalues = allvalues(:);
    nonanvalues = vectorvalues(index,1);
    %aggregate.(names{n}) = nonanvalues(aggregate.index);
    aggregate_test.(names{n}) = nonanvalues(aggregate_test.index);
    aggregate_train.(names{n}) = nonanvalues(aggregate_train.index);
    
end

aggregate_train=aggregate_train(~any(ismissing(aggregate_train),2),:);
aggregate_test=aggregate_test(~any(ismissing(aggregate_test),2),:);

%% test plots
clear ind_train ind_test SWE_train_c SWE_train_r SWE_test_c SWE_test_r

ind_train = find(aggregate_train.slope<=maxslope);
SWE_train_c = ((aggregate_train.x(ind_train) - DEM.ref.XWorldLimits(1))./DEM.cellsize) + 1/2; %interp
SWE_train_r = ((DEM.ref.YWorldLimits(2) - aggregate_train.y(ind_train))./DEM.cellsize) + 1/2;

ind_test = find(aggregate_test.slope<=maxslope);
SWE_test_c = ((aggregate_test.x(ind_test) - DEM.ref.XWorldLimits(1))./DEM.cellsize) + 1/2; %interp
SWE_test_r = ((DEM.ref.YWorldLimits(2) - aggregate_test.y(ind_test))./DEM.cellsize) + 1/2;

if strcmp(plots, 'on')==1
for n= 1:length(names)
figure ();
colormap(jet)
scatter(SWE_train_c,SWE_train_r,3,GPR.(names{n}))
axis ij;
axis image;
colorbar;
caxis([0,swemax]);
xlabel('east');
ylabel('north');
text(700, 900, labels{n});
end
end

%% build tables for regression
clear xySWE NormalxySWE SWEmu SWEsigma
xySWE_train = aggregate_train.SWE ./ cosd(aggregate_train.slope); %corrects for slope
[NormalxySWE_train SWEmu_train SWEsigma_train] = zscore(xySWE_train); % standardized

TerrainTable_train = table(NormalxySWE_train(ind_train));
TerrainTable_train.Properties.VariableNames{'Var1'} = 'SWE';
for n= 1:length(names)
    [add mu sigma] = zscore(aggregate_train.(names{n})); %standardized
    TerrainTable_train.(names{n}) = add(ind_train);
    muV_train.(names{n}) = mu;
    sigmaV_train.(names{n})= sigma;
end

%Test dataset is normalized by the values (mean, std) from the training dataset

xySWE_test = aggregate_test.SWE ./ cosd(aggregate_test.slope); %corrects for slope

NormalxySWE_test = (xySWE_test-SWEmu_train)./SWEsigma_train;

TerrainTable_test = table(NormalxySWE_test(ind_test));
TerrainTable_test.Properties.VariableNames{'Var1'} = 'SWE';
for n= 1:length(names)
    TerrainTable_test.(names{n}) = (aggregate_test.(names{n})-muV_train.(names{n}))./sigmaV_train.(names{n});
end

%% Regression
%'Criterion' allows you to choose what type of cuttoff gets used for including params. choices are 'Variance', 'SSE','AIC','BIC', 'RSquared',or 'AdjRSquared'.
%'PEnter' allows you to choose the cutoff for the designated criterion
%'Upper' allows you to limit model complexity, where 'linear' means you will only include linear terms rather than interaction terms, etc.
mdl_{i} = stepwiselm(TerrainTable_train,'ResponseVar','SWE','Criterion','AIC','PEnter',0,'PRemove',0.01,'upper','linear','verbose',0,'PredictorVars',{'Z','Sb','northness','eastness'});

%% Loop to find Curvature
if strcmp(findOptimumCurvature, 'yes')==1 
 wb = waitbar(1/100, 'OptimizeCurvature Loop');  
 [CurvatureValues] = OptimizeCurvature(DEM.Z, DEM.cellsize, TerrainTable, gprIndicies, ind, wb);  
end

%% Loop to find Sx
if strcmp(findOptimumSx, 'yes')==1 
  wb = waitbar(1/100, 'OptimizeSx Loop');  
 [SxValues] = OptimizeSx(DEM.Z, DEM.cellsize, TerrainTable, gprIndicies, ind, wb);   
end

%%
%figure showing what survey lines were included and not included in the
%regression results

%calculate dimensions for figure
if itr==5
    xx=1;
    yy=5;
elseif itr==10
    xx=2;
    yy=5;
elseif itr==50
    xx=5;
    yy=10;
elseif itr==100
    xx=10;
    yy=10;

end

subplot(xx,yy,i)
hold on, box on, grid on,
scatter(SWE_train_c,SWE_train_r,5,'k'); 
hold on
scatter(SWE_test_c,SWE_test_r,5,'r'); 
title(['Subset ' num2str(i)]);
%caxis([0 6]);
xlim([200 850]);
ylim([200 1050]);
set(gca,'Ydir','reverse')

%% Overlay residuals on hillshade
do=2;
if do==1

%this will produce a figure for each iteration - potentially a lot!
figure()

ax1 = axes;
imagesc(DEM.shade, 'alphadata', ~isnan(DEM.Z));hold on
axis([0 DEM.ref.RasterExtentInWorldX 0 DEM.ref.RasterExtentInWorldY]);
axis ij;
axis image;

ax2 = axes;
linkaxes([ax1,ax2])
scatter(SWE_train_c,SWE_train_r,5,mdl.Residuals.Raw); 
hold on
scatter(SWE_test_c,SWE_test_r,'.k');

axis ij;
axis image;
caxis([-0.5 0.5]);
ax2.Visible = 'off';
ax2.XTick = [];
ax2.YTick = [];
xlim([200 850]);
ylim([200 1050]);
set(gca,'Fontsize',24); 
colormap(ax1,'gray')
colormap(ax2,'parula')
set([ax1,ax2],'Position',[.10 .12 .7 .78]);
cb2 = colorbar(ax2,'Position',[.88 .12 .05 .78]);
title(ax1,'Model Residuals');
xlabel(ax1,'east [m]');
ylabel(ax1,'north [m]');
end
%%
clear ypred
%calculate SWE for Terrain_test dataset
ypred = predict(mdl_{i},TerrainTable_test);

%bring back to real values
ypred = ypred .*SWEsigma_train + SWEmu_train;

%calculate RMSE
%corrected SWE estimate
RMSE_{i} = sqrt(mean((ypred - xySWE_test).^2));  % Root Mean Squared Error
%calculate R2 and RMSE
R2_{i} = corr(xySWE_test,ypred)^2;

%make new terrain table that covers the entire glacier area but has been
%normalized by the training dataset since that is what the regression tree
%is based on

Z=((reshape(DEM.Z,[],1))-muV_train.Z)./sigmaV_train.Z;
Sb=((reshape(DEM.Sb,[],1))-muV_train.Sb)./sigmaV_train.Sb;
%Sx=((reshape(DEM.Sx,[],1))-muV_train.Sx)./sigmaV_train.Sx;
curvature=((reshape(DEM.curvature,[],1))-muV_train.curvature)./sigmaV_train.curvature;
northness=((reshape(DEM.northness,[],1))-muV_train.northness)./sigmaV_train.northness;
eastness=((reshape(DEM.eastness,[],1))-muV_train.eastness)./sigmaV_train.eastness;
slope=((reshape(DEM.slope,[],1))-muV_train.slope)./sigmaV_train.slope;
DEM_TerrainTable=table(Z,Sb,curvature,northness,eastness,slope);
yHatStandardized = predict(mdl_{i},DEM_TerrainTable); %use the model to predict the full field
yHatLinear = yHatStandardized .* SWEsigma_train + SWEmu_train; %returns it from standardized to actual
XSWE.all_{i} = reshape(yHatLinear,[1316,1095]); 
clear yHatStandardized yHatLinear

%%
%calculating glacier wide averages for the mvr

XSWE.on_{i} = NaN(size(XSWE.all_{i}));
XSWE.on_{i}(GL.aligned==1) = XSWE.all_{i}(GL.aligned==1);
XSWE.on_{i}(XSWE.on_{i}<0)=0;

%GL.SWE_mean = nansum(nansum(XSWE.on))/GL.area;
GL.SWE_{i} = nansum(nansum(XSWE.on_{i}))/GL.area;

clear p Terrain_train Terrain_test ypred SWE3 temp SWE_test_est
%%
%this section runs the random forest regression tree after the
%elevation trend has been removed

%find linear fit between SWE and elevation for the training dataset
temp_SWE=xySWE_train;
temp_elev=aggregate_train.Z;
zt = polyfit(temp_elev,temp_SWE,1);

%this detrends the training dataset using the training dataset elevation trend
temp_SWE_train2=temp_SWE-(temp_elev*zt(1,1)+zt(1,2));

%standardize the detrended dataset
[temp_SWE_train3, rt_mu, rt_sig]=zscore(temp_SWE_train2);

%build new training terrain table without Z
TerrainTable_train_detrend=table(temp_SWE_train3,TerrainTable_train.Sb, TerrainTable_train.curvature,TerrainTable_train.northness, TerrainTable_train.eastness, TerrainTable_train.slope);
TerrainTable_train_detrend.Properties.VariableNames = {'SWE' 'Sb' 'curvature' 'northness' 'eastness' 'slope'};

%this de-elevation trends the test dataset
temp_SWE_test2=aggregate_test.SWE-(aggregate_test.Z*zt(1,1)+zt(1,2));

%standardize the detrended dataset by the training set mean/st dev
temp_SWE_test3=(temp_SWE_test2-rt_mu)./rt_sig;

%build new training terrain table without Z
TerrainTable_test_detrend=table(temp_SWE_test3,TerrainTable_test.Sb,TerrainTable_test.curvature,TerrainTable_test.northness, TerrainTable_test.eastness, TerrainTable_test.slope);
TerrainTable_test_detrend.Properties.VariableNames = {'SWE' 'Sb' 'curvature' 'northness' 'eastness' 'slope'};

%%
%make template
t = templateTree('NumPredictorsToSample','all',...
    'PredictorSelection','interaction-curvature','Surrogate','on');
rng(1); % For reproducibility

%find optimal regression tree
Mdl_detrend_{i} = fitrensemble(TerrainTable_train_detrend,'SWE','Method','bag','NumLearningCycles',100,'Learners',t);

%Estimate predictor importance measures by permuting out-of-bag observations. Compare the estimates using a bar graph.
RT_param_sig_{i} = oobPermutedPredictorImportance(Mdl_detrend_{i});

% now apply this model to the rest of the observational data (test)
yHat=predict(Mdl_detrend_{i},TerrainTable_test_detrend);

%bring back into real space - minus elevation trend
yHat = yHat .*rt_sig + rt_mu;

%now add back in elevation trend
yHat=yHat + (aggregate_test.Z*zt(1,1) + zt(1,2));

%calculate R2 and RMSE
%this is the cosd modified SWE estimate
R2_rt_{i} = corr(xySWE_test,yHat)^2;
RMSE_rt_{i}=sqrt(mean((xySWE_test-yHat).^2));

%%
%extrapolate regression tree results to entire glacier area

%predict SWE at location and then bring back into real value
    yHat2=predict(Mdl_detrend_{i},DEM_TerrainTable);
    yHat2 = yHat2 .*rt_sig + rt_mu;

%reshape DEM.Z
Z2=reshape(DEM.Z,[],1);

%now add back in elevation trend
yHat2=yHat2 + (Z2*zt(1,1) + zt(1,2));

%reshape SWE vector into map view
SWE_rt_{i}=reshape(yHat2,[1316,1095]);

%calculate mean
SWE_DEM_mean_{i}.on = NaN(size(SWE_rt_{i}));
SWE_DEM_mean_{i}.on(GL.aligned==1) = SWE_rt_{i}(GL.aligned==1);
SWE_DEM_mean_{i}.on(SWE_DEM_mean_{i}.on<0)=0;
SWE_rt_mean_{i} = nansum(nansum(SWE_DEM_mean_{i}.on))/GL.area;


end

%%
%compile and boxplot glacier wide averages from the multi-variable
%regression and decision tree
for i=1:itr
    GL.SWE(i,1)=GL.SWE_{1,i};
end

%create boxplot of the decision tree mean SWEs
%add the dtree results in the MVR file on row 2
GL.SWE(:,2)=cell2mat(SWE_rt_mean_);

figure,
hold on, box on, grid on
bh=boxplot(GL.SWE(:,:));
set(gca,'Fontsize',30); 
ylabel('Glacier-wide Mean SWE [m]');
set(bh,'linewidth',4);
ylim([mean(median(GL.SWE(:,:)))-.75 mean(median(GL.SWE(:,:)))+.75])
hold on
set(gca,'XTickLabel',{'MVR','Decision Tree'});

%%
%calculate spatial pattern of median value from both the MVR and Regression
%Tree

%create multi-dimensional array and then take median value and plot
for i=1:itr
GL_array(:,:,i)=XSWE.on_{1,i};
end

%take median value of MVR
GL_array_med=median(GL_array(:,:,:),3);
GL_array_stdev=std(GL_array(:,:,:),0,3);

%now do for regression tree 
for i=1:itr
GL_array_rt(:,:,i)=SWE_DEM_mean_{1,i}.on;
end

%take median value
GL_array_med_rt=median(GL_array_rt(:,:,:),3);
GL_array_stdev_rt=std(GL_array_rt(:,:,:),0,3);

%plot median distributed value
figure
subplot(2,3,1)
imagesc(GL_array_med, 'alphadata', GL.aligned);hold on
axis([0 DEM.ref.RasterExtentInWorldX 0 DEM.ref.RasterExtentInWorldY]);
axis ij;
axis image;
h1 = gca;
caxis([0,6]);
%colormap(h1,jet); 
hold on, grid on

c=colorbar('eastoutside');
xlim([200 850]);
ylim([200 1050]);
set(gca,'Fontsize',16);
xlabel('east');
ylabel('north');
ylabel(c,'SWE [m]');
title('MV Regression');

subplot(2,3,2)
imagesc(GL_array_stdev, 'alphadata', GL.aligned);hold on
axis([0 DEM.ref.RasterExtentInWorldX 0 DEM.ref.RasterExtentInWorldY]);
axis ij;
axis image;
h1 = gca;
caxis([0,0.5]);
%colormap(h1,jet); 
hold on, grid on

c=colorbar('eastoutside');
xlim([200 850]);
ylim([200 1050]);
set(gca,'Fontsize',16);
xlabel('east');
ylabel('north');
ylabel(c,'SWE [m]');
title('St. Dev MVR');

subplot(2,3,4)
imagesc(GL_array_med_rt, 'alphadata', GL.aligned);hold on
axis([0 DEM.ref.RasterExtentInWorldX 0 DEM.ref.RasterExtentInWorldY]);
axis ij;
axis image;
h1 = gca;
caxis([0,6]);
%colormap(h1,jet); 
hold on, grid on

c=colorbar('eastoutside');
xlim([200 850]);
ylim([200 1050]);
set(gca,'Fontsize',16);
xlabel('east');
ylabel('north');
ylabel(c,'SWE [m]');
title('Regression Tree');

subplot(2,3,5)
imagesc(GL_array_stdev_rt, 'alphadata', GL.aligned);hold on
axis([0 DEM.ref.RasterExtentInWorldX 0 DEM.ref.RasterExtentInWorldY]);
axis ij;
axis image;
h1 = gca;
caxis([0,0.5]);
%colormap(h1,jet); 
hold on, grid on

c=colorbar('eastoutside');
xlim([200 850]);
ylim([200 1050]);
set(gca,'Fontsize',16);
xlabel('east');
ylabel('north');
ylabel(c,'SWE [m]');
title('St. Dev R. Tree');

subplot(1,3,3)
%calculate difference
diff=GL_array_med-GL_array_med_rt;

imagesc(diff, 'alphadata', GL.aligned);hold on
axis([0 DEM.ref.RasterExtentInWorldX 0 DEM.ref.RasterExtentInWorldY]);
axis ij;
axis image;
h1 = gca;
caxis([-1,1]);
%colormap(h1,jet); 
hold on, grid on

c=colorbar('eastoutside');
xlim([200 850]);
ylim([200 1050]);
set(gca,'Fontsize',16);
xlabel('east');
ylabel('north');
ylabel(c,'SWE [m]');
title('Difference');



%% plot GPR observations overlaid on hillshade
figure
%set(gca,'Fontsize',24); 
ax1 = axes;
imagesc(DEM.shade, 'alphadata', ~isnan(DEM.Z));hold on
axis([0 DEM.ref.RasterExtentInWorldX 0 DEM.ref.RasterExtentInWorldY]);
axis ij;
axis image;

ax2 = axes;
linkaxes([ax1,ax2])
scatter(SWE_entire_c,SWE_entire_r,10,aggregate.SWE); hold on
axis ij;
axis image;
%caxis([-swemax,swemax]);
caxis([0,6]);
ax2.Visible = 'off';
ax2.XTick = [];
ax2.YTick = [];
xlim([200 850]);
ylim([200 1050]);

colormap(ax1,'gray')
cmap=colormap(ax2,('parula'));
cmap=flipud(cmap);
colormap(ax2,cmap);

set([ax1,ax2],'Position',[.1 .10 .7 .75]);
cb2 = colorbar(ax2,'Position',[.85 .10 .05 .75]);
set(gca,'Fontsize',24); 

%CHANGE DATES HERE
title(ax1,'SWE Observations 2017');
xlabel(ax1,'east [m]');
ylabel(ax1,'north [m]');

clear ax1 ax2

toc
