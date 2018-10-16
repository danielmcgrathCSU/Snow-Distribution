function [ CurvatureValues ] = OptimizeCurvature(DEMZ, DEMcellsize, ObTable, gprIndicies, ind, wb)
%OptimizeCurvature simply loops through a range of radius values.
%It then replaces the old Curvature value in the table and runs a new stepwise
%regression.
%It outputs a table that gives the radius, Curvature coefficient
%estimate, and the total model R2. NaN values in the R2 and Coefficient 
%estimate mean that Curvature did not appear in the model given that radius. 
%Larger coefficients indicate that the Curvature value is explaining 
%more of the variance, but that can be confusing as the magnitude of Curvature can 
%also vary. So the model R2 value is the most important measure of how 
%useful Curvature is. 

CurvatureValues = [];

for n = 1:50;
    radius = 10*n;
    waitbar(n/50, wb, 'OptimizeCurvature is slow - get some coffee');
    index = ~isnan(DEMZ(:));
         allvalues = Curvature(DEMZ, DEMcellsize, radius, 'method', 'xy');
         vectorvalues = allvalues(:);
        nonanvalues = vectorvalues(index,1);
        rvalues= nonanvalues(gprIndicies);
        [ObTable.curvature] = zscore(rvalues(ind));
        P1 = find(strcmp(ObTable.Properties.VariableNames,'Z'));
        P2 = find(strcmp(ObTable.Properties.VariableNames,'curvature'));
        mdl = stepwiselm(ObTable,'ResponseVar','SWE','PredictorVars',[P1,P2],'upper', 'linear'); %
        curvatureLoc = strcmp(mdl.PredictorNames,'curvature');
        if sum(curvatureLoc) > 0;
            CurvatureValues = [CurvatureValues;[radius, mdl.Coefficients.Estimate(curvatureLoc), mdl.Rsquared.Ordinary]];
        else
            CurvatureValues = [CurvatureValues;[radius, NaN, NaN]];
        end
end

figure()
scatter(CurvatureValues(:,1),CurvatureValues(:,3),3,CurvatureValues(:,2))
xlabel('radius');
ylabel('R2');
title('Curvature Radius Estimate');

close(wb)

end

