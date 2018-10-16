function [ SxValues ] = OptimizeSx(DEMZ, DEMcellsize, ObTable, gprIndicies, ind, wb)
%OptimizeSx simply loops through a range of distance and direction values.
%It then replaces the old Sx value in the table and runs a new stepwise
%regression.
%It outputs a table that gives the direction, distance, Sx coefficient
%estimate, and the total model R2. NaN values in the R2 and Coefficient 
%estimate mean that Sx did not appear in the model given that distance and 
%direction. Larger coefficients indicate that the Sx value is explaining 
%more of the variance, but that can be confusing as the magnitude of Sx can 
%also vary. So the model R2 value is the most important measure of how 
%useful Sx is. 


SxValues = [];

for n = 1:5;
    direction = 2*n + 34;
    waitbar(n/5, wb, 'OptimizeSx is slower - come back tomorrow');
    index = ~isnan(DEMZ(:));
    for m = 1:5;
        distance = 200*m + 1400;
         allvalues = Sx(DEMZ, DEMcellsize, distance, direction);
         vectorvalues = allvalues(:);
        nonanvalues = vectorvalues(index,1);
        rvalues= nonanvalues(gprIndicies);
        [ObTable.Sx] = zscore(rvalues(ind));
        P1 = find(strcmp(ObTable.Properties.VariableNames,'Z'));
        P2 = find(strcmp(ObTable.Properties.VariableNames,'Sx'));
        mdl = stepwiselm(ObTable,'ResponseVar','SWE','upper', 'linear'); %'PredictorVars',[P1,P2],
        SxLoc = strcmp(mdl.PredictorNames,'Sx');
        if sum(SxLoc) > 0;
            SxValues = [SxValues;[direction, distance, mdl.Coefficients.Estimate(SxLoc), mdl.Rsquared.Ordinary]];
        else
            SxValues = [SxValues;[direction, distance, NaN, NaN]];
        end

    end
end

figure()
scatter(SxValues(:,1),SxValues(:,2),10,SxValues(:,3))
xlabel('direction')
ylabel('distance')
title('Sx Coeficient Estimate')

figure()
scatter(SxValues(:,1),SxValues(:,2),10,SxValues(:,4))
xlabel('direction')
ylabel('distance')
title('Model R2')

close(wb)

end

