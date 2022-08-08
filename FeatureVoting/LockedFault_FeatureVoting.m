function [votingMatrix] = LockedFault_FeatureVoting(Set0,Set1,posE,posN,AZs,VFs,LDs, V_th)
%LockedFault_FeatureVoting computes the feature voting for a locked fault model, ona single postion
%for the fault fixing
%   Set0, Set1 are the measurment data sets
%   posE, posN - is the fixing location for the fault
%   AZs, VFs, LDs - are the values to be voted in the process for Azimuths[radians], Velocities[mm] and Lock-Depths[m]
%   V_th is the threshold for velocty acceptance criteria
%
% Omer Bar 2020 Jan, version 1.0
%% 
a  = zeros(numel(AZs),numel(VFs),numel(LDs));
[ai,vi,li]=ind2sub(size(votingMatrix),1:numel(votingMatrix));

% vectors speed
dn = (Set1.vctrs(:,1)-Set0.vctrs(:,1));
de = (Set1.vctrs(:,2)-Set0.vctrs(:,2)); 
vs = sqrt((dn.*dn) +...      
          (de.*de)      );  
% vs = sqrt((Set1.vctrs(:,1)-Set0.vctrs(:,1)).*(Set1.vctrs(:,1)-Set0.vctrs(:,1)) +...      %dn2
%     (Set1.vctrs(:,2)-Set0.vctrs(:,2)).*(Set1.vctrs(:,2)-Set0.vctrs(:,2))      );   %de2
vs = vs * 1000; %mm per year
i1 = (abs(dn) > de) .* (dn < 0); % where dn < 0
i2 = (abs(de) > dn) .* (de < 0); % where de < 0

i = find((i1+i2) >0);  % where (de < 0) or (dn < 0)
vs(i) = -vs(i);

for k = 1:size(vs)% moving on each vector_size_change
    pAi = find(string(Set0.VectorsAndVCVs{k,1})==Set0.Points{:,1}); % Vector's FromPoint
    pBi = find(string(Set0.VectorsAndVCVs{k,2})==Set0.Points{:,1}); % Vector's   ToPoint
    % Distance of FromPoint from simulated fault
    lpA = round(DistLinePoint(  posE-1e6.*sin(AZs(ai))      ,posN-1e6.*cos(AZs(ai)),...
                                posE+1e6.*sin(AZs(ai))      ,posN+1e6.*cos(AZs(ai)),...
                        ones(1,numel(ai)).*Set0.crds(pAi,2) ,ones(1,numel(ai)).*Set0.crds(pAi,1)));
    % Left or Right from the fault
    l_r_A =PointSideFromLine(   posE-1e6.*sin(AZs(ai)),posN-1e6.*cos(AZs(ai)),...
                                posE+1e6.*sin(AZs(ai)),posN+1e6.*cos(AZs(ai)),...
        ones(1,numel(ai)).*Set0.crds(pAi,2),ones(1,numel(ai)).*Set0.crds(pAi,1));
    % -- %
    % Distance of ToPoint from simulated fault
    lpB = round(DistLinePoint(  posE-1e6.*sin(AZs(ai)),posN-1e6.*cos(AZs(ai)),...
                                posE+1e6.*sin(AZs(ai)),posN+1e6.*cos(AZs(ai)),...
        ones(1,numel(ai),1).*Set0.crds(pBi,2),ones(1,numel(ai)).*Set0.crds(pBi,1)));
    % Left or Right from the fault
    l_r_B =PointSideFromLine(   posE-1e6.*sin(AZs(ai)),posN-1e6.*cos(AZs(ai)),...
                                posE+1e6.*sin(AZs(ai)),posN+1e6.*cos(AZs(ai)),...
        ones(1,numel(ai)).*Set0.crds(pBi,2),ones(1,numel(ai)).*Set0.crds(pBi,1));
    
    % computing vector change rate by senarios (indexes by pixels of H3)
    l_r_A = l_r_A'; l_r_B = l_r_B';
    tempv = (VFs(vi)./pi) .* atan(lpB./LDs(li)).*l_r_B  -  (VFs(vi)./pi) .* atan(lpA./LDs(li)).*l_r_A ;  %             tempv = -((VFs(vi)./1000)./pi) .* atan(lpA.*l_r_A./LDs(li)) + ((VFs(vi)./1000)./pi) .* atan(lpB.*l_r_B./LDs(li));
    % tempv = (VFs(vi)./pi) .* atan(   (LDs(li).*(lpB - lpA)) ./ (LDs(li)+(lpB .* lpA))    );
    
    % selecting matching senarios to measured data (to a treshold)
    i = abs(tempv-vs(k)) < V_th;
    
    % voting
    votingMatrix(i) = votingMatrix(i)+1;
end % for k = 1:size(t0Set.crds,1)

votingMatrix = votingMatrix / numel(vs); % precentage only!
end