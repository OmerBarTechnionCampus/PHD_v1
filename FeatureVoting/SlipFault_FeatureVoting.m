function [votingMatrix5] = SlipFault_FeatureVoting(Set0,Set1,posE,posN,AZs,VFs,V_th)
%LockedFault_FeatureVoting computes the feature voting for a locked fault model, ona single postion
%for the fault fixing
%   Set0, Set1 are the measurment data sets
%   posE, posN - is the fixing location for the fault
%   AZs, VFs - are the values to be voted in the process for Azimuths[radians] and Velocities[mm] 
%   V_th is the threshold for velocty acceptance criteria
%
% Omer Bar 2020 Jan, version 1.0
%% 
votingMatrix5  = zeros(numel(AZs),numel(VFs),numel(posE),numel(posN));
[ai,vi,pei,pni]=ind2sub(size(votingMatrix5),1:numel(votingMatrix5));

% vectors speed
dn = (Set1.vctrs(:,1)-Set0.vctrs(:,1));
de = (Set1.vctrs(:,2)-Set0.vctrs(:,2)); 
vs = sqrt((dn.*dn) +...      
          (de.*de)      ) ./ (Set1.Time - Set0.Time);   %   ============== sqrt(de^2+dn^2)./ dt for multiepoch
az = round(Azimuth(Set1.vctrs(:,2),Set1.vctrs(:,1)),6);
% vs = sqrt((Set1.vctrs(:,1)-Set0.vctrs(:,1)).*(Set1.vctrs(:,1)-Set0.vctrs(:,1)) +...      %dn2
%     (Set1.vctrs(:,2)-Set0.vctrs(:,2)).*(Set1.vctrs(:,2)-Set0.vctrs(:,2))      );   %de2
vs = vs * 1000; %mm per year
i1 = (abs(dn) > de) .* (dn < 0); % where dn < 0
i2 = (abs(de) > dn) .* (de < 0); % where de < 0

i = find((i1+i2) >0);  % where (de < 0) or (dn < 0)
vs(i) = -vs(i);  % minus sign for shortening vectors

for k = 1:size(vs)% moving on each vector_size_change
    dts = datetime('now');
    pAi = find(string(Set0.VectorsAndVCVs{k,1})==Set0.Points{:,1}); % Vector's FromPoint
    pBi = find(string(Set0.VectorsAndVCVs{k,2})==Set0.Points{:,1}); % Vector's   ToPoint
%     % Distance of FromPoint from simulated fault
%     %round
%     lpA = (DistLinePoint(  posE(pei)-1e6.*sin(AZs(ai))      ,posN(pni)-1e6.*cos(AZs(ai)),...
%                            posE(pei)+1e6.*sin(AZs(ai))      ,posN(pni)+1e6.*cos(AZs(ai)),...
%                         ones(1,numel(ai)).*Set0.crds(pAi,2) ,ones(1,numel(ai)).*Set0.crds(pAi,1)));
    % Left or Right from the fault
    l_r_A =PointSideFromLine(   posE(pei)-1e6.*sin(AZs(ai)),posN(pni)-1e6.*cos(AZs(ai)),...
                                posE(pei)+1e6.*sin(AZs(ai)),posN(pni)+1e6.*cos(AZs(ai)),...
        ones(1,numel(ai)).*Set0.crds(pAi,2),ones(1,numel(ai)).*Set0.crds(pAi,1));
    % -- %
%     % Distance of ToPoint from simulated fault
%     %round
%     lpB = (DistLinePoint(  posE(pei)-1e6.*sin(AZs(ai)),posN(pni)-1e6.*cos(AZs(ai)),...
%                                 posE(pei)+1e6.*sin(AZs(ai)),posN(pni)+1e6.*cos(AZs(ai)),...
%         ones(1,numel(ai),1).*Set0.crds(pBi,2),ones(1,numel(ai)).*Set0.crds(pBi,1)));
    % Left or Right from the fault
    l_r_B =PointSideFromLine(   posE(pei)-1e6.*sin(AZs(ai)),posN(pni)-1e6.*cos(AZs(ai)),...
                                posE(pei)+1e6.*sin(AZs(ai)),posN(pni)+1e6.*cos(AZs(ai)),...
        ones(1,numel(ai)).*Set0.crds(pBi,2),ones(1,numel(ai)).*Set0.crds(pBi,1));
    
    % computing vector change rate by senarios (indexes by pixels of H3)
    l_r_A = l_r_A'; l_r_B = l_r_B';
    tempv = VFs(vi)./2.*l_r_B  -  VFs(vi)./2.*l_r_A ;  
    
    % selecting matching senarios to measured data (to a treshold)
    i = abs(tempv-vs(k)) < V_th;

    % voting    
    votingMatrix5(i) = votingMatrix5(i)+1;
    if (k==1)
        fprintf('ended vector numer %d (time = %f seconds)\n',k,seconds(datetime('now')-dts));
    end
end % for k = 1:size(t0Set.crds,1)

% votingMatrix5 = votingMatrix5 ./ numel(vs) .*100 ; % precentage only!
end
