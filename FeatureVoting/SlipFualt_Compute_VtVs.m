function [v_mat] = SlipFualt_Compute_VtVs(Set0,Set1,posE,posN,AZs,VFs)
%Compute_VtVs is a function to describe the dis-simmilarity between a
%   parametrized model of a LockedFualt to the measuremtts taken
%   Set0, Set1 are the measurment data sets
%   posE, posN - is the fixing location for the fault
%   AZs, VFs, LDs - are the values to be voted in the process for Azimuths[radians], Velocities[mm] and Lock-Depths[m]
%   V_th is the threshold for velocty acceptance criteria
%
%Omer Bar,2020 Jan., version 1.0
%
%% preparing the correct form of data
% if (size(posE,1) == 1); posE = posE';end;
% if (size(posN,1) == 1); posN = posN';end;
% if (size(AZs,1) == 1);  AZs = AZs';  end;
% if (size(VFs,1) == 1);  VFs = VFs';  end;
% if (size(LDs,1) == 1);  LDs = LDs';  end;

%%
dt = Set1.Time - Set0.Time;
% vectors speed
dn = (Set1.vctrs(:,1)-Set0.vctrs(:,1)); %dn
de = (Set1.vctrs(:,2)-Set0.vctrs(:,2)); %de
vs = sqrt((dn.*dn) +...      
          (de.*de)      );  
      vs = vs *1000 ./ dt; %mm per year           ./ dt for multiepoch
i1 = (abs(dn) > de) .* (dn < 0); % where dn < 0
i2 = (abs(de) > dn) .* (de < 0); % where de < 0
i = find((i1+i2) > 0);   % where (de < 0) or (dn < 0)
vs(i) = -vs(i);


v_mat = zeros(numel(AZs),size(Set0.VectorsAndVCVs,1));
for k = 1:numel(vs) % runnig on vectors
    % DistLinePoint(e1,n1,e2,n2,ep,np)
    pAi = find(string(Set0.VectorsAndVCVs{k,1})==Set0.Points{:,1}); % vectors FromPoint
    pBi = find(string(Set0.VectorsAndVCVs{k,2})==Set0.Points{:,1}); % vectors   ToPoint
    
%     % Dist from fualt
%     lpA = round(DistLinePoint(  posE-1e6.*sin(AZs),posN-1e6.*cos(AZs),...
%                                 posE+1e6.*sin(AZs),posN+1e6.*cos(AZs),...
%         ones(1,numel(AZs)).*Set0.crds(pAi,2),ones(1,numel(AZs)).*Set0.crds(pAi,1)));
    % left or right from fualt
    l_r_A =PointSideFromLine(   posE-1e6.*sin(AZs),posN-1e6.*cos(AZs),...
                                posE+1e6.*sin(AZs),posN+1e6.*cos(AZs),...
        ones(1,numel(AZs)).*Set0.crds(pAi,2),ones(1,numel(AZs)).*Set0.crds(pAi,1));
%     % Dist from fualt
%     lpB = round(DistLinePoint(  posE-1e6.*sin(AZs),posN-1e6.*cos(AZs),...
%                                 posE+1e6.*sin(AZs),posN+1e6.*cos(AZs),...
%         ones(1,numel(AZs),1).*Set0.crds(pBi,2),ones(1,numel(AZs)).*Set0.crds(pBi,1)));
    % left or right from fualt
    l_r_B =PointSideFromLine(   posE-1e6.*sin(AZs),posN-1e6.*cos(AZs),...
                                posE+1e6.*sin(AZs),posN+1e6.*cos(AZs),...
        ones(1,numel(AZs)).*Set0.crds(pBi,2),ones(1,numel(AZs)).*Set0.crds(pBi,1));
    % % current vector change in size and direction due to fualt scenarios
    % tempv = (VFs(vi)./pi) .* atan(   (LDs(li).*(lpB - lpA)) ./ (LDs(li)+(lpB .* lpA))    );
    
    % deviation of the simultion from the measured change in the vector
%      l_r_A = l_r_A'; l_r_B = l_r_B';
     
%     if (size(lpB,1) == 1); lpB = lpB';end;
%     if (size(lpA,1) == 1); lpA = lpA';end;
    if (size(VFs,1) == 1);  VFs = VFs';  end;
%     if (size(LDs,1) == 1);  LDs = LDs';  end;
    if (size(l_r_A,1) == 1);  l_r_A = l_r_A';  end;
    if (size(l_r_B,1) == 1);  l_r_B = l_r_B';  end;

    
    v_mat(:,k) =( (VFs)./2.*l_r_B - (VFs)./2.*l_r_A ) -  vs(k); % .*tan(deltaAz)
%     %%%%% the angle beween the 2 in needed
%     dotP = 0;
end %for
end %function
