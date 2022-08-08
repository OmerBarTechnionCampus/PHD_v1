function [newSets] = moveSet3D_slipFault(dt,vFhz,vFv,azF,posF,Set0,noiseHz,noiseV,nTimes)
%%  2019.03.05  - tested ok +- 1mm

% updated 2019.03.05
newSets = {Set0}; newSets{nTimes} = NaN;

% minE = floor(min(Set0.crds(:,2))) - 3000;% minN = floor(min(Set0.crds(:,1))) - 3000;
% maxE = floor(max(Set0.crds(:,2))) + 3000;% maxN = floor(max(Set0.crds(:,1))) + 3000;
% m = tan(azF);
% e1 = minE; n1 = m*(posF(1) - e1) + posF(2);
% e2 = maxE; n2 = m*(posF(1) - e2) + posF(2);
e1 = posF(1) -10000.*sin(azF);   e2 = posF(1) +10000.*sin(azF);
n1 = posF(2) -10000.*cos(azF);   n2 = posF(2) +10000.*cos(azF);

for t = 1 : (nTimes - 1)
    set_new = newSets{t};                 
    set_new.Time = set_new.Time + dt;
    
    for nms = 1:size(Set0.Points,1)
        for i = 1:size(set_new.Points,1)
            if (strcmpi(set_new.Points{i,1}{1},Set0.Points{nms,1}{1}) == true)
                e = set_new.crds(i,2); n = set_new.crds(i,1);
                l_r = PointSideFromLine(e1,n1,e2,n2,e,n);  % left or right from line
%                 l  = DistLinePoint(e1,n1,e2,n2,e,n);

                v = vFhz / 2;
                de = v*dt*sin(azF)       * l_r;
                dn = v*dt*cos(azF)       * l_r;
                du = vFv * dt            * l_r;
            
                %updating crds
                set_new.crds(i,2) =  set_new.crds(i,2) + de;
                set_new.crds(i,1) =  set_new.crds(i,1) + dn;
                set_new.crds(i,3) =  set_new.crds(i,3) + du;
                
                % rectifying vectors 
                % when from Point is moving
                p = string(set_new.VectorsAndVCVs{:,1})==Set0.Points{nms,1}{1};
                set_new.VectorsAndVCVs{ (p == 1) ,3:5} = set_new.VectorsAndVCVs{ (p == 1) ,3:5} - [dn,de,du];
                
                % when  to  Point is moving
                p = string(set_new.VectorsAndVCVs{:,2})==Set0.Points{nms,1}{1};
                set_new.VectorsAndVCVs{ (p == 1) ,3:5} = set_new.VectorsAndVCVs{ (p == 1) ,3:5} + [dn,de,du];
                continue;
            end %if
        end %for i
    end %for nms
    
    % updaing the entire DataSet
    set_new.Points{:,2:4} = set_new.crds;
    set_new.vctrs         = set_new.VectorsAndVCVs{:,3:5};
    % setting up error estimates - adding zero-avarage noise
    set_new.vcvs = rand(size(set_new.vcvs)).*erf(mean(mean(abs(set_new.vcvs))));
    
    % adding noise
    rng('shuffle');
    set_new.vcvs         = set_new.vcvs         + mean([noiseHz,noiseV]) .* rand(size(set_new.vcvs));                   % ------ add the v-Noise to the rnadomization
    set_new.vctrs(:,1:2) = set_new.vctrs(:,1:2) +      noiseHz           .* rand(size(set_new.vctrs,1),2);
    set_new.vctrs(:,3)   = set_new.vctrs(:, 3 ) +      noiseV            .* rand(size(set_new.vctrs,1),1);
    set_new.VectorsAndVCVs{:,3:end} = [set_new.vctrs,set_new.vcvs];
    newSets{t+1} = set_new;
end % for Times

end% end function