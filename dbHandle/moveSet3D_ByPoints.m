function [newSets] = moveSet3D_ByPoints(PointNames,dt,vFhz,vFv,azF,Set0,noiseHz,noiseV,nTimes)
% updated 2019.03.01
if (numel(azF) == 1)
    vFhz    = repmat(vFhz(1)   ,size(PointNames));
    vFv     = repmat(vFv(1)    ,size(PointNames));
    azF     = repmat(azF(1)    ,size(PointNames));
%     noiseHz = repmat(noiseHz(1),size(PointNames));
%     noiseV  = repmat(noiseV(1) ,size(PointNames));
end

newSets = {Set0}; newSets{nTimes} = NaN;
for t = 1 : (nTimes - 1)
    set_new = Set0;                 
    set_new.Time = set_new.Time + dt;
    
    for nms = 1:length(PointNames)
        for i = 1:size(set_new.Points,1)
            de = vFhz(nms) * dt * sin(azF(nms));
            dn = vFhz(nms) * dt * cos(azF(nms));
            du = vFv(nms)  * dt           ;
            if (strcmpi(set_new.Points{i,1}{1},PointNames{nms}) == true)
                %updating crds
                set_new.crds(i,2) =  set_new.crds(i,2) + de;
                set_new.crds(i,1) =  set_new.crds(i,1) + dn;
                set_new.crds(i,3) =  set_new.crds(i,3) + du;
                
                % rectifying vectors 
                % when from Point is moving
                p = string(set_new.VectorsAndVCVs{:,1})==string(PointNames{nms});
                set_new.VectorsAndVCVs{ (p == 1) ,3:5} = set_new.VectorsAndVCVs{ (p == 1) ,3:5} - [dn,de,du];
                
                % when to Point is moving
                p = string(set_new.VectorsAndVCVs{:,2})==string(PointNames{nms});
                set_new.VectorsAndVCVs{ (p == 1) ,3:5} = set_new.VectorsAndVCVs{ (p == 1) ,3:5} + [dn,de,du];
            end %if
        end %for i
    end %for nms
    
    % updaing the entire DataSet
    set_new.Points{:,2:4} = set_new.crds;
    set_new.vctrs         = set_new.VectorsAndVCVs{:,3:5};
    
    % adding noise
    rng('shuffle');  %rng(1234); 
    set_new.vcvs         = set_new.vcvs         + mean([noiseHz,noiseV]) .* rand(size(set_new.vcvs));                   % ------ add the v-Noise to the rnadomization
    set_new.vctrs(:,1:2) = set_new.vctrs(:,1:2) +      noiseHz           .* rand(size(set_new.vctrs,1),2);
    set_new.vctrs(:,3)   = set_new.vctrs(:, 3 ) +      noiseV            .* rand(size(set_new.vctrs,1),1);
    set_new.VectorsAndVCVs{:,3:end} = [set_new.vctrs,set_new.vcvs];
    newSets{t+1} = set_new;
end % for Times
end% end function