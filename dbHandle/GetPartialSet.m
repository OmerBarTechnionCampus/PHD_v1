function [partialSet] = GetPartialSet(Set,PointNames)
% version 2
    partialSet.Time = Set.Time;
    vs = zeros(size(Set.VectorsAndVCVs,1),1);
    
    partialSet.Points(1,:) = Set.Points(1,:);
    partialSet.Points(1,:).Name{1} = 'NaN'; partialSet.Points{1,2:end} = NaN.*ones(1,size(partialSet.Points,2)-1);
    for i = 1:numel(PointNames)
        name = PointNames{i};
        if (isempty(name) == true); continue; end;
        partialSet.Points(end+1,:) = Set.Points(find(string(Set.Points{:,1}) == name,1),:);
        
        f11 = find(string(Set.VectorsAndVCVs{:,1})==name) ; f12 = find(string(Set.VectorsAndVCVs{:,2})==name);
        vs(f11) = vs(f11) +1;  vs(f12) = vs(f12) +1;
    end
    partialSet.Points = partialSet.Points(2:end,:);
    partialSet.crds = partialSet.Points{:,2:end};
    
    partialSet.VectorsAndVCVs(1:numel(find(vs>1)),:) = Set.VectorsAndVCVs(vs>1,:);
    partialSet.vcvs = partialSet.VectorsAndVCVs{:,6:end};
    partialSet.vctrs = partialSet.VectorsAndVCVs{:,3:5};
end

%{
version 1

function [partialSet] = GetPartialSet(Set,namePnt1,namePnt2,namePnt3)
    partialSet.Time = Set.Time;
    partialSet.Points(1,:) = Set.Points(find(string(Set.Points{:,1}) == namePnt1,1),:);
    partialSet.Points(2,:) = Set.Points(find(string(Set.Points{:,1}) == namePnt2,1),:);
    partialSet.Points(3,:) = Set.Points(find(string(Set.Points{:,1}) == namePnt3,1),:);
    partialSet.crds = partialSet.Points{:,2:end};
    
    f11 = find(string(Set.VectorsAndVCVs{:,1})==namePnt1) ; f12 = find(string(Set.VectorsAndVCVs{:,2})==namePnt1);
    f21 = find(string(Set.VectorsAndVCVs{:,1})==namePnt2) ; f22 = find(string(Set.VectorsAndVCVs{:,2})==namePnt2);
    f31 = find(string(Set.VectorsAndVCVs{:,1})==namePnt3) ; f32 = find(string(Set.VectorsAndVCVs{:,2})==namePnt3);
    
    vs = zeros(size(Set.VectorsAndVCVs,1),1);
    vs(f11) = vs(f11) +1;  vs(f12) = vs(f12) +1;
    vs(f21) = vs(f21) +1;  vs(f22) = vs(f22) +1;
    vs(f31) = vs(f31) +1;  vs(f32) = vs(f32) +1;
    
    if (sum(vs > 1) < 3)
        partialSet = NaN;
        return;
    else
        partialSet.VectorsAndVCVs(1:3,:) = Set.VectorsAndVCVs(vs>1,:);
        partialSet.vcvs = partialSet.VectorsAndVCVs{:,6:end};
        partialSet.vctrs = partialSet.VectorsAndVCVs{:,3:5};
    end
end
%}