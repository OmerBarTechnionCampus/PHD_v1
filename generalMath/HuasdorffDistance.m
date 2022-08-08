%% biDirectional Huasdorff Distance
function [FinalHuasdorffRes] = HuasdorffDistance(setA,setB)
% calculates the bi-directional hausdorff distance
%
%setA = [[x1,y1]
%		 [x2,y2]
%			...
%		 [xn,yn]
%
%setB = [[x1,y1]
%		 [x2,y2]
%			...
%		 [xn,yn]
%

FinalHuasdorffRes = single_huasdorff(setA(1,:),setB);

for i = 2 : 1 : length(setA)
    res = single_huasdorff(setA(i,:),setB);
    if (res > FinalHuasdorffRes) %max
        FinalHuasdorffRes = res;
    end
end

for i = 1 : 1 : length(setA)
    res = single_huasdorff(setB(i,:),setA);
    if (res > FinalHuasdorffRes) %max
        FinalHuasdorffRes = res;
    end
end

end %return FinalHuasdorffRes;

%% Single Direction Huasdorff Distance

function [res] = single_huasdorff(point, set)
% calculates the single-directional hausdorff distance
%point = [x,y]
%setA = [[x1,y1]
%		 [x2,y2]
%			...
%		 [xn,yn]
res = norm(point' - set(1,:)');
for i = 2 : 1 : length(set)
    d = norm(point' - set(i,:)');
    if (res > d) %% min
        res = d;
    end
end

end %return res;