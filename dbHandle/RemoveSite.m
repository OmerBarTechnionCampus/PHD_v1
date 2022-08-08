function [set_new] = RemoveSite(set,siteName2Remove)
% this function will keep only sites, and vectors which 
% not consist of Site given
%
% Omer Bar 2021, April;    version 1

ip = ~strcmpi(set.Points{:,1},siteName2Remove);
set_new.Points = set.Points(ip,:);
set_new.crds = set_new.Points{:,2:end};

icFrom = ~strcmpi(set.VectorsAndVCVs{:,1},siteName2Remove);
icTo   = ~strcmpi(set.VectorsAndVCVs{:,2},siteName2Remove);
ic = icFrom + icTo;
set_new.VectorsAndVCVs = set.VectorsAndVCVs(ic==2,:);
set_new.vcvs = set.VectorsAndVCVs{ic==2,end-6:end};
set_new.vctrs = set.VectorsAndVCVs{ic==2,3:5};
if (numel(set.Time) > 1)
    % for multiEpochData
    set_new.Time = set.Time{ic==2};
else
    % on single epoch
    set_new.Time = set.Time;
end

end %function