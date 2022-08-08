function [set_new] = RemoveSites(set,siteNames2Remove)
% this function will keep only sites, and vectors which 
% not consist of sites (plural) given
%
% siteNames2Remove is a cell array
set_new = set;
for i = 1:numel(siteNames2Remove)
    set_new =  RemoveSite(set_new,siteNames2Remove{i});
end % for

end %function