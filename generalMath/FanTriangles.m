function [t] = FanTriangles(pnts)
% this function will compute set of triangle to a fan spread
%
%
% Omer Bar 2021-01
%% appprox. center of gravity
meanP = mean(pnts);

%% dividing points to two groups - convex_hull and internal points
% convex hull
ch = convhull(pnts(:,1:2));   % column vector
% hull = pnts(ch,:); % extracting the convex hull points 
% keeping points within the hull
i = ones(numel(pnts),1);
i(ch) = 0; %
i = find(i);
% inner = pnts(i,:);

%% getting mid point to create trinagles from
if numel(i) > 0
    % when there are points within the convex hull
    P2 = pnts; P2(ch,:) = NaN;
    diff = P2 - meanP;
    diff = diff(:,1).*diff(:,1) + diff(:,2).*diff(:,2) + diff(:,3).*diff(:,3);
    diff = sqrt(diff);

    mp = find(diff == min(diff),1);
else
    % when the entire point set is the convex hull
    diff = pnts - meanP;
    diff = diff(:,1).*diff(:,1) + diff(:,2).*diff(:,2) + diff(:,3).*diff(:,3);
    diff = sqrt(diff);

    mp = find(diff == min(diff),1);
    
    % removing point from hull
    ch(mp) = NaN;
    ch = ch(~isnan(ch));
    if (  (numel(pnts,1) - numel(ch) ) == 1  )
        % when all points are the hull , numel(ch) === numel(pnts(:,1)) + 1
        % if (checkup) true, meaning initial hull point removed (twice - from start and end)
        % need to fix the indexes
        ch(end+1) = ch(1);
    end
end %if
% midP = pnts(mp,:);

%% building the triangles from the mid point
ch = ch(1:end-1);
ch2 = [ch(end);ch(1:end-1)];
MP = ones(size(ch)).*mp;
t = [ch,ch2,MP];

t = sortrows(sort(t,2)); % arranging numbers of points (small to large, in each row, and than all rows)
end % function