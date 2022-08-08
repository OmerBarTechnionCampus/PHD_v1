function [pos] = PointPositionToLine(e1,n1,e2,n2,ep,np)
% 0 - on line, -1 left, +1 right
% based upon https://math.stackexchange.com/questions/274712/calculate-on-which-side-of-a-straight-line-is-a-given-point-located
pos = sign( (ep-e1).*(n2-n1) - (np-n1).*(e2-e1) );
end