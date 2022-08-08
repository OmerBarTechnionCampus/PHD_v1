function [x,y,z,radii] = GetBoundingCircle(x1,y1,z1,x2,y2,z2,x3,y3,z3)
%  assuming that the edges of the triangle are: 1-2,2-3,3-1
%
% example
%
%
% Omer Bar v1.0 (NOvember 2017)
%% prepinn the data
x1 = reshape(x1,[max(size(x1)),1]);
x2 = reshape(x2,[max(size(x2)),1]);
x3 = reshape(x3,[max(size(x3)),1]);

y1 = reshape(y1,[max(size(y1)),1]);
y2 = reshape(y2,[max(size(y2)),1]);
y3 = reshape(y3,[max(size(y3)),1]);

z1 = reshape(z1,[max(size(z1)),1]);
z2 = reshape(z2,[max(size(z2)),1]);
z3 = reshape(z3,[max(size(z3)),1]);

if isnan(z1)
    z1 = zeros(size(x1));
end
if isnan(z2)
    z2 = zeros(size(x1));
end
if isnan(z3)
    z3 = zeros(size(x1));
end

%% entry tests
if ( (size(x1,1)~=size(x2,1)) || (size(x3,1)~=size(x2,1)) || (size(x1,1)~=size(x3,1)) )
    error('x data is not at same sizes');
end
if ( (size(y1,1)~=size(y2,1)) || (size(y3,1)~=size(y2,1)) || (size(y1,1)~=size(y3,1)) )
    error('y data is not at same sizes');
end
if ( (size(z1,1)~=size(z2,1)) || (size(z3,1)~=size(z2,1)) || (size(z1,1)~=size(z3,1)) )
    error('z data is not at same sizes');
end
if ( (size(x1,1)~=size(y1,1)) || (size(x1,1)~=size(z1,1)) || (size(y1,1)~=size(z1,1)) )
    error('point 1 data is not at same size');
end
if ( (size(x2,1)~=size(y2,1)) || (size(x2,1)~=size(z2,1)) || (size(y2,1)~=size(z2,1)) )
    error('point 1 data is not at same size');
end
if ( (size(x3,1)~=size(y3,1)) || (size(x3,1)~=size(z3,1)) || (size(y3,1)~=size(z3,1)) )
    error('point 1 data is not at same size');
end

%% actual calculations
try
    p1 = [x1,y1,z1];
    p2 = [x2,y2,z2];
    p3 = [x3,y3,z3];
    [center,rad,~,~,msgs] = circlefit3d(p1,p2,p3);
    if ((rad < 0) || (isempty(msgs)==false))
        x       = ones(size(x1,1),1).*NaN;
        y       = x;
        z       = x;
        radii   = x;
        return;
    end
    
    x       = center(:,1);
    y       = center(:,2);
    z       = center(:,3);
    radii   = rad;
catch %excp
    x       = ones(size(x1,1),1).*NaN;
    y       = x;
    z       = x;
    radii   = x;
end
end % function    GetBoundingCircle