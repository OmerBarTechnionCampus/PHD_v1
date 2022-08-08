function [nidx] = ind2sub_wNeightbours(siz,ndx,ngh)
%IND2SUB Multiple subscripts from linear index.
%   IND2SUB is used to determine the equivalent subscript values
%   corresponding to a given single index into an array.
%
%   [I,J] = IND2SUB(SIZ,IND) returns the arrays I and J containing the
%   equivalent row and column subscripts corresponding to the index
%   matrix IND for a matrix of size SIZ.  
%   For matrices, [I,J] = IND2SUB(SIZE(A),FIND(A>5)) returns the same
%   values as [I,J] = FIND(A>5).
%
%   [I1,I2,I3,...,In] = IND2SUB(SIZ,IND) returns N subscript arrays
%   I1,I2,..,In containing the equivalent N-D array subscripts
%   equivalent to IND for an array of size SIZ.
%
%   Class support for input IND:
%      float: double, single
%      integer: uint8, int8, uint16, int16, uint32, int32, uint64, int64
%
%   See also SUB2IND, FIND.
 
%   Copyright 1984-2015 The MathWorks, Inc. 
%   modified by Omer Bar 2020-09

nout = max(numel(siz),1);
siz = double(siz);
lensiz = length(siz);

if lensiz < nout
    siz = [siz ones(1,nout-lensiz)];
elseif lensiz > nout
    siz = [siz(1:nout-1) prod(siz(nout:end))];
end

if nout > 2
    k = cumprod(siz);
    for i = nout:-1:3
        vi = rem(ndx-1, k(i-1)) + 1;
        vj = (ndx - vi)/k(i-1) + 1;
        var_out{i} = double(vj);
        ndx = vi;
    end
end

if nout >= 2
    vi = rem(ndx-1, siz(1)) + 1;
    var_out{2} = double((ndx - vi)/siz(1) + 1);
    var_out{1} = double(vi);
end

% adding neightbours
% each row is a sub-index, each coloumn is the value
nidx = zeros(nout,2*ngh+1);
nidx(:,ngh+1) = cell2mat(var_out);
for i = 1:ngh
    nidx(:,ngh+1+i) = nidx(:,ngh+1) + i;
    nidx(:,ngh+1-i) = nidx(:,ngh+1) - i;
end

% removig entire colouns for partial existance of nieghbour (sub-index less then 1)
nidx = nidx(:,sum(nidx>0)     == size(nidx,1));
% removig entire colouns for partial existance of nieghbour (sub-index greater than max index)
nidx = nidx(:,sum((nidx <= siz') == 0) == 0);

nidx = nidx';
end %function