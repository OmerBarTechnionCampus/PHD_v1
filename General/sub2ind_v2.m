function ndx = sub2ind_v2(siz,varargin)
%SUB2IND Linear index from multiple subscripts.
%   SUB2IND is used to determine the equivalent single index
%   corresponding to a given set of subscript values.
%
%   IND = SUB2IND(SIZ,I,J) returns the linear index equivalent to the
%   row and column subscripts in the arrays I and J for a matrix of
%   size SIZ. 
%
%   IND = SUB2IND(SIZ,I1,I2,...,IN) returns the linear index
%   equivalent to the N subscripts in the arrays I1,I2,...,IN for an
%   array of size SIZ.
%
%   I1,I2,...,IN must have the same size, and IND will have the same size
%   as I1,I2,...,IN. For an array A, if IND = SUB2IND(SIZE(A),I1,...,IN)),
%   then A(IND(k))=A(I1(k),...,IN(k)) for all k.
%
%   Class support for inputs I,J: 
%      float: double, single
%      integer: uint8, int8, uint16, int16, uint32, int32, uint64, int64
%
%   See also IND2SUB.

%   Copyright 1984-2015 The MathWorks, Inc.

siz = double(siz);
lensiz = length(siz);
if lensiz < 2
    error(message('MATLAB:sub2ind:InvalidSize'));
end
numOfIndInput = nargin-1;
if numOfIndInput==1
    if isa(varargin{1},'double')== true
        idx = varargin{1};
        for k=1:size(idx,2)
            varargin{k} = idx(:,k);
        end
    else
        if isa(varargin{1},'cell')== true
           idx = varargin{1};
            for k=1:size(idx,2)
                varargin{k} = idx{:,k};
            end 
        else
            error(message('MATLAB:sub2ind:InvalidDataType - data type not handled'));
        end
    end
    % re omput size of input data after type conversion 
    numOfIndInput = numel(varargin);
end

% if lensiz < numOfIndInput
%     %Adjust for trailing singleton dimensions
%     siz = [siz, ones(1,numOfIndInput-lensiz)];
% elseif lensiz > numOfIndInput
%     %Adjust for linear indexing on last element
%     siz = [siz(1:numOfIndInput), prod(siz(numOfIndInput:end))];
% end

ndx = double(varargin{1}); % first index_argument_in
s = size(varargin{1});     % first index_argument_in size

%Compute linear indices
k = [0, cumprod(siz)];


for i = 1:numOfIndInput
    v1 = varargin{i};
    if any(min(v1(:)) < 1) || any(max(v1(:)) > siz(i))
        %Verify subscripts are within range
        error(message('MATLAB:sub2ind:IndexOutOfRange'));
    end %if 
    if ~isequal(s,size(v1))
        error(message('MATLAB:sub2ind:SubscriptVectorSize'));
    end %if
    
    %Compute linear indices
     ndx = ndx + (double(v1)-1)*k(i);
end  % for



end % function