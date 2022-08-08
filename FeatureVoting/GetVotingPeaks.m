function [votingPeaks_index] = GetVotingPeaks(H, maxPeaks)
% this function computes the K-peaks for a voting matrix
% Voting values are always >=0
%    example:   Cp_idx = GetVotingPeaks(H5, 10);
% Omr bar, 2020.09.29, version 1.0
%% Computation
% setting nan value
votingPeaks_index = NaN;
%% removing un-wanted data
Hsz = size(H);
maxH = max(H(:));% H = H(1:numel(H)); 

while max(H(:)) > 256
    H = H-256;
end
H(H<=0) = 0;

thresh = graythresh(H(1:numel(H)));  %the reshape is for the graythresh function only
% fixing thresh for the normalization...
thresh = thresh * max(H(:));

H(H<=thresh) = 0;

if (max(H(:))==0)
    return;
end

%% getting peaks iteratively
% stdev = std(H(H>0))  / 10  ;
neighbours = 3;
votingPeaks_index = NaN.*ones(maxPeaks,1);

iter = 0;
max_iter = maxPeaks * 2;
while (sum(find(isnan(votingPeaks_index))) > 0) % counting the NaN values
    iter = iter +1;
    [v,i] = max(H(:));
    if (v == 0)
        break;
    end
    
    t = find(isnan(votingPeaks_index),1);
    votingPeaks_index(t) =   i;
    
    %removing peak and its relative close values
       %H(H>=v-stdev) = 0;    
    sub_idx = ind2sub_wNeightbours(Hsz,i,neighbours);
    [idx_wNghbr] = sub2ind_v2(Hsz, sub_idx );
    H(idx_wNghbr) = 0;
    
    if (iter == max_iter) ; break; end;
end % while

%removing nan values
votingPeaks_index=votingPeaks_index(~isnan(votingPeaks_index));
votingPeaks_index = unique(votingPeaks_index);
end % function








% function [Is] = getCloserIndexs(idx,sze,nghbr)
% % this function will get the closest neighbors to the idx sent
% switch numel(sze)
%     case 2
%         [i1,i2] = ins2sub(sze,idx);
%         i1 = i1-ngbr:1:i1+ngbr;
%         i2 = i2-ngbr:1:i2+ngbr;
%         
%         Is = sub2ind(sze,i1,i2);
%     case 3
%         [i1,i2,i3] = ins2sub(sze,idx);
%         i1 = i1-ngbr:1:i1+ngbr;
%         i2 = i2-ngbr:1:i2+ngbr;
%         i3 = i3-ngbr:1:i3+ngbr;
%         
%         Is = sub2ind(sze,i1,i2,i3);
%     case 4
%         [i1,i2,i3,i4] = ins2sub(sze,idx);
%         i1 = i1-ngbr:1:i1+ngbr;
%         i2 = i2-ngbr:1:i2+ngbr;
%         i3 = i3-ngbr:1:i3+ngbr;
%         i4 = i4-ngbr:1:i4+ngbr;
%         
%         Is = sub2ind(sze,i1,i2,i3,i4);
%     case 5
%         [i1,i2,i3,i4,i5] = ins2sub(sze,idx);
%         i1 = i1-ngbr:1:i1+ngbr;
%         i2 = i2-ngbr:1:i2+ngbr;
%         i3 = i3-ngbr:1:i3+ngbr;
%         i4 = i4-ngbr:1:i4+ngbr;
%         i5 = i5-ngbr:1:i5+ngbr;
%         
%         Is = sub2ind(sze,i1,i2,i3,i4,i5);
%     case 6
%         [i1,i2,i3,i4,i5,i6] = ins2sub(sze,idx);
%         i1 = i1-ngbr:1:i1+ngbr;
%         i2 = i2-ngbr:1:i2+ngbr;
%         i3 = i3-ngbr:1:i3+ngbr;
%         i4 = i4-ngbr:1:i4+ngbr;
%         i5 = i5-ngbr:1:i5+ngbr;
%         i6 = i6-ngbr:1:i5+ngbr;
%         
%         Is = sub2ind(sze,i1,i2,i3,i4,i5,i6);
%     otherwise
%         error('getCloserIndexs::un handeled demention size');
% end %switch
% end %function
