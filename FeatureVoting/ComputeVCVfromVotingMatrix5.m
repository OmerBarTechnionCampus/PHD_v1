function [vcv_stat] = ComputeVCVfromVotingMatrix5(mat,index1,index2,v1,v2,v3)
% this function will comvute v-cv's / stdev's from the voting matrix distribution
% by the input index numbers

mat_size = size(mat);
dims = numel(size(mat));

% %%% version 1
% p = 1:dims;
% s1 = p(1); p(1) = index1; p(index1) = s1;
% s1 = p(2); p(2) = index2; p(index2) = s1;
% 
% mat = permute(mat,p);
% mat2 = mat(:,:,round(mean(v1)),round(mean(v2)),round(mean(v3)));


v1 = floor(mean(v1));
v2 = floor(mean(v2));
v3 = floor(mean(v3));

% mat2 = nan.*ones(mat_size(index1),mat_size(index2)); % pre-allocate memory

% p is a 0/1 valued matrix for cells which are needed for extraction -
% index matrix
p = zeros(numel(mat),dims);
p(:,   1:dims==index1   ) = 1;
p(:,   1:dims==index2   ) = 1;

p2 = nan*ones(1,dims);
p2(index1) = 0;
p2(index2) = 0;

p2(find(isnan(p2),1)) = v1;
p2(find(isnan(p2),1)) = v2;
p2(find(isnan(p2),1)) = v3;

for i1 = 1:mat_size(index1)
    for i2 = 1:mat_size(index2)
        p2(index1)=i1; p2(index2)=i2;
        n = sub2ind(mat_size,p2(1),p2(2),p2(3),p2(4),p2(5)); % this method works only on matrices with 5 dimensions
        %p(ind2sub(mat_size,n)) = p(ind2sub(mat_size,n)) +1 ;
        p(n,:) = p(n,:) + 1;
    end
end

n = sum(p,2);     % sum of each row
n = n==max(n(:)); % index of max values 

mat2 = mat(  ind2sub(mat_size,find(n==1))  );
mat2 = reshape(mat2,    [mat_size(index1),mat_size(index2)]    );

%transpose for imagesc
mat2 = mat2';

% var = sum(mat2(:))/numel(mat2);
bw = mat2./max(mat2(:));% > 0.95;
figure, imagesc(bw);
figure(gcf), colorbar;

bw = bw > 0.95;
figure, imagesc(bw);
vcv_stat = regionprops('table',bw,'Centroid','MajorAxisLength','MinorAxisLength','Orientation');

vcv_stat.dX = vcv_stat.MajorAxisLength.*cos(vcv_stat.Orientation);
vcv_stat.dY = vcv_stat.MajorAxisLength.*sin(vcv_stat.Orientation);
end





% run example 
% mat = abs(round(100*randn(3,5,6,10,3)));
% ComputeVCVfromVotingMatrix5(mat,2,3,mod(1,3),mod(5,10),mod(2,3))
% dbquit