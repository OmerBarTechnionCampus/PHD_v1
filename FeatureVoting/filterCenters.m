function [Cf]=filterCenters(C,idx, strMode, Az_th, V_th, L_th, pos_th)
% C are centers from kmeans clustering
% column order is Az, V, depth, posE, posN
% idx is the conectivity vector - each point belongs to which center
%
% Cf column order: 
% Az, V, depth, posE, posN, popularity
%% computations
uidx = unique(idx); % centers indexes
for k = 1:numel(uidx)
    uidx(k) = sum(idx == uidx(k));
end
puidx = uidx / numel(idx) * 100; % precentages ["popularity"]

% adding the puidx to the Centers matrix
% then sorting Centers in descding order by puidx
C = sortrows([C , puidx], -(size(C,2)+1));  

if strcmpi(strMode,'SingleFault')
    % filtering by presentage by the mean
    C = C( C(:,end)>= mean(C(:,end)),:); 
    
    % checking for similar geographic locations    
    % adding first solution (the best one)
    Cf = NaN .* ones(size(C));
    
    for k = 1:size(C,1)
        if isnan(C(k,1)) ; continue; end;
        kf = find(isnan(Cf(:,1)),1);
        Cf(kf,:) = C(k,:); % set the first NaN row value as the curret suggested center
        C(k,:) = NaN;

        i = find( ~isnan(C(:,1)) );
        az = Azimuth(C(i,end-1) - Cf(kf,end-1) , C(i,end) - Cf(kf,end));
        daz = Cf(kf,1) - az;
        i = i .* ( abs(daz) <= Az_th ) ;
        if (sum(i) > 0) % centers were found
            C(i==1,:) = NaN; % marking for deleting the solutions
        end
    end % for k = 2:size(C,1)
    
    Cf = Cf(  ~isnan(Cf(:,1)) ,:); % removing the NaN values , and the last column ["popularity"]
else % MultiFault
    Cf = NaN;
    error('filterCenters :: thi section is unhandeled')
end % if - mode



end % function