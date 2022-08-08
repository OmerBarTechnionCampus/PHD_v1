function [set0, set1] = ArrangeMultiEpochSets(Sets)
% this function is designed to aggaine the multi-epoch data 
% to a structure similar to a single set
% changes will (a) sort the vector to a specific order, and 
% (b) for each vector to determine the time diffs by changing the ".Time"
% to a multi scalar value (yyyy.doy)
%
%
% Omer Bar 2021-02, version 1.0

%%
% initializing empty data sets
set0 = Sets{1};
%             Points: [8ª4 table]
%               crds: [8ª3 double]
%     VectorsAndVCVs: [28ª11 table]
%               Time: 2015.08767123288
%               vcvs: [28ª6 double]
%              vctrs: [28ª3 double]
% clearing data
set0.Points = []; set0.crds = []; set0.Time = []; set0.vcvs = []; set0.vctrs = []; 
set0.VectorsAndVCVs = set0.VectorsAndVCVs(1,:);
set0.VectorsAndVCVs{:,3:end} = NaN;  set0.VectorsAndVCVs{:,1:2} = {'', ''};

set1 = set0;

for k0 = 1:(numel(Sets) - 1)
    for k1 = (k0+1):numel(Sets)
        % current epochs in this iteration
        tmpSet0 = Sets{k0};   tmpSet1 = Sets{k1};
        
        % ========== working on points =============
        % gross initialization
        if ((k0 == 1) && (k1==2)) 
            set0.Points = tmpSet0.Points; 
            set1.Points = tmpSet1.Points;
        else
            % new epoch
            for c = 1:size(tmpSet0.crds,1)
                % if point from tmpSet0 is not in set0 list (new point from current epoch)
                if isPointInList(set0.Points{:,1},tmpSet0.Points{c,1}) == 0  
                    % add it to the list - will assist in comparing un-equal sets
                    set0.Points(end+1,:) = tmpSet0.Points(c,:); 
                end
                % if point is not in set1 list (new point from current epoch)
                if isPointInList(set1.Points{:,1},tmpSet0.Points{c,1}) == 0  
                    % add it to the list - will assist in comparing un-equal sets
                    set1.Points(end+1,:) = tmpSet0.Points(c,:); 
                end
            end % for c - for each point
            
            for c = 1:size(tmpSet1.crds,1)
                % if point is not in set0 list
                if isPointInList(set0.Points{:,1},tmpSet1.Points{c,1}) == 0  
                    % add it to the list - will assist in comparing un-equal sets
                    set0.Points(end+1,:) = tmpSet1.Points(c,:); 
                end
                % if point is not in set1 list
                if isPointInList(set1.Points{:,1},tmpSet1.Points{c,1}) == 0  
                    % add it to the list - will assist in comparing un-equal sets
                    set1.Points(end+1,:) = tmpSet1.Points(c,:); 
                end
            end % for c - for each point
        end;
       
        % updating data structrue - points
        set0.crds = set0.Points{:, 2:end};
        set1.crds = set1.Points{:, 2:end};
        if (numel(set0.crds) ~= numel(set1.crds))
            % cheking if all points in set0 exsits in set1
            for c = 1:numel(set0.Points{:, 1})
                % if point is not in set1 list
                if (isPointInList(set1.Points{:,1},set0.Points{c,1}) == 0 )
                    % add it to the list - will assist in comparing un-equal sets
                    set1.Points(end+1,:) = set0.Points(c,:); 
                end
            end
            
            % cheking if all points in set0 exsits in set1
            for c = 1:numel(set1.Points{:, 1})
                % if point is not in set1 list
                if (isPointInList(set0.Points{:,1},set1.Points{c,1}) == 0 )
                    % add it to the list - will assist in comparing un-equal sets
                    set0.Points(end+1,:) = set1.Points(c,:); 
                end
            end
        end
        
        set0.Points = sortrows(set0.Points,1);
        set1.Points = sortrows(set1.Points,1);
        
        % ======== working on vectors ============ 
        % getting them by relying on the point names
        % reason - vector might not exists bewteen two epochs (natural process)
        for c1 = 1:numel(set0.Points{:,1})
            for c2 = (c1+1):numel(set0.Points{:,1})
                if (c1==c2); continue; end;
                
                % tmpSet0
                % getting indexes for vectors which include c1 and c2 points
                c1f0 = strcmpi(tmpSet0.VectorsAndVCVs{:,1},set0.Points{c1,1}); % as from point
                c1t0 = strcmpi(tmpSet0.VectorsAndVCVs{:,2},set0.Points{c1,1}); % as  to  point
                c2f0 = strcmpi(tmpSet0.VectorsAndVCVs{:,1},set0.Points{c2,1}); % as from point
                c2t0 = strcmpi(tmpSet0.VectorsAndVCVs{:,2},set0.Points{c2,1}); % as from point
                v0 = find((c1f0 + c1t0 + c2f0 + c2t0) == 2,1);
                
                % tmpSet1
                % getting indexes for vectors which include c1 and c2 points
                c1f1 = strcmpi(tmpSet1.VectorsAndVCVs{:,1},set0.Points{c1,1}); % as from point
                c1t1 = strcmpi(tmpSet1.VectorsAndVCVs{:,2},set0.Points{c1,1}); % as  to  point
                c2f1 = strcmpi(tmpSet1.VectorsAndVCVs{:,1},set0.Points{c2,1}); % as from point
                c2t1 = strcmpi(tmpSet1.VectorsAndVCVs{:,2},set0.Points{c2,1}); % as from point
                v1 = find((c1f1 + c1t1 + c2f1 + c2t1) == 2,1);
                
                if (numel(v0)*numel(v1) == 0); continue; end; % vector not found on one of the epochs
                
                set0.VectorsAndVCVs(end+1,:) = tmpSet0.VectorsAndVCVs(v0,:);
                set1.VectorsAndVCVs(end+1,:) = tmpSet1.VectorsAndVCVs(v1,:);
                set0.Time(end+1,1) = tmpSet0.Time;
                set1.Time(end+1,1) = tmpSet1.Time;
                
                if (isnan(set0.VectorsAndVCVs{1,3})  && (size(set0.VectorsAndVCVs,1))>1)
                    % meaning the first row is an empty data row and needs
                    % to be deleted
                    set0.VectorsAndVCVs = set0.VectorsAndVCVs(2:end,:);
                    set1.VectorsAndVCVs = set1.VectorsAndVCVs(2:end,:);
                    % reoder the naming convenvtion
                    set1.VectorsAndVCVs(end,1:2) = set0.VectorsAndVCVs(end,1:2);
                end
                
                % checking if vector solution is backwards
                if ( ( sign(set0.VectorsAndVCVs{end,3}) ~= sign(set1.VectorsAndVCVs{end,3}) ) ...
                  || ( sign(set0.VectorsAndVCVs{end,4}) ~= sign(set1.VectorsAndVCVs{end,4}) ) ... 
                  || ( sign(set0.VectorsAndVCVs{end,5}) ~= sign(set1.VectorsAndVCVs{end,5}) ) )
                    % flipping one of the vectors        
                    set1.VectorsAndVCVs{end,3:5} = set1.VectorsAndVCVs{end,3:5} .* (-1);
                end
                
                
            end % for c1
        end % for c2
        
        % updating data structrue - vectors
        set0.vctrs = set0.VectorsAndVCVs{:,3:5};
        set0.vcvs  = set0.VectorsAndVCVs{:,6:end};
        set1.vctrs = set1.VectorsAndVCVs{:,3:5};
        set1.vcvs  = set1.VectorsAndVCVs{:,6:end};
    end % for k1
end %for k0
end %function


%% internal functions
function [tf] = isPointInList(strArrNames,name)
% return 0 if point is not in list, 1 if it is
    temp = strcmpi(strArrNames,name);
    tf = sum(temp) > 0; % if empty then will return 0 (false), else 1 (true)
end 


