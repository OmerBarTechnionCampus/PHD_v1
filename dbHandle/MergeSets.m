function [mrgSet] = MergeSets(Sets, isDoAVG)
% merging sets....
% a boolean is there is we wish to preform avaraged between the data
%
% Omer Bar, 2021-05, version 1.0
if numel(Sets) == 1
    mrgSet = Sets;return;
end

warning('MergeSets :: time of the new sets are being disregarded - first set time is the time');
mrgSet = Sets{1};
for s = 2 : numel(Sets)
    curSet = Sets{s};
    for p = 1 : size(curSet.Points,1)
%         disp(p); disp(string(curSet.Points{p,1}));
        pPM = find(string(curSet.Points{p,1})==mrgSet.Points{:,1}); % searching point name
        
        if isempty(pPM) % a new point to the data set
            % point handeling
            mrgSet.Points(end+1,:) = curSet.Points(p,:);           
            
            % vector handeling
            pA = (string(curSet.VectorsAndVCVs{:,1})==curSet.Points{p,1}); % Vector's FromPoint
            pB = (string(curSet.VectorsAndVCVs{:,2})==curSet.Points{p,1}); % Vector's   ToPoint
        
            pVc = find(  (pA + pB) > 0  ); % vectors to add to the new dataset [one of the points is the new point]
            mrgSet.VectorsAndVCVs((end+1):(end+numel(pVc)),: ) = curSet.VectorsAndVCVs(pVc,:);
            % zeroizing the currently added vectors - not need to be treating them again
            curSet.VectorsAndVCVs{pVc,1} = {''};curSet.VectorsAndVCVs{pVc,2} = {''};curSet.VectorsAndVCVs{pVc,3:end} = NaN;
        else % point is already in the dataset
            % vector handeling
            pA = (string(curSet.VectorsAndVCVs{:,1})==mrgSet.Points{pPM,1}); % Vector's FromPoint
            pB = (string(curSet.VectorsAndVCVs{:,2})==mrgSet.Points{pPM,1}); % Vector's   ToPoint
        
            pVc = find(  (pA + pB) > 0  ); % vectors to add to the new dataset [one of the points is the new point]
            for v = 1:numel(pVc)
                % getting existing vectors 
                pA = (string(mrgSet.VectorsAndVCVs{:,1})==curSet.VectorsAndVCVs{pVc(v),1}); % Vector's FromPoint
                pB = (string(mrgSet.VectorsAndVCVs{:,2})==curSet.VectorsAndVCVs{pVc(v),2}); % Vector's   ToPoint
                pVm = find(  (pA+pB)>1  ,1); % same vectors as current to be added vector
                
                if ~isempty(pVm)
                    % if vector allready is in - check if the direction is same.... (from -> to and not other way around)
                    % if backwards -> flip vector direction
                    if ~strcmpi(curSet.VectorsAndVCVs{pVc(v),1} ,mrgSet.VectorsAndVCVs{pVm,1}) % checking names - if same "from point" than they are at same direction
                        % need to reverse the vector
                            temp = curSet.VectorsAndVCVs{pVc(v),1};
                        curSet.VectorsAndVCVs{pVc(v),1} = curSet.VectorsAndVCVs{pVc(v),2};
                        curSet.VectorsAndVCVs{pVc(v),2} = temp;
                        curSet.VectorsAndVCVs{pVc(v),3:5} = -curSet.VectorsAndVCVs{pVc(v),3:5};
                    end
                else
                    % this is a new vector (was not avialable beofre - nothing to check...
                end
            end % for
            mrgSet.VectorsAndVCVs(end+1:end+numel(pVc),: ) = curSet.VectorsAndVCVs(pVc,:);
            % zeroizing the currently added vectors - not need to be treating them again
            curSet.VectorsAndVCVs{pVc,1} = {''};curSet.VectorsAndVCVs{pVc,2} = {''};curSet.VectorsAndVCVs{pVc,3:end} = NaN;
        end
        
    end % points
end % for s = 2 : numel(numel(sets))
%%
if (isDoAVG == true)
    % mean will be computed for the VCVs and a weighted mean for the vectors themselves
    for v = 1:size(mrgSet.VectorsAndVCVs,1)
        if isnan(mrgSet.VectorsAndVCVs{v,3})
            continue; 
        end; % data was already taken cared of
        
        pA = (string(mrgSet.VectorsAndVCVs{:,1})== mrgSet.VectorsAndVCVs{v,1}); % Vector's FromPoint 
        pB = (string(mrgSet.VectorsAndVCVs{:,2})== mrgSet.VectorsAndVCVs{v,2}); % Vector's   ToPoint 
        p = ( (pA+pB)>1 ); 
        i = find(p,1);
%         disp(find(p)); disp(mrgSet.VectorsAndVCVs(p,:));
        
        % only 1 vector like that - no need to compute AVG
        if (numel(find(p)) == 1)
            continue; 
        end;  
        
        vcvs  = mrgSet.VectorsAndVCVs{p,end-5:end}; % order: covNN covNE covNU covEE covEU covUU
        vctrs = mrgSet.VectorsAndVCVs{p,  3  :  5}; % order: dN, dE, dU
        mrgSet.VectorsAndVCVs{p,3:end} = NaN;   % reducing the data
            p(i) = 0; % to not zeroize the one needed vector names
        mrgSet.VectorsAndVCVs{p,1} = {''}; mrgSet.VectorsAndVCVs{p,2} = {''};  % reducing the data
        
        mrgSet.VectorsAndVCVs{i,end-5:end} = mean(vcvs); % mean of weights
        % computing weighted sum
        p = [vcvs(:,1),vcvs(:,end-2),vcvs(:,end)];
        pPM = p.*vctrs;
        mrgSet.VectorsAndVCVs{i,3:5} = sum(pPM)./sum(p); % weighted sum
    end % for v = 1:size(mrgSet.VectorsAndVCVs,1)
    
    n = ~isnan(mrgSet.VectorsAndVCVs{:,3});
    mrgSet.VectorsAndVCVs = mrgSet.VectorsAndVCVs(n==1,:);
end % if (isDoAVG == true)
%%
% fixing non table data
mrgSet.crds = mrgSet.Points{:,2:4};
mrgSet.vctrs = mrgSet.VectorsAndVCVs{:,3:5};
mrgSet.vcvs = mrgSet.VectorsAndVCVs{:,6:end};
end % --- function ----