function [Sets] = ReadDataSetsFromBerneseFiles(inCRDpaths, inVCVpaths)
% example:
% [Sets] = ReadDataSetsFromBerneseFiles({'C:\GPSDATA\CAMPAIGN52\G1\STA\RED_G1.CRD'}, {'C:\GPSDATA\CAMPAIGN52\G1\OUT\G1FULL.COV'})

Sets{numel(inCRDpaths)} = {}; %pre-allocation of memory
centralSiteNumber = -1;
for i = 1:length(inCRDpaths)
    points = table({'NaN'},NaN,NaN,NaN,'VariableNames',{ 'Name' 'Grid_Northing' 'Grid_Easting' 'Height'});
    vectors = table({'NaN'},{'NaN'},NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,'VariableNames',{ 'Point_From' 'Point_To' 'dN' 'dE' 'dHt' 'covNN' 'covNE' 'covNU' 'covEE' 'covEU' 'covUU'});
    
    % ---- working on sites ----
    [sites] = ReadDataFromCRD(inCRDpaths{i});

    for k = 1:length(sites)
        if strcmpi('csar',sites{k}.Name)
            centralSiteNumber = k; 
            break;
        end
    end
    % preparing the transformation matrix...
    lat = sites{centralSiteNumber}.N*pi/180;
    lon = sites{centralSiteNumber}.E*pi/180;
    S = calcS(lat,lon);

    points(length(sites),end) = {-1}; %pre-allocate memory
    t = 0;
    for k = 1:length(sites)
        if ~( strcmpi('w',sites{k}.Flag) || strcmpi('a',sites{k}.Flag) )
            continue; % move to next point if this point was not estimated / computed
                      % it is just listed.....
        else 
            t = t + 1;
        end
        points{t,1} = {sites{k}.Name};
        if (k == centralSiteNumber)
             points{t,2:4} = [0,0,0];
        else
            uen = S * ([sites{k}.X;sites{k}.Y;sites{k}.Z] - [sites{centralSiteNumber}.X;sites{centralSiteNumber}.Y;sites{centralSiteNumber}.Z]);
            points(t,2:4) = {uen(3),uen(2),uen(1)}; % N E U
        end
    end %for i
    points = points(1:t,:);
    Set.Points = points;
    Set.crds = points{:,2:end};  
    
    % ---- working on vectors ----
    % reading the data
%     [vcv, rmsunit, obsnum, unkowns, index] = BuildVarCoVarFromCOVfile(inVCVpaths{i});
    [vcv, ~, ~, ~, index] = BuildVarCoVarFromCOVfile(inVCVpaths{i});
    
    n = length(points.Name);
    n = n*(n-1)/2; %number of max vectors between all sites
    vectors{1:n,3} = NaN*ones(n,1); %pre-allocate memory
    
    for j = 1:(length(points.Name)-1)
        for k = (j+1):length(points.Name)
            c = find(isnan(vectors{:,3}),1);
            if isempty(c)
                %pre-allocate for additional memory
                vectors{n+1:n+1000,3} = NaN*ones(numel(n+1:n+1000),1); %pre-allocate more memory to the table
                c = find(isnan(vectors{:,3}),1);
            end
            vectors{c,1}={points(j,:).Name{1}};  vectors{c,2}={points(k,:).Name{1}};
%             dxyz = [sites{k}.X;sites{k}.Y;sites{k}.Z] - [sites{j}.X;sites{j}.Y;sites{j}.Z];
%             duen = S * dxyz;
%             vectors{c,3:5} = duen';
            vectors{c,3:5} = points{k,2:4} - points{j,2:4};
            
            p1 = find(string(index{:,2})=={points(j,:).Name{1}},1);
            p2 = find(string(index{:,2})=={points(k,:).Name{1}},1);
            if (isempty(p1)) || (isempty(p2)); continue; end;
            
            tmpVCV = vcv( p1:(p1+2) , p2:(p2+2) );
            tmpVCV = S*tmpVCV; % U E N
            
                             % 'covNN'       'covNE'      'covNU'       'covEE'      'covEU'       'covUU'
            vectors(c,6:end) = {tmpVCV(3,3), tmpVCV(3,2),  tmpVCV(3,1), tmpVCV(2,2), tmpVCV(2,1), tmpVCV(1,1)};
        end % for k
    end % for j
    
    % removing "blank" lines [NaN lines]
    c = find(isnan(vectors{:,3}),1);
    if ~isempty(c)
        vectors = vectors(1:(c-1),:);
    end
    
    Set.VectorsAndVCVs = vectors;
    Set.vctrs = vectors{:,3:5};
    Set.vcvs = vectors{:,6:end};
    
    % ------ setting the epoch time -------
    Set.Time = sites{1}.Year + sites{1}.DOY/366.25;
    
    Sets{i} = Set;
end % for i = 1: length(in_paths)


end % function