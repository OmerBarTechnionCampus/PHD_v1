function [Sets] = ReadDataSetsFromGeolbFiles(in_paths)
Sets = {};
for i = 1:length(in_paths)
    points = table({'NaN'},NaN,NaN,NaN,'VariableNames',{ 'Name' 'Grid_Northing' 'Grid_Easting' 'Height'});
    vectors = table({'NaN'},{'NaN'},NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,'VariableNames',{ 'Point_From' 'Point_To' 'dN' 'dE' 'dHt' 'covNN' 'covNE' 'covNU' 'covEE' 'covEU' 'covUU'});
    
    % crds = readGeoLabCoordinatesFile(in_paths{i});
    % vctrs = readVectorsFromGeolabFile(in_paths{i});
    % fullCOV = readFullCOVfile_GeoLab(in_paths{i});
    [coordsWGS84,vectorsNEU] = readVectorsFromGeolabFile(in_paths{i});
    
%     pos = find(string(coordsWGS84{:,1}) == 'CSAR');
    pos = find(string(coordsWGS84{:,1}) == coordsWGS84.Name{2}); % same as 'CSAR' for main example
    lat = coordsWGS84{pos,2}*pi/180;
    lon = coordsWGS84{pos,3}*pi/180;
    S = calcS(lat,lon);
    
    [x0,y0,z0] = lla2ecef(lat,lon,coordsWGS84{pos,4});
    crds = zeros(size(coordsWGS84,1),3);
    for k = 1:size(crds,1)
        points{k,1} = coordsWGS84{k,1};
        if (k ~= pos)
            [tmpx,tmpy,tmpz] = lla2ecef(coordsWGS84{k,2}*pi/180,coordsWGS84{k,3}*pi/180,coordsWGS84{k,4});
            uen = S * ([tmpx;tmpy;tmpz] - [x0;y0;z0]);
            crds(k,:) = [uen(3),uen(2),uen(1)];
            points(k,2:4) = {uen(3),uen(2),uen(1)}; % N E U
        else
            crds(k,:) = [0,0,0];
            points(k,2:4) = {0,0,0};
        end
    end
    
    vectors(1:size(vectorsNEU,1),1:size(vectors,2)-6) = vectorsNEU(:,1:end-1);
    for k = 1:size(vectors,1)
        mat = cell2mat( vectorsNEU{k,end});
        vectors(k,end-5:end) = {mat(1,1),mat(1,2),mat(1,3),mat(2,2),mat(2,3),mat(3,3)};
    end
    
    Set.Points = points;
    Set.crds = crds;  
    Set.VectorsAndVCVs = vectors;
        str = in_paths{i};
        pos1= strfind(str,'\');
        pos2= strfind(str,'.');
        str = str(pos1(end)+1:pos2(end)-1);
        pos1 = strfind(str,'_');
        year = str2double(str(1:(pos1(1)-1)));
        day  = str2double(str(pos1(1)+1:(pos1(2)-1)));
    Set.Time = year + day/366.25; 
    Set.vcvs = vectors{:,6:end};
    Set.vctrs = vectors{:,3:5};
        Set.vcvs = sqrt(abs(Set.vcvs)).*sign(Set.vcvs);     % moving to milimeters
        Set.VectorsAndVCVs{:,end-5:end} = Set.vcvs;         % updating    
    Sets{end+1} = Set;
end % for i = 1 : length(in_paths)

end