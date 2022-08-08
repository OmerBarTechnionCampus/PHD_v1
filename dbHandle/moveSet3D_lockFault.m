function [newSets] = moveSet3D_lockFault(dt,ldF,vFhz,vFv,azF,posF,Set0,noiseHz,noiseV,nTimes)
%{
% updated 2019.07.28    -    wrong error
newSets = {Set0}; newSets{nTimes} = NaN;

minE = floor(min(Set0.crds(:,2))) - 3000;% minN = floor(min(Set0.crds(:,1))) - 3000;
maxE = floor(max(Set0.crds(:,2))) + 3000;% maxN = floor(max(Set0.crds(:,1))) + 3000;
m = tan(azF);
e1 = minE; n1 = m*(posF(1) - e1) + posF(2);
e2 = maxE; n2 = m*(posF(1) - e2) + posF(2);

for t = 1 : (nTimes - 1)
    set_new = Set0;                 
    set_new.Time = set_new.Time + dt;
    
    for nms = 1:size(Set0.Points,1)
        for i = 1:size(set_new.Points,1)
            e = set_new.crds(i,2); n = set_new.crds(i,1);
            l_r = PointSideFromLine(e1,n1,e2,n2,e,n);  % left or right from line
            l  = DistLinePoint(e1,n1,e2,n2,e,n);
            

            v = vFhz / pi * atan(l/ldF); % 4.4 formula - from PhD tyhsis Lior Shahar
            de = v*dt*sin(azF)       * l_r;
            dn = v*dt*cos(azF)       * l_r;
            du = vFv * dt            * l_r;
            
            if (strcmpi(set_new.Points{i,1}{1},Set0.Points{nms,1}{1}) == true)
                %updating crds
                set_new.crds(i,2) =  set_new.crds(i,2) + de;
                set_new.crds(i,1) =  set_new.crds(i,1) + dn;
                set_new.crds(i,3) =  set_new.crds(i,3) + du;
            end %if
        end %for i
    end %for nms
    
    % updaing the PointsDataSet
    set_new.Points{:,2:4} = set_new.crds;
    
    % updaing the VectorsDataSet
    for i = 1:size(set_new.VectorsAndVCVs,1)
        name1 = set_new.VectorsAndVCVs{i,1} {1};
        name2 = set_new.VectorsAndVCVs{i,2} {1};
        
        i1 = 0; i2= 0;
        for k = 1:size(set_new.Points,1)
            if (i1 >0) && (i2>0);                 break;                  end;
            if (strcmpi(name1,set_new.Points{k,1}{1})); i1 = k; continue; end;
            if (strcmpi(name2,set_new.Points{k,1}{1})); i2 = k; continue; end;
        end
        
        set_new.VectorsAndVCVs{i,3:5} = set_new.Points{i2,[3,2,4]} - set_new.Points{i1,[3,2,4]};
    end
    set_new.vctrs = set_new.VectorsAndVCVs{:,3:5};
    
    
    % adding noise
    rng('shuffle');  %rng(1234); 
    set_new.vcvs         = set_new.vcvs         + mean([noiseHz,noiseV]) .* rand(size(set_new.vcvs));                   % ------ add the v-Noise to the rnadomization
    set_new.vctrs(:,1:2) = set_new.vctrs(:,1:2) +      noiseHz           .* rand(size(set_new.vctrs,1),2);
    set_new.vctrs(:,3)   = set_new.vctrs(:, 3 ) +      noiseV            .* rand(size(set_new.vctrs,1),1);
    set_new.VectorsAndVCVs{:,3:end} = [set_new.vctrs,set_new.vcvs];
    newSets{t+1} = set_new;
end % for Times
%}


%%  2020.02.18  - tested ok +- 1mm
%{

% updated 2020.02.18
newSets = {Set0}; newSets{nTimes} = NaN;

e1 = posF(1) -10000.*sin(azF);   e2 = posF(1) +10000.*sin(azF);
n1 = posF(2) -10000.*cos(azF);   n2 = posF(2) +10000.*cos(azF);

for t = 1 : (nTimes - 1)
    set_new = Set0;                 
    set_new.Time = set_new.Time + dt;
    
    for nms = 1:size(Set0.Points,1)
        for i = 1:size(set_new.Points,1)
            if (strcmpi(set_new.Points{i,1}{1},Set0.Points{nms,1}{1}) == true)
                e = set_new.crds(i,2); n = set_new.crds(i,1);
                l_r = PointSideFromLine(e1,n1,e2,n2,e,n);  % left or right from line
                l  = DistLinePoint(e1,n1,e2,n2,e,n);


                v = vFhz / pi * atan(l/ldF); % 4.4 formula - from PhD tyhsis Lior Shahar
                de = v*dt*sin(azF)       * l_r;
                dn = v*dt*cos(azF)       * l_r;
                du = vFv * dt / (l*l)    * l_r;
            
                %updating crds
                set_new.crds(i,1) =  set_new.crds(i,1) + dn;
                set_new.crds(i,2) =  set_new.crds(i,2) + de;
                set_new.crds(i,3) =  set_new.crds(i,3) + du;
                continue;
            end %if
        end %for i
    end %for nms
    
    % rectifying vectors 
    for k = 1:size(Set0.VectorsAndVCVs,1)
        pAi = find(string(Set0.VectorsAndVCVs{k,1})==Set0.Points{:,1}); % vectors FromPoint
        pBi = find(string(Set0.VectorsAndVCVs{k,2})==Set0.Points{:,1}); % vectors   ToPoint
        
        % FromPoint - Dist from fualt
        lpA = round(DistLinePoint(  e1,n1,...
                                    e2,n2,...
                                    Set0.crds(pAi,2),Set0.crds(pAi,1)));
        % FromPoint - left or right from fualt
        l_r_A =PointSideFromLine(   e1,n1,...
                                    e2,n2,...
                                    Set0.crds(pAi,2),Set0.crds(pAi,1));
        % ToPoint -   Dist from fualt
        lpB = round(DistLinePoint(  e1,n1,...
                                    e2,n2,...
                                    Set0.crds(pBi,2),Set0.crds(pBi,1))); % e, n
        % ToPoint -   left or right from fualt
        l_r_B =PointSideFromLine(   e1,n1,...
                                    e2,n2,...
                                    Set0.crds(pBi,2),Set0.crds(pBi,1)); % e, n
        
        % total vector change rate per time
%         v = (vFhz./pi) .* atan(lpB./ldF).*l_r_B - (vFhz./pi) .* atan(lpA./ldF).*l_r_A;
        vscr = abs( (vFhz./pi) .* atan(lpB./ldF) - (vFhz./pi) .* atan(lpA./ldF) );  % vector size change rate
        
        if (l_r_B * l_r_A  == 0) 
            error('moveSet3D_lockFault :: unable to handle a píint on the line');
        end
        
        mltpr = NaN; theta = NaN;
        if (l_r_B * l_r_A  == 1)  % on same side of the fault
            % az of the vector is from the closest point to the fault to the farest
            if (lpB > lpA) % same as vector
                dn1 = Set0.VectorsAndVCVs{k,3};
                de1 = Set0.VectorsAndVCVs{k,4};
            else %(lpA >= lpB) -  flip vector
                dn1 = -Set0.VectorsAndVCVs{k,3};
                de1 = -Set0.VectorsAndVCVs{k,4};
            end
            vctrAz = Azimuth(de1 , dn1); %de, dn    %[Az] = Azimuth(de,dn)
            theta = (vctrAz - azF) * 180 / pi;
            
            % vectorssize un-chnaged -> vector is parallel to the fault      //------//
            if (    (theta == 000) || (theta == 180) || (theta == 360) ) 
                mltpr = 0;
            end
            % vector is shrinking   >------<
            if (    ((theta > 000) && (theta < 090))  ||...
                    ((theta > 180) && (theta < 270))       )  
                mltpr = -1;
            end
            % vector is growing   <------>
            if (    ((theta >=090) && (theta < 180))  ||...
                    ((theta >=270) && (theta < 360))       )  
                mltpr = 1;
            end
            
        else % (l_r_B * l_r_A  == -1)  - on opposite sides of the fault
            if (lpB > lpA) % same as vector
                dn1 = Set0.VectorsAndVCVs{k,3};
                de1 = Set0.VectorsAndVCVs{k,4};
            else %(lpA >= lpB) -  flip vector
                dn1 = -Set0.VectorsAndVCVs{k,3};
                de1 = -Set0.VectorsAndVCVs{k,4};
            end
            vctrAz = Azimuth(de1 , dn1); %de, dn    %[Az] = Azimuth(de,dn)
            theta = (vctrAz - azF) * 180 / pi;

            % vectorssize un-chnaged -> vector is parallel to the fault 
            if (    (theta == 000) || (theta == 180) || (theta == 360) ) 
                mltpr = 0;
            end
            % vector is growing   <------>
            if (    ((theta > 000) && (theta < 090))  ||...
                    ((theta > 180) && (theta < 270))       )  
                mltpr = -1;
            end
            % vector is shrinking  >------<
            if (    ((theta > 090) && (theta < 180))  ||...
                    ((theta > 270) && (theta < 360))       )  
                mltpr = 1;
            end
        end
        
        de = vscr*dt*sin( azF ) * mltpr;
        dn = vscr*dt*cos( azF ) * mltpr;
        
        warning('moveSet3D_lockFault :: currently up-veocity is unavailable');
            du = vFv * ( (l_r_B / (lpB*lpB)) - (l_r_A / (lpA*lpA)) );   % similar to the horizntal space
        du = 0; %vFv * dt;
        
        set_new.VectorsAndVCVs{k ,3:5} = set_new.VectorsAndVCVs{k,3:5} + [dn,de,du];
    end
    
    % updaing the entire DataSet
    set_new.Points{:,2:4} = set_new.crds;
    set_new.vctrs         = set_new.VectorsAndVCVs{:,3:5};
    
    % adding noise
    rng('shuffle');  %rng(1234); 
    set_new.vcvs         = set_new.vcvs         + mean([noiseHz,noiseV]) .* rand(size(set_new.vcvs));                   % ------ add the v-Noise to the rnadomization
    set_new.vctrs(:,1:2) = set_new.vctrs(:,1:2) +      noiseHz           .* rand(size(set_new.vctrs,1),2);
    set_new.vctrs(:,3)   = set_new.vctrs(:, 3 ) +      noiseV            .* rand(size(set_new.vctrs,1),1);
    set_new.VectorsAndVCVs{:,3:end} = [set_new.vctrs,set_new.vcvs];
    newSets{t+1} = set_new;
end % for Times

end% end function
%}

%%  2019.03.05  - tested ok +- 1mm

% updated 2019.03.05
newSets = {Set0}; newSets{nTimes} = NaN;

% minE = floor(min(Set0.crds(:,2))) - 3000;% minN = floor(min(Set0.crds(:,1))) - 3000;
% maxE = floor(max(Set0.crds(:,2))) + 3000;% maxN = floor(max(Set0.crds(:,1))) + 3000;
% m = tan(azF);
% e1 = minE; n1 = m*(posF(1) - e1) + posF(2);
% e2 = maxE; n2 = m*(posF(1) - e2) + posF(2);
e1 = posF(1) -10000.*sin(azF);   e2 = posF(1) +10000.*sin(azF);
n1 = posF(2) -10000.*cos(azF);   n2 = posF(2) +10000.*cos(azF);

for t = 1 : (nTimes - 1)
    set_new = newSets{t};                 
    set_new.Time = set_new.Time + dt;
    
    for nms = 1:size(Set0.Points,1)
        for i = 1:size(set_new.Points,1)
            if (strcmpi(set_new.Points{i,1}{1},Set0.Points{nms,1}{1}) == true)
                e = set_new.crds(i,2); n = set_new.crds(i,1);
                l_r = PointSideFromLine(e1,n1,e2,n2,e,n);  % left or right from line
                l  = DistLinePoint(e1,n1,e2,n2,e,n);


                v = vFhz / pi * atan(l/ldF); % 4.4 formula - from PhD tyhsis Lior Shahar
                de = v*dt*sin(azF)       * l_r;
                dn = v*dt*cos(azF)       * l_r;
                du = vFv * dt            * l_r;
            
                %updating crds
                set_new.crds(i,2) =  set_new.crds(i,2) + de;
                set_new.crds(i,1) =  set_new.crds(i,1) + dn;
                set_new.crds(i,3) =  set_new.crds(i,3) + du;
                
                % rectifying vectors 
                % when from Point is moving
                p = string(set_new.VectorsAndVCVs{:,1})==Set0.Points{nms,1}{1};
                set_new.VectorsAndVCVs{ (p == 1) ,3:5} = set_new.VectorsAndVCVs{ (p == 1) ,3:5} - [dn,de,du];
                
                % when  to  Point is moving
                p = string(set_new.VectorsAndVCVs{:,2})==Set0.Points{nms,1}{1};
                set_new.VectorsAndVCVs{ (p == 1) ,3:5} = set_new.VectorsAndVCVs{ (p == 1) ,3:5} + [dn,de,du];
                continue;
            end %if
        end %for i
    end %for nms
    
    % updaing the entire DataSet
    set_new.Points{:,2:4} = set_new.crds;
    set_new.vctrs         = set_new.VectorsAndVCVs{:,3:5};
    % setting up error estimates - adding zero-avarage noise
    set_new.vcvs = rand(size(set_new.vcvs)).*erf(mean(mean(abs(set_new.vcvs))));
    
    % adding noise
    rng('shuffle');
    set_new.vcvs         = set_new.vcvs         + mean([noiseHz,noiseV]) .* rand(size(set_new.vcvs));                   % ------ add the v-Noise to the rnadomization
    set_new.vctrs(:,1:2) = set_new.vctrs(:,1:2) +      noiseHz           .* rand(size(set_new.vctrs,1),2);
    set_new.vctrs(:,3)   = set_new.vctrs(:, 3 ) +      noiseV            .* rand(size(set_new.vctrs,1),1);
    set_new.VectorsAndVCVs{:,3:end} = [set_new.vctrs,set_new.vcvs];
    newSets{t+1} = set_new;
end % for Times

end% end function



%% Old version - big error
%{

function [newSets] = moveSet3D_lockFault(dt,ldF,vFhz,vFv,azF,posF,Set0,noiseHz,noiseV,nTimes)
newSets = {Set0}; newSets{nTimes} = NaN;
for t = 1 : (nTimes - 1)
set_new = Set0;
set_new.Time = set_new.Time + dt;
minE = floor(min(Set0.crds(:,2))) - 3000;% minN = floor(min(Set0.crds(:,1))) - 3000;
maxE = floor(max(Set0.crds(:,2))) + 3000;% maxN = floor(max(Set0.crds(:,1))) + 3000;
m = tan(azF);
e1 = minE; n1 = m*(posF(1) - e1) + posF(2);
e2 = maxE; n2 = m*(posF(1) - e2) + posF(2);



for i = 1:length(set_new.crds)
    e = set_new.crds(i,2); n = set_new.crds(i,1);
    l_r = PointSideFromLine(e1,n1,e2,n2,e,n);
    l  = DistLinePoint(e1,n1,e2,n2,e,n);
    
    v = vFhz / pi * atan(l/ldF); % 4.4 formula - from PhD tyhsis Lior Shahar
    de = v*dt*sin(azF)       * l_r;
    dn = v*dt*cos(azF)       * l_r;
    du = vFv * dt            * l_r;
    
    set_new.crds(i,2) =  e + de;
    set_new.crds(i,1) =  n + dn;
    set_new.crds(i,3) =  set_new.crds(i,3) + du;
end
%updating
rng(1234);
set_new.Points{:,2:4} = set_new.crds;
set_new.vcvs = set_new.vcvs + (noiseHz-0)*rand(3,6);  % ------ add the v-Noise to the rnadomization

set_new.VectorsAndVCVs{:,6:end} = set_new.vcvs;
for i = 1:size(set_new.VectorsAndVCVs,1)
    p1 = find(string(set_new.VectorsAndVCVs{i,1})==set_new.Points{:,1},1);
    p2 = find(string(set_new.VectorsAndVCVs{i,2})==set_new.Points{:,1},1);
    dneu = set_new.crds(p2,:) - set_new.crds(p1,:);
    set_new.VectorsAndVCVs(i,3:5) = {dneu(1),dneu(2),dneu(3)};
end
set_new.vctrs = set_new.VectorsAndVCVs{:,3:5};

newSets{t+1} = set_new;
end % for Times
end% end function
%}
