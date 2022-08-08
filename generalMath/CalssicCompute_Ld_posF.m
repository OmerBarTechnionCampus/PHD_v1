function [ld,pose,posn, sigpost2] = CalssicCompute_Ld_posF(Set0,Set1,v0,az0,pose0,posn0,ld0)
%% this function will compute the locked-fault model parameters using Adjstment

%% computing vectors speeds [mm/year]
dn = (Set1.vctrs(:,1)-Set0.vctrs(:,1));
de = (Set1.vctrs(:,2)-Set0.vctrs(:,2)); 
vs = sqrt((dn.*dn) +...      
          (de.*de)      ) / (Set1.Time - Set0.Time);  
% az = round(Azimuth(Set1.vctrs(:,2),Set1.vctrs(:,1)),6);
% vs = sqrt((Set1.vctrs(:,1)-Set0.vctrs(:,1)).*(Set1.vctrs(:,1)-Set0.vctrs(:,1)) +...      %dn2
%     (Set1.vctrs(:,2)-Set0.vctrs(:,2)).*(Set1.vctrs(:,2)-Set0.vctrs(:,2))      );   %de2
vs = vs * 1000; %mm per year
i1 = (abs(dn) > de) .* (dn < 0); % where dn < 0
i2 = (abs(de) > dn) .* (de < 0); % where de < 0

i = find((i1+i2) >0);  % where (de < 0) or (dn < 0)
vs(i) = -vs(i);  % minus sign for shortening vectors

%% preparing matrices for adjustment
lb = reshape(vs,[numel(vs),1]);
X = [pose0;posn0;ld0];
dX = [99999;99999;99999];

thPos = 100; %[m]
thLD = 1000; % [m]
while (norm(dX(1:2)) > thPos) && (abs(dX(3)) > thLD)
    pose0 = X(1); posn0 = X(2); ld0 = X(3);
    A = zeros(numel(vs),3);
    L0 = zeros(numel(vs),1);

    % compute A, L0
    for k = 1:numel(vs)
        pAi = find(string(Set0.VectorsAndVCVs{k,1})==Set0.Points{:,1}); % Vector's FromPoint
        pBi = find(string(Set0.VectorsAndVCVs{k,2})==Set0.Points{:,1}); % Vector's   ToPoint
        
        nf = Set0.Points{pAi,2}; ef = Set0.Points{pAi,3};
        nt = Set0.Points{pBi,2}; et = Set0.Points{pBi,3};
        % Distance of FromPoint from simulated fault
        %round
        lpA = (DistLinePoint(       pose0-10.*sin(az0)      ,posn0-10.*cos(az0),...
                                    pose0+10.*sin(az0)      ,posn0+10.*cos(az0),...
                                    ef ,nf));

        l_r_A =PointSideFromLine(   pose0-10.*sin(az0),posn0-10.*cos(az0),...
                                    pose0+10.*sin(az0),posn0+10.*cos(az0),...
                                    ef,nf);

        lpB = (DistLinePoint(       pose0-10.*sin(az0)      ,posn0-10.*cos(az0),...
                                    pose0+10.*sin(az0)      ,posn0+10.*cos(az0),...
                                    et ,nt));

        l_r_B =PointSideFromLine(   pose0-10.*sin(az0),posn0-10.*cos(az0),...
                                    pose0+10.*sin(az0),posn0+10.*cos(az0),...
                                    et,nt);


        % vector (edge) change rate (velocity) due to end-points (nodes) movement....
        % [locked fault model]
        v_v = (v0./pi) .* atan(lpB./ld0).*l_r_B  -  (v0./pi) .* atan(lpA./ld0).*l_r_A ;
        
        % deriviatives using symbolic matlab....
        diff_V_V_pose0 = (40*v0*atan(abs((posn0 - 10*cos(az0))*(pose0 + 10*sin(az0)) - (posn0 + 10*cos(az0))*(pose0 - 10*sin(az0)) + 20*ef*cos(az0) - 20*nf*sin(az0))/(20*ld0*(cos(az0)^2 + sin(az0)^2)^(1/2)))*dirac(20*cos(az0)*(ef - pose0 + 10*sin(az0)) - 20*sin(az0)*(nf - posn0 + 10*cos(az0)))*cos(az0))/pi - (40*v0*atan(abs((posn0 - 10*cos(az0))*(pose0 + 10*sin(az0)) - (posn0 + 10*cos(az0))*(pose0 - 10*sin(az0)) + 20*et*cos(az0) - 20*nt*sin(az0))/(20*ld0*(cos(az0)^2 + sin(az0)^2)^(1/2)))*dirac(20*cos(az0)*(et - pose0 + 10*sin(az0)) - 20*sin(az0)*(nt - posn0 + 10*cos(az0)))*cos(az0))/pi + (v0*sign(20*cos(az0)*(ef - pose0 + 10*sin(az0)) - 20*sin(az0)*(nf - posn0 + 10*cos(az0)))*cos(az0)*sign((posn0 - 10*cos(az0))*(pose0 + 10*sin(az0)) - (posn0 + 10*cos(az0))*(pose0 - 10*sin(az0)) + 20*ef*cos(az0) - 20*nf*sin(az0)))/(ld0*pi*(cos(az0)^2 + sin(az0)^2)^(1/2)*(abs((posn0 - 10*cos(az0))*(pose0 + 10*sin(az0)) - (posn0 + 10*cos(az0))*(pose0 - 10*sin(az0)) + 20*ef*cos(az0) - 20*nf*sin(az0))^2/(400*ld0^2*(cos(az0)^2 + sin(az0)^2)) + 1)) - (v0*sign(20*cos(az0)*(et - pose0 + 10*sin(az0)) - 20*sin(az0)*(nt - posn0 + 10*cos(az0)))*cos(az0)*sign((posn0 - 10*cos(az0))*(pose0 + 10*sin(az0)) - (posn0 + 10*cos(az0))*(pose0 - 10*sin(az0)) + 20*et*cos(az0) - 20*nt*sin(az0)))/(ld0*pi*(cos(az0)^2 + sin(az0)^2)^(1/2)*(abs((posn0 - 10*cos(az0))*(pose0 + 10*sin(az0)) - (posn0 + 10*cos(az0))*(pose0 - 10*sin(az0)) + 20*et*cos(az0) - 20*nt*sin(az0))^2/(400*ld0^2*(cos(az0)^2 + sin(az0)^2)) + 1));
        diff_V_V_posn0 = (40*v0*atan(abs((posn0 - 10*cos(az0))*(pose0 + 10*sin(az0)) - (posn0 + 10*cos(az0))*(pose0 - 10*sin(az0)) + 20*et*cos(az0) - 20*nt*sin(az0))/(20*ld0*(cos(az0)^2 + sin(az0)^2)^(1/2)))*dirac(20*cos(az0)*(et - pose0 + 10*sin(az0)) - 20*sin(az0)*(nt - posn0 + 10*cos(az0)))*sin(az0))/pi - (40*v0*atan(abs((posn0 - 10*cos(az0))*(pose0 + 10*sin(az0)) - (posn0 + 10*cos(az0))*(pose0 - 10*sin(az0)) + 20*ef*cos(az0) - 20*nf*sin(az0))/(20*ld0*(cos(az0)^2 + sin(az0)^2)^(1/2)))*dirac(20*cos(az0)*(ef - pose0 + 10*sin(az0)) - 20*sin(az0)*(nf - posn0 + 10*cos(az0)))*sin(az0))/pi - (v0*sign(20*cos(az0)*(ef - pose0 + 10*sin(az0)) - 20*sin(az0)*(nf - posn0 + 10*cos(az0)))*sin(az0)*sign((posn0 - 10*cos(az0))*(pose0 + 10*sin(az0)) - (posn0 + 10*cos(az0))*(pose0 - 10*sin(az0)) + 20*ef*cos(az0) - 20*nf*sin(az0)))/(ld0*pi*(cos(az0)^2 + sin(az0)^2)^(1/2)*(abs((posn0 - 10*cos(az0))*(pose0 + 10*sin(az0)) - (posn0 + 10*cos(az0))*(pose0 - 10*sin(az0)) + 20*ef*cos(az0) - 20*nf*sin(az0))^2/(400*ld0^2*(cos(az0)^2 + sin(az0)^2)) + 1)) + (v0*sign(20*cos(az0)*(et - pose0 + 10*sin(az0)) - 20*sin(az0)*(nt - posn0 + 10*cos(az0)))*sin(az0)*sign((posn0 - 10*cos(az0))*(pose0 + 10*sin(az0)) - (posn0 + 10*cos(az0))*(pose0 - 10*sin(az0)) + 20*et*cos(az0) - 20*nt*sin(az0)))/(ld0*pi*(cos(az0)^2 + sin(az0)^2)^(1/2)*(abs((posn0 - 10*cos(az0))*(pose0 + 10*sin(az0)) - (posn0 + 10*cos(az0))*(pose0 - 10*sin(az0)) + 20*et*cos(az0) - 20*nt*sin(az0))^2/(400*ld0^2*(cos(az0)^2 + sin(az0)^2)) + 1));
        diff_V_V_ld0   = (v0*sign(20*cos(az0)*(ef - pose0 + 10*sin(az0)) - 20*sin(az0)*(nf - posn0 + 10*cos(az0)))*abs((posn0 - 10*cos(az0))*(pose0 + 10*sin(az0)) - (posn0 + 10*cos(az0))*(pose0 - 10*sin(az0)) + 20*ef*cos(az0) - 20*nf*sin(az0)))/(20*ld0^2*pi*(cos(az0)^2 + sin(az0)^2)^(1/2)*(abs((posn0 - 10*cos(az0))*(pose0 + 10*sin(az0)) - (posn0 + 10*cos(az0))*(pose0 - 10*sin(az0)) + 20*ef*cos(az0) - 20*nf*sin(az0))^2/(400*ld0^2*(cos(az0)^2 + sin(az0)^2)) + 1)) - (v0*sign(20*cos(az0)*(et - pose0 + 10*sin(az0)) - 20*sin(az0)*(nt - posn0 + 10*cos(az0)))*abs((posn0 - 10*cos(az0))*(pose0 + 10*sin(az0)) - (posn0 + 10*cos(az0))*(pose0 - 10*sin(az0)) + 20*et*cos(az0) - 20*nt*sin(az0)))/(20*ld0^2*pi*(cos(az0)^2 + sin(az0)^2)^(1/2)*(abs((posn0 - 10*cos(az0))*(pose0 + 10*sin(az0)) - (posn0 + 10*cos(az0))*(pose0 - 10*sin(az0)) + 20*et*cos(az0) - 20*nt*sin(az0))^2/(400*ld0^2*(cos(az0)^2 + sin(az0)^2)) + 1));
        
        A(k,:) = [diff_V_V_pose0, diff_V_V_posn0, diff_V_V_ld0];
        L0(k,1) = v_v;
    end %for k = 1:numel(vs)
    
    L = L0-lb;
    dX = (A' * A) \ (A' * L);
    X = X + dX;
    v = A * dX - L;
end

sigpost2 = v'*v / ( size(A,1) - size(A,2) );
pose = X(1); posn = X(2); ld = X(3);
end % function


%{
syms Vhz Az posE posN Ld nf ef nt et



lpA = (DistLinePoint(  posE-10.*sin(Az)      ,posN-10.*cos(Az),...
    posE+10.*sin(Az)      ,posN+10.*cos(Az),...
    ef ,nf));

l_r_A =PointSideFromLine(   posE-10.*sin(Az),posN-10.*cos(Az),...
    posE+10.*sin(Az),posN+10.*cos(Az),...
    ef,nf);

lpB = (DistLinePoint(  posE-10.*sin(Az)      ,posN-10.*cos(Az),...
    posE+10.*sin(Az)      ,posN+10.*cos(Az),...
    et ,nt));

l_r_B =PointSideFromLine(   posE-10.*sin(Az),posN-10.*cos(Az),...
    posE+10.*sin(Az),posN+10.*cos(Az),...
    et,nt);


% vector (edge) change rate (velocity) due to end-points (nodes) movement....
% [locked fault model]
v_v = (Vhz./pi) .* atan(lpB./Ld).*l_r_B  -  (Vhz./pi) .* atan(lpA./Ld).*l_r_A ;

diff(v_v, posE)
 
ans =
 
(40*Vhz*atan(abs((posN - 10*cos(Az))*(posE + 10*sin(Az)) - (posN + 10*cos(Az))*(posE - 10*sin(Az)) + 20*ef*cos(Az) - 20*nf*sin(Az))/(20*Ld*(cos(Az)^2 + sin(Az)^2)^(1/2)))*dirac(20*cos(Az)*(ef - posE + 10*sin(Az)) - 20*sin(Az)*(nf - posN + 10*cos(Az)))*cos(Az))/pi - (40*Vhz*atan(abs((posN - 10*cos(Az))*(posE + 10*sin(Az)) - (posN + 10*cos(Az))*(posE - 10*sin(Az)) + 20*et*cos(Az) - 20*nt*sin(Az))/(20*Ld*(cos(Az)^2 + sin(Az)^2)^(1/2)))*dirac(20*cos(Az)*(et - posE + 10*sin(Az)) - 20*sin(Az)*(nt - posN + 10*cos(Az)))*cos(Az))/pi + (Vhz*sign(20*cos(Az)*(ef - posE + 10*sin(Az)) - 20*sin(Az)*(nf - posN + 10*cos(Az)))*cos(Az)*sign((posN - 10*cos(Az))*(posE + 10*sin(Az)) - (posN + 10*cos(Az))*(posE - 10*sin(Az)) + 20*ef*cos(Az) - 20*nf*sin(Az)))/(Ld*pi*(cos(Az)^2 + sin(Az)^2)^(1/2)*(abs((posN - 10*cos(Az))*(posE + 10*sin(Az)) - (posN + 10*cos(Az))*(posE - 10*sin(Az)) + 20*ef*cos(Az) - 20*nf*sin(Az))^2/(400*Ld^2*(cos(Az)^2 + sin(Az)^2)) + 1)) - (Vhz*sign(20*cos(Az)*(et - posE + 10*sin(Az)) - 20*sin(Az)*(nt - posN + 10*cos(Az)))*cos(Az)*sign((posN - 10*cos(Az))*(posE + 10*sin(Az)) - (posN + 10*cos(Az))*(posE - 10*sin(Az)) + 20*et*cos(Az) - 20*nt*sin(Az)))/(Ld*pi*(cos(Az)^2 + sin(Az)^2)^(1/2)*(abs((posN - 10*cos(Az))*(posE + 10*sin(Az)) - (posN + 10*cos(Az))*(posE - 10*sin(Az)) + 20*et*cos(Az) - 20*nt*sin(Az))^2/(400*Ld^2*(cos(Az)^2 + sin(Az)^2)) + 1))
 
diff(v_v, posN)
 
ans =
 
(40*Vhz*atan(abs((posN - 10*cos(Az))*(posE + 10*sin(Az)) - (posN + 10*cos(Az))*(posE - 10*sin(Az)) + 20*et*cos(Az) - 20*nt*sin(Az))/(20*Ld*(cos(Az)^2 + sin(Az)^2)^(1/2)))*dirac(20*cos(Az)*(et - posE + 10*sin(Az)) - 20*sin(Az)*(nt - posN + 10*cos(Az)))*sin(Az))/pi - (40*Vhz*atan(abs((posN - 10*cos(Az))*(posE + 10*sin(Az)) - (posN + 10*cos(Az))*(posE - 10*sin(Az)) + 20*ef*cos(Az) - 20*nf*sin(Az))/(20*Ld*(cos(Az)^2 + sin(Az)^2)^(1/2)))*dirac(20*cos(Az)*(ef - posE + 10*sin(Az)) - 20*sin(Az)*(nf - posN + 10*cos(Az)))*sin(Az))/pi - (Vhz*sign(20*cos(Az)*(ef - posE + 10*sin(Az)) - 20*sin(Az)*(nf - posN + 10*cos(Az)))*sin(Az)*sign((posN - 10*cos(Az))*(posE + 10*sin(Az)) - (posN + 10*cos(Az))*(posE - 10*sin(Az)) + 20*ef*cos(Az) - 20*nf*sin(Az)))/(Ld*pi*(cos(Az)^2 + sin(Az)^2)^(1/2)*(abs((posN - 10*cos(Az))*(posE + 10*sin(Az)) - (posN + 10*cos(Az))*(posE - 10*sin(Az)) + 20*ef*cos(Az) - 20*nf*sin(Az))^2/(400*Ld^2*(cos(Az)^2 + sin(Az)^2)) + 1)) + (Vhz*sign(20*cos(Az)*(et - posE + 10*sin(Az)) - 20*sin(Az)*(nt - posN + 10*cos(Az)))*sin(Az)*sign((posN - 10*cos(Az))*(posE + 10*sin(Az)) - (posN + 10*cos(Az))*(posE - 10*sin(Az)) + 20*et*cos(Az) - 20*nt*sin(Az)))/(Ld*pi*(cos(Az)^2 + sin(Az)^2)^(1/2)*(abs((posN - 10*cos(Az))*(posE + 10*sin(Az)) - (posN + 10*cos(Az))*(posE - 10*sin(Az)) + 20*et*cos(Az) - 20*nt*sin(Az))^2/(400*Ld^2*(cos(Az)^2 + sin(Az)^2)) + 1))
 
diff(v_v, Ld)
 
ans =
 
(Vhz*sign(20*cos(Az)*(ef - posE + 10*sin(Az)) - 20*sin(Az)*(nf - posN + 10*cos(Az)))*abs((posN - 10*cos(Az))*(posE + 10*sin(Az)) - (posN + 10*cos(Az))*(posE - 10*sin(Az)) + 20*ef*cos(Az) - 20*nf*sin(Az)))/(20*Ld^2*pi*(cos(Az)^2 + sin(Az)^2)^(1/2)*(abs((posN - 10*cos(Az))*(posE + 10*sin(Az)) - (posN + 10*cos(Az))*(posE - 10*sin(Az)) + 20*ef*cos(Az) - 20*nf*sin(Az))^2/(400*Ld^2*(cos(Az)^2 + sin(Az)^2)) + 1)) - (Vhz*sign(20*cos(Az)*(et - posE + 10*sin(Az)) - 20*sin(Az)*(nt - posN + 10*cos(Az)))*abs((posN - 10*cos(Az))*(posE + 10*sin(Az)) - (posN + 10*cos(Az))*(posE - 10*sin(Az)) + 20*et*cos(Az) - 20*nt*sin(Az)))/(20*Ld^2*pi*(cos(Az)^2 + sin(Az)^2)^(1/2)*(abs((posN - 10*cos(Az))*(posE + 10*sin(Az)) - (posN + 10*cos(Az))*(posE - 10*sin(Az)) + 20*et*cos(Az) - 20*nt*sin(Az))^2/(400*Ld^2*(cos(Az)^2 + sin(Az)^2)) + 1))
 


%}