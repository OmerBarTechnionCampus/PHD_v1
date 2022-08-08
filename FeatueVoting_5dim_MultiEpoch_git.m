% Feature Voting for locked fault model
%
%
% Omer Bar
close all; clc;
pack; clear all; pack; 
clc;
%% perps
% addpath ('D:\GoogleDrive\_Coding\Matlab');
% addpath ('D:\GoogleDrive\_Coding\Matlab\Geodesy');
% addpath ('D:\GoogleDrive\_Coding\Matlab\IO');
% AddPaths('D:\GoogleDrive\_Coding\Matlab\IO');
% addpath ('D:\GoogleDrive\_Coding\Matlab\GUI');
% addpath ('D:\GoogleDrive\_Coding\Matlab\Math');
% addpath ('D:\GoogleDrive\_Coding\Matlab\GPS\')
% AddPaths('D:\GoogleDrive\_Coding\Matlab\GPS\');
% AddPaths('D:\GoogleDrive\PHD\coding\'); AddPaths(cd);
% clear all; format long g; close all; fclose all;  clc;

% CSAR  IG05/12  E:189919.643, N:710505.891?	

mPath2outputFolder = '\\coding\save_run_files\';
curRun = datestr(datetime('now'),'YYYY_mm_dd_HH_MM_SS');
mkdir([mPath2outputFolder,curRun]);
diary([mPath2outputFolder,curRun,'\diary.txt']);
%% verbose options
isVerbose = true;
isShowFigs = false; %#ok 
isWraning = false;  %#ok


        trueAz = NaN; %[deg]
        trueHz_Vel = NaN; %[mm/year]
        trueLd = NaN; %[deg]
        truePosE = NaN;
        truePosN = NaN;
        trueDeltaTime = NaN; %[years]
         
if (isVerbose == false); clc; end;

%% work mode
% Single Fault Mode
strMode = 'SingleFault';
% % Multi Fault Mode
% strMode = 'MultiFault';

% Fault Types
% strFaultType = 'Slip';
strFaultType = 'Locked';

% Pertrubation Method
% prtbMethod = 'SAMPLE_STD';
prtbMethod = 'tiny';

% mInputMethod = 'Synthetic';
mInputMethod = 'Bernese';
% mInputMethod = 'GiladTests';

fprintf('Mode = %s\nFaultType = %s\nPertrubation Method = %s\nInput = %s\n',strMode,strFaultType,prtbMethod,mInputMethod);
switch (mInputMethod)
    case 'Synthetic'
        % --- for data in GEOLAB format
%         in_paths = {'\\coding\SyntheticData\G0_data4_triangle_DB\2015_032_geolab.txt',...
%                     '\\coding\SyntheticData\G0_data4_triangle_DB\2016_001_geolab.txt',...
%                     '\\coding\SyntheticData\G0_data4_triangle_DB\2017_001_geolab.txt',...
%                     '\\coding\SyntheticData\G0_data4_triangle_DB\2018_001_geolab.txt'     };

        in_paths = {'\\coding\SyntheticData\G0_data4_triangle_DB\2015_032_geolab.txt',};
        [Sets] = ReadDataSetsFromGeolbFiles(in_paths);
        
        
        
        t0Set = Sets{1};
        
        trueAz = 0.3; %[deg]
        trueHz_Vel = 8.2; %[mm/year]
        trueLd = 12015; %[m]
        truePosE = 62281; % easting between MRAV - KATZ - DSF
        truePosN = 50000;
                
        
        trueDeltaTime = 1.0; %[years]
        % ----- noise ----
        
        % - vertical noise
        noiseV  = 0; % [mm]
        
%         acc = set.vcvs(:,end).*set.vcvs(:,end);
%         acc = sqrt(acc);
%         maxVAcc = max(acc); %[meters]
%         noiseV = maxVAcc * 1000;%[mili-meters]
        
        % - horizontal noise - 
        % handeling true-measured noise
%         acc = t0Set.vcvs(:,1).*t0Set.vcvs(:,1) + t0Set.vcvs(:,2).*t0Set.vcvs(:,2) + 2.*t0Set.vcvs(:,3);
%         acc = sqrt(sqrt(acc));
%         maxHzAcc = max(acc); %[meters]
%             % to set the proper error estimate between a scientific program and a commercial program;
%             maxHzAcc = maxHzAcc / 10;
%         noiseHz = maxHzAcc * 1000; %[mili-meters]

        noiseHz = trueHz_Vel * 0.08; % [mm] - max presentage working on synthetic data - 10%
        % noiseHz = sqrt(  max(abs(t0Set.vcvs(:)))  ) * 1 * 1000; % [mm]
        
%         noiseHz = 0; % [m]
%         noiseHz = 0.01; % [mm] = 0.00001m
        
        
        
        
        mAddedEpochs = 2;

        switch (strFaultType)
            case 'Locked'
                    % SiteMove_Names = {'HRMN' , 'KROM',  'KATZ', 'CSAR' , 'BSHM' , 'KABR' , 'NZRT' , 'MRAV'};
                    % moveSet3D_lockFault(dt,ldF,vFhz,vFv,azF,posF,Set0,noiseHz,noiseV,nTimes)
                fprintf('\n**** This is For Locked Fualt... ****\n');
                fprintf('\nTrue Locked-Fault parameters\nAz = %.3f, vHz = %.1fmm, Depth = %.1f\n',trueAz,trueHz_Vel,trueLd);
                fprintf('TruePos = E:%.1f; N:%.1f\n',truePosE,truePosN);
                fprintf('Noise = Hz:%.1f; Up:%.1f\n',noiseHz,noiseV);
                synthSets = moveSet3D_lockFault(trueDeltaTime,trueLd,...           % deltaYear = 1.0 ; FaultDepth=10000;
                    trueHz_Vel/1000,0,trueAz.*pi/180,...       % Hz_vel = 0.006; vt_vel = 0;  Az = 0;
                    [truePosE,truePosN],t0Set,noiseHz/1000,noiseV/1000,mAddedEpochs+1); %[East = 65000, North = 50000] , on original set + single epoch
                t1Set = synthSets{end};
            case 'Slip'
                fprintf('\n**** This is For Slip Fualt... ****\n');
                fprintf('\nTrue Slip-Fault parameters\nAz = %.3f, vHz = %.1fmm\n',trueAz,trueHz_Vel);
                fprintf('TruePos = E:%.1f; N:%.1f\n',truePosE,truePosN); 
                fprintf('Noise = Hz:%.1f; Up:%.1f\n',noiseHz,noiseV);
                synthSets = moveSet3D_slipFault(trueDeltaTime,...           % deltaYear = 1.0 ;
                    trueHz_Vel/1000,0,trueAz.*pi/180,...       % Hz_vel = 0.006; vt_vel = 0;  Az = 0;
                    [truePosE,truePosN],t0Set,noiseHz/1000,noiseV/1000,mAddedEpochs+1); %[East = 65000, North = 50000] , on original set + single epoch
                t1Set = synthSets{end};
        end %switch
        Sets = synthSets;
        fprintf('==============Synthetic Data Created =========================\n');
    case 'Bernese'
        
        % true data on Bernese format
        fprintf('.\n'); fprintf('..\n'); fprintf('...\n');
        fprintf('...loading true data.....\n');
        

        %--- for data from BERNESE CAMPAINGS
        % reads daily data and merge it manually - less accurate
        years = 2010:2021; years =years-2000;
        Sets{numel(years)} = {};
        path2BSWsol = '\\coding\SyntheticData\BerneseCampaigns\';
        for yy = 1:numel(years)
            dCRD = dir([path2BSWsol,'RED',num2str(years(yy)),'*.CRD']);  %dCRD = dir([path2BSWsol,'RED11*.CRD']);
            dCOV = dir([path2BSWsol,'RED_',num2str(years(yy)),'*.COV']); %dCOV = dir([path2BSWsol,'RED_11*.COV']);

            in_paths_CRD{numel(dCRD)} = {};in_paths_COV{numel(dCRD)} = {};
            for k = 1:numel(dCRD)
                in_paths_CRD{k} = [dCRD(k).folder,'\',dCRD(k).name];
                in_paths_COV{k} = [dCOV(k).folder,'\',dCOV(k).name];
            end

            % in_paths_CRD = {'d:\GoogleDrive\PHD\coding\SyntheticData\BerneseCampaigns\RED10002.CRD', ...
            %                 'd:\GoogleDrive\PHD\coding\SyntheticData\BerneseCampaigns\RED21001.CRD'};
            % in_paths_COV = {'d:\GoogleDrive\PHD\coding\SyntheticData\BerneseCampaigns\RED_10002.COV', ...
            %                 'd:\GoogleDrive\PHD\coding\SyntheticData\BerneseCampaigns\RED_21001.COV'};
            [SetsYY] = ReadDataSetsFromBerneseFiles(in_paths_CRD,in_paths_COV);
            Sets{yy} = MergeSets(SetsYY,true);
        end
       
        
        % removing sites from data
        for k = 1:length(Sets)
            [set_new] = RemoveSites(Sets{k},{'ADIS','ANKR','BHR2','ISBA','MDVJ','NICO','TEHN','YIBL','ZECK'});
            [set_new] = RemoveSites(set_new,{'ALON','DRAG','DSEA','ELAT','JSLM','KLHV','NIZN','NRIF','RAMO','SPIR','SLOM','KSLM','TEL0','YOSH','YRCM','TELA'});
            Sets{k} = set_new;
        end
        % removing HRMN  
        for i = 1: numel(Sets)
            ts = Sets{i};
            ts = RemoveSite(ts,'HRMN');
            
            Sets{i} = ts;
        end % for

        Sets = {Sets{1} , Sets{5}, Sets{6}, Sets{10} };
        fprintf('==============True Data Loaded =========================\n');
    case 'GiladTests'
        % % % % % ---------- gilad examples ---------- % % % % % 
        %
        % % fprintf('Gilad Test Slip Fault');
        % load('\\test_runs\2019_10\Gilad_fixed_case1.mat');
        % t0Set = t0SetG; t1Set = t1SetG;
        %
        % fprintf('Gilad Test Locked Fault\n');
        load('\\test_runs\2019_10\Gilad_fixed_case2.mat');
        t0Set = t0SetG; t1Set = t1SetG;
        Sets = {t0Set , t1Set};
        %
        % % % % %fixing the error_estimates(once, then re-saved to mat files)
        % % % % rng('shuffle'); % setting randomization
        % % % % errScale  = randn([size(t1Set.VectorsAndVCVs,1), 3])./10; % creating random numbers for each vector
        % % % % errScale2 = randn([size(t1Set.VectorsAndVCVs,1), 3])./10; % creating random numbers for each vector
        % % % %     full_err = [errScale(:,1),errScale2(:,1:2),errScale(:,2),errScale2(:,3),errScale(:,3)]*mean(mean(t1Set.vcvs));% % % % % ---------- gilad examples ---------- % % % % % 
        % % % % t1Set.VectorsAndVCVs{:,6:end} = full_err;
end %switch


%% Arranging Vectors 
%  to be all in same order
%  and have dTime data embeded for each vector

 [t0Set, t1Set] = ArrangeMultiEpochSets(Sets); % from data read from sets
%  [t0Set, t1Set] = ArrangeMultiEpochSets(synthSets); % from synthetic data
 
 InitialVectorSet = size(t0Set.vcvs,1);

%%
% rouding the true data in....
SigDig = 4;
fprintf('rounding the original data to %d significant digits\n',SigDig);
t0Set.crds  = round(t0Set.crds,1);
t0Set.vctrs = round(t0Set.vctrs,SigDig);
t0Set.Points{:,2:4} = round(t0Set.Points{:,2:4},1);
t0Set.VectorsAndVCVs{:,3:5} = round(t0Set.VectorsAndVCVs{:,3:5},SigDig);

t1Set.crds  = round(t1Set.crds,1);
t1Set.vctrs = round(t1Set.vctrs,SigDig);
t1Set.Points{:,2:4} = round(t1Set.Points{:,2:4},1);
t1Set.VectorsAndVCVs{:,3:5} = round(t1Set.VectorsAndVCVs{:,3:5},SigDig);


%% adding pretrubed data
% adding the noise of the measurement itself
% [t0Set] = AddPertrubedVectors(t0Set,1000);
% [t1Set] = AddPertrubedVectors(t1Set,1000);
% [t0Set1] = AddPertrubedVectors(t0Set,100);
% [t1Set1] = AddPertrubedVectors(t1Set,100);

% adding zeros mean random noise to the vector sizes
orig_t0Set = t0Set;
orig_t1Set = t1Set;

switch (strFaultType)
    case 'Locked'
        mltplr = 100; %1000
    case 'Slip'
        mltplr = 20;
end % switch
rng('shuffle'); % re-setting randomization
fprintf('adding %d pertrubed vectors (mupltipler of %d)\n',mltplr*size(t0Set.VectorsAndVCVs,1),mltplr);
[t0Set] = AddPertrubedVectors_Random_MultiEpoch(t0Set,mltplr,prtbMethod);  % floor(size(t0Set.VectorsAndVCVs,1)/2)
[t1Set] = AddPertrubedVectors_Random_MultiEpoch(t1Set,mltplr,prtbMethod);  % floor(size(t1Set.VectorsAndVCVs,1)/2)

% % % % dual pertrubation
% % % prtbMethod = 'SAMPLE_STD';  % using vcv's
% % % [t0Set] = AddPertrubedVectors_Random_MultiEpoch(t0Set,mltplr*0.05,prtbMethod);  % floor(size(t0Set.VectorsAndVCVs,1)/2)
% % % [t1Set] = AddPertrubedVectors_Random_MultiEpoch(t1Set,mltplr*0.05,prtbMethod);  % floor(size(t1Set.VectorsAndVCVs,1)/2)
% % % 
% % % prtbMethod = 'tiny';    % using pertrubation
% % % [t0Set] = AddPertrubedVectors_Random_MultiEpoch(t0Set,mltplr*1.00,prtbMethod);  % floor(size(t0Set.VectorsAndVCVs,1)/2)
% % % [t1Set] = AddPertrubedVectors_Random_MultiEpoch(t1Set,mltplr*1.00,prtbMethod);  % floor(size(t1Set.VectorsAndVCVs,1)/2)


% % % % % t0Set(end+1:end+size(t0Set1,1)-size(t0Set,1),:) = t0Set1(size(t0Set,1)+1:end,:);
% % % % % t0Set(end+1:end+size(t0Set2,1)-size(t0Set,1),:) = t0Set2(size(t0Set,1)+1:end,:);
% % % % % t1Set(end+1:end+size(t1Set1,1)-size(t1Set,1),:) = t1Set1(size(t1Set,1)+1:end,:);
% % % % % t1Set(end+1:end+size(t1Set2,1)-size(t1Set,1),:) = t1Set2(size(t1Set,1)+1:end,:);


% % % % rouding....
% % % SigDig = 7;
% % % fprintf('rounding the data to %d significant digits\n',SigDig);
% % % t0Set.crds  = round(t0Set.crds,1);
% % % t0Set.vctrs = round(t0Set.vctrs,SigDig);
% % % t0Set.Points{:,2:4} = round(t0Set.Points{:,2:4},1);
% % % t0Set.VectorsAndVCVs{:,3:5} = round(t0Set.VectorsAndVCVs{:,3:5},SigDig);
% % % 
% % % t1Set.crds  = round(t1Set.crds,1);
% % % t1Set.vctrs = round(t1Set.vctrs,SigDig);
% % % t1Set.Points{:,2:4} = round(t1Set.Points{2,4},1);
% % % t1Set.VectorsAndVCVs{:,3:5} = round(t1Set.VectorsAndVCVs{:,3:5},SigDig);


%% ---------------------- normalization ----------------------- %
fprintf('normalizing data\n');
% simple divid due to fact that one of the sites is "0,0,0"
%
nrmFctr = max( [reshape(abs(t0Set.crds),[numel(t0Set.crds),1]);reshape(abs(t0Set.vctrs),[numel(t0Set.vctrs),1]);...
                reshape(abs(t1Set.crds),[numel(t1Set.crds),1]);reshape(abs(t1Set.vctrs),[numel(t1Set.vctrs),1])]    );
r2 = t0Set.vctrs(:,1).* t0Set.vctrs(:,1) + t0Set.vctrs(:,2).* t0Set.vctrs(:,2) + t0Set.vctrs(:,3).* t0Set.vctrs(:,3) ;

nrmFctr = ceil( max(nrmFctr, sqrt(max(r2))) );

% nrmFctr = 1; % overide to test threshold computations

fprintf('normalization Factor = %f\n',nrmFctr);
% nrmFctr = 2 * nrmFctr;            
% % % % % % % % for not normalizing         
% % % % % % % nrmFctr = 1; % for not normalizing

t0Set.crds = t0Set.crds./nrmFctr;
t0Set.vctrs = t0Set.vctrs./nrmFctr;   
t0Set.vcvs = t0Set.vcvs./nrmFctr;
t0Set.Points{:,2:4} = t0Set.Points{:,2:4}./nrmFctr;
t0Set.VectorsAndVCVs{:,3:end} = t0Set.VectorsAndVCVs{:,3:end}./nrmFctr;

t1Set.crds = t1Set.crds./nrmFctr;
t1Set.vctrs = t1Set.vctrs./nrmFctr;
t1Set.vcvs = t1Set.vcvs./nrmFctr;
t1Set.Points{:,2:4} = t1Set.Points{:,2:4}./nrmFctr;
t1Set.VectorsAndVCVs{:,3:end} = t1Set.VectorsAndVCVs{:,3:end}./nrmFctr;

%----- keeping aside the original set --------
t0Set_0.crds = t0Set.crds;
t0Set_0.Points = t0Set.Points;
t0Set_0.vctrs = t0Set.vctrs(1:InitialVectorSet,:);
t0Set_0.vcvs = t0Set.vcvs(1:InitialVectorSet,:);
t0Set_0.VectorsAndVCVs = t0Set.VectorsAndVCVs(1:InitialVectorSet,:);
t0Set_0.Time = t0Set.Time(1:InitialVectorSet,:);

t1Set_0.crds = t1Set.crds;
t1Set_0.Points = t1Set.Points;
t1Set_0.vctrs = t1Set.vctrs(1:InitialVectorSet,:);
t1Set_0.vcvs = t1Set.vcvs(1:InitialVectorSet,:);
t1Set_0.VectorsAndVCVs = t1Set.VectorsAndVCVs(1:InitialVectorSet,:);
t1Set_0.Time = t1Set.Time(1:InitialVectorSet,:);

%----- keeping aside the optimal set - max DeltaTime -------- <- for initial vaules
dt = t1Set_0.Time - t0Set_0.Time;
i = find(dt==max(dt));
t0Set_00.crds = t0Set_0.crds;
t0Set_00.Points = t0Set_0.Points;
t0Set_00.vctrs = t0Set_0.vctrs(i,:);
t0Set_00.vcvs = t0Set_0.vcvs(i,:);
t0Set_00.VectorsAndVCVs = t0Set_0.VectorsAndVCVs(i,:);
t0Set_00.Time = t0Set_0.Time(i,:);

t1Set_00.crds = t1Set_0.crds;
t1Set_00.Points = t1Set_0.Points;
t1Set_00.vctrs = t1Set_0.vctrs(i,:);
t1Set_00.vcvs = t1Set_0.vcvs(i,:);
t1Set_00.VectorsAndVCVs = t1Set_0.VectorsAndVCVs(i,:);
t1Set_00.Time = t1Set_0.Time(i,:);

% enlarging data set
[t0Set_00] = AddPertrubedVectors_Random_MultiEpoch(t0Set_00,2,prtbMethod);
[t1Set_00] = AddPertrubedVectors_Random_MultiEpoch(t1Set_00,2,prtbMethod);

%% Thresholds and stopping conditions
% [V_th,Az_th,L_th] = computeTresholds_LockedFault(orig_t0Set,orig_t1Set);
% 
% V_th = V_th / nrmFctr;
% L_th = L_th / nrmFctr;
% % Az_th no need to reduce to nrmFctr
% pos_th = L_th;


V_th = 0.1/nrmFctr; %mm per year        % V_th = 0.1; %mm per year
L_th = 500/nrmFctr; %meters    % used to be 5000
Az_th = 1; %degrees
pos_th = 500/nrmFctr; %meters


% v_sc = V_th/2; %sc=stopping condition 
%% Feature Voting on selected fixed fault-locking positions
% posE = [64000,65000,66000]; posN = [49000,50000,51000]; % works great
% posE = [45000,55000,65000,75000,85000]; posN = [30000,40000,50000,60000,70000]; % works great!

% dpos = mean([(max(t0Set.crds(:,1)) - min(t0Set.crds(:,1))),(max(t0Set.crds(:,2)) - min(t0Set.crds(:,2)))]) / 10;
% n = min( [floor((max(t0Set.crds(:,1)) - min(t0Set.crds(:,1))) / dpos), floor((max(t0Set.crds(:,2)) - min(t0Set.crds(:,2))) / dpos)] );
% posN = linspace(min(t0Set.crds(:,1)) + dpos,max(t0Set.crds(:,1)) -dpos ,n);
% posE = linspace(min(t0Set.crds(:,2)) + dpos,max(t0Set.crds(:,2)) -dpos ,n);
% posE = floor(posE); posN = floor(posN);
%                                           older values
azStep  = 05;             % azStep = 0;      - 5
vStep   = 0.5/nrmFctr;    % vStep  = 0;      - 0.2
posStep = 2500/nrmFctr;   % step size [m]    - 2500
ldStep  = 2000/nrmFctr;   % ldStep = 0;      - 2000
LDmax   = 50000/nrmFctr;  % ------------------ max curstul depth



% % % % -------- setting up a structure for solutions ---------------
% % % solutions = table( NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,...
% % %     'VariableNames',{'posE' 'posN' 'az' 'v' 'ld' 'vtv' 'v_step' 'az_step' 'ld_step' 'v_th' 'l_th' 'az_th' } );
% % % allSolutios = solutions;
% % % allFilteredsSolutios = solutions; 



%% Initial Values:
%  ------- fault speed ------- 
% vs = sqrt((t1Set.vctrs(:,1)-t0Set.vctrs(:,1)).*(t1Set.vctrs(:,1)-t0Set.vctrs(:,1)) +...      %dn2
%           (t1Set.vctrs(:,2)-t0Set.vctrs(:,2)).*(t1Set.vctrs(:,2)-t0Set.vctrs(:,2))      );   %de2
      
vs = sqrt((t1Set.vctrs(1:InitialVectorSet,1)-t0Set.vctrs(1:InitialVectorSet,1)).*(t1Set.vctrs(1:InitialVectorSet,1)-t0Set.vctrs(1:InitialVectorSet,1)) +...      %dn2
          (t1Set.vctrs(1:InitialVectorSet,2)-t0Set.vctrs(1:InitialVectorSet,2)).*(t1Set.vctrs(1:InitialVectorSet,2)-t0Set.vctrs(1:InitialVectorSet,2))      );   %de2
vs = vs *1000 ./ (t1Set.Time(1:InitialVectorSet) - t0Set.Time(1:InitialVectorSet)); %mm per year        ./ DeltaTime - for Multi Epoch
% vs = ceil(vs*1e5)/1e5;

% v0 = 2*mean( vs ) * nrmFctr;    % more proper for true data with errors
% v0 = mean( vs ) * nrmFctr;    
v0 = max( vs ) * nrmFctr;    

v0 = v0 * 0.8; % to reduce some potential noise

% VFmax = 1.5 * ( v0 ) / nrmFctr;
% VFmin = 0.5 * ( v0 ) / nrmFctr;
VFmax = (v0 + 4*vStep*nrmFctr)/ nrmFctr;
VFmin = (v0 - 4*vStep*nrmFctr)/ nrmFctr;
% VFmax = min( 1.5 * ( v0 ) / nrmFctr, (v0 + 5*vStep)/ nrmFctr );
% VFmin = max( 0.5 * ( v0 ) / nrmFctr, (v0 - 5*vStep)/ nrmFctr );

VFs =  VFmin  : vStep      :  VFmax ;   %VFs = round(VFs,4);
disp(unique(vs.*nrmFctr)); disp(v0); disp(VFs.*nrmFctr);


% trying to compute fualt azimuth can help computing a single fault - but
% can cause a dis-converjing for multi fault
% az = Azimuth(t0Set.vctrs(:,2),t0Set.vctrs(:,1));    % Azimuth(de,dn)

% ldStep = ldStep * 3;    azStep = azStep * 3; vStep = vStep * 2;

LDs =     0   : ldStep     :  LDmax ;
AZs =   -90   : azStep     :  90    ;   AZs = pi./180.*AZs;   % used to be [-180,180]
% AZs = -180     : azStep     : 180   ;   AZs = pi./180.*AZs;   % used to be [-180,180]
VFs =  VFmin  : vStep      :  VFmax ;   %VFs = round(VFs,4);

%  ------- "fan-lock"-Position ------- 
% sInitPosTypeComp = 'covarage';
% sInitPosTypeComp = 'speed';
sInitPosTypeComp = 'SpecialCases';
switch sInitPosTypeComp
    case 'covarage'
        % by netwrok covarage
        mdE = round(mean(t0Set.Points(:,:).Grid_Easting) /posStep)*posStep;
        % mnE = floor(min( t0Set.Points(:,:).Grid_Easting) /posStep)*posStep;
        mxE = ceil( max( t0Set.Points(:,:).Grid_Easting) /posStep)*posStep;

        mdN = round(mean(t0Set.Points(:,:).Grid_Northing)/posStep)*posStep;
        % mnN = floor(min( t0Set.Points(:,:).Grid_Northing)/posStep)*posStep;
        mxN = ceil( max( t0Set.Points(:,:).Grid_Northing)/posStep)*posStep;

        tk = min(  [  (mxE-mdE)/posStep ,...
                     (mxN-mdN)/posStep] ) + 3;
        posE = mdE-tk*posStep : posStep : mdE+tk*posStep;    %posE = posE + 5000;
        posN = mdN-tk*posStep : posStep : mdN+tk*posStep; 

    case 'speed'
% % %         % by vectors max speed / deformation
% % %         vs0 = sqrt((t1Set.vctrs(:,1)-t0Set.vctrs(:,1)).*(t1Set.vctrs(:,1)-t0Set.vctrs(:,1)) +...      %dn2
% % %           (t1Set.vctrs(:,2)-t0Set.vctrs(:,2)).*(t1Set.vctrs(:,2)-t0Set.vctrs(:,2))      );   %de2
% % %         vs0 = vs0 *1000 ./ (t1Set.Time - t0Set.Time); %mm per year        ./ DeltaTime - for Multi Epoch
% % % 
% % %         [~,i] = min(vs0(1:InitialVectorSet).*nrmFctr);
% % %         s1 = t1Set.VectorsAndVCVs{i,1};
% % %         s2 = t1Set.VectorsAndVCVs{i,2};
% % %         
% % %         i1 = find(strcmpi(t1Set.Points{:,1},s1),1);
% % %         i2 = find(strcmpi(t1Set.Points{:,1},s2),1);
% % %         
% % %         ne0 = [mean(t1Set.crds([i1,i2],1)),   mean(t1Set.crds([i1,i2],2))];
% % %         
% % %         tk = 5;
% % %         posE = ne0(2)-tk*posStep : posStep : ne0(2)+tk*posStep;    %posE = posE + 5000;
% % %         posN = ne0(1)-tk*posStep : posStep : ne0(1)+tk*posStep; 
    case 'SpecialCases'
        % this is for center of earthquakes in Sea of Galile between 2010-2021
        
        % CSAR  IG05/12  E:189919.643, N:710505.891
        mdE = (254440 - 189919.643) / nrmFctr ;
        mdN = (751180 - 710505.891) / nrmFctr ;
        posStep = 1 / nrmFctr; %meter
            tk = 2;
        posE = mdE-tk*posStep : posStep : mdE+tk*posStep;    %posE = posE + 5000;
        posN = mdN-tk*posStep : posStep : mdN+tk*posStep; 
        
        az_0 = 0;
        azStep = 0.0002; %degrees
        AZs =   az_0-tk*azStep   : azStep     :  az_0+tk*azStep    ;   AZs = pi./180.*AZs;
end % switch


% posE = floor(min(t0Set.Points(:,:).Grid_Easting))  : posStep :  ceil(max(t0Set.Points(:,:).Grid_Easting)) ;
% posN = floor(min(t0Set.Points(:,:).Grid_Northing)) : posStep :  ceil(max(t0Set.Points(:,:).Grid_Northing));
% posE = mnE : posStep : mxE; 
% posN = mnN : posStep : mxE; 

fprintf('LDs\n');disp(LDs.*nrmFctr);
fprintf('Azs\n');disp(AZs.*180./pi); 
fprintf('VFs\n');disp(VFs.*nrmFctr); 
fprintf('posNs\n');disp(posN.*nrmFctr);
fprintf('posEs\n');disp(posE.*nrmFctr);
% disp(numel(LDs)); disp(numel(AZs)); disp(numel(VFs)); disp(numel(posN)); disp(numel(posE));
% disp(prod([numel(LDs); numel(AZs); numel(VFs); numel(posN); numel(posE)]));

% % % % % % % % % % % % % % % % % % % % VFs = sort([VFs, trueHz_Vel/nrmFctr]);
% % % % % % % % % % % % % % % % % % % % AZs = sort([AZs,trueAz*pi/180]);
% % % % % % % % % % % % % % % % % % % % LDs = sort([LDs, trueLd/nrmFctr]);
% % % % % % % % % % % % % % % % % % % % posE = sort([posE,truePosE/nrmFctr]);
% % % % % % % % % % % % % % % % % % % % posN = sort([posN,truePosN/nrmFctr]);
% fprintf('%f : %f :%f \n',VFmax/4,vStep,VFmax);
% fprintf('normalization Factor = %f\n',nrmFctr);
% disp(cumprod([numel(LDs),numel(AZs),numel(VFs),numel(posE),numel(posN)]));

% VFs = 7; vStep = 0;
% LDs = 10000; ldStep = 0;
% VFs = -10 : 1: 10;
% AZs = -10         : azStep     :10   ; AZs = pi./180.*AZs;
% VFs = 0    : vStep      : VFmax ;

        
        
%     % this is for sanity check
%     AZs = 0;
%     VFs = 6;
%     LDs = 10000;


%  curSol = table( NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,...
%           'VariableNames',{'posE' 'posN' 'az' 'v' 'ld' 'vtv' 'v_step' 'az_step' 'ld_step' 'v_th' 'l_th' 'az_th' } );
%      % steps and threshlods
dtStart_FulllCompute = datetime('now');dtStart_FulllCompute.Format = 'yyyy-MM-dd HH:mm:ss';
iteration = -1;

try %main try
while ( (posStep > pos_th) || (vStep > V_th) || (ldStep > L_th) || (azStep > Az_th)) % thersholds for values % ################ posStep 
%     % ------ check for stopping condition
%     if (vStep <= v_sc); break; end;
%     if (iteration > 7); break; end;
    iteration = iteration+1;
    fprintf('Iteration # %03i : \n',iteration);
    if (iteration>0); save([mPath2outputFolder,curRun,'\iter_',num2str(iteration),'.mat']); end;

    % ------ computing 
    dtNow = datetime('now');  dtNow.Format = 'yyyy-MM-dd HH:mm:ss';
    fprintf('%s - computing voting matrix....\n',string(dtNow));
    % H5 = NaN. *ones(numel(Azs),numel(VFs),numel(LDs),numel(posE),numel(posN)) ;
    switch (strFaultType)
        case 'Locked'
            if (iteration == 0)
                % for the initial guess sending the initial vectors only
                H5 = LockedFault_FeatureVoting5(t0Set_0,t1Set_0,posE(1:numel(posE)),posN(1:numel(posN)),AZs,VFs,LDs, V_th);
            else
                % for rest of computations - running with all original + pertrubed
                H5 = LockedFault_FeatureVoting5(t0Set,t1Set,posE(1:numel(posE)),posN(1:numel(posN)),AZs,VFs,LDs, V_th);
            end
        case 'Slip'
            H5 = SlipFault_FeatureVoting(t0Set,t1Set,posE(1:numel(posE)),posN(1:numel(posN)),AZs,VFs, V_th);
        otherwise
            error('hough_tests_5dim::wrong Fault Type');
    end %switch
    
    
    dtNow2 = datetime('now');  dtNow2.Format = 'yyyy-MM-dd HH:mm:ss';
    fprintf('time to create voting matrix [5-D] = %.4f seconds\n',seconds(dtNow2-dtNow));
    fprintf('numer of elements(H5) = %d\n',numel(H5));
 %{
         %strDrawSign = '';
%         switch (rem(iteration,5))
%             case 0
%                 strDrawSign = '+r';
%             case 1
%                 strDrawSign = '+g';
%             case 2
%                 strDrawSign = '+b';
%             case 3
%                 strDrawSign = '+c';
%             case 4
%                 strDrawSign = '+m';
%             otherwise
%                 strDrawSign = '.k';
%     %             iteration = 0;
%         end %switch (iteration)

    
%     if (isShowFigs == true)
%         colormap('jet');
%         figure(p)  , scatter3(aaa,vvv,lll,H3(fi)*numel(vs),H3(fi)./max(H3(fi)) );
%         figure(gcf), hold on;
% 
%         figure(gcf),xlabel('Azimuth [degrees]');
%         figure(gcf),ylabel('Fault Velocity [mm per year]');
%         figure(gcf),zlabel('Lock Depth [meters]');
%         figure(gcf),title(['Pos = [e=',num2str(posE(p)),',n=',num2str(posN(p)),'], maxCnocenzus = ',num2str( floor(max(H3(:))*100) ),'%']);
%     figure(gcf), view(-65,60);
%     end
%}            
            %% pass on solution possible solutions for consenzus testings (v'Pv)
    if strcmpi(strMode ,'SingleFault')
        %%%%%%%%% ------------------------------
        %             Single Fault Mode
        %%%%%%%%% ------------------------------
        % computing diffrences between multi-models and all-vectors in
        % mm-per-time [mm/year]
%% 
%{
        %%%%%%%%% ------------------------------
        %     Filter using VtV on csenarios
        %%%%%%%%% ------------------------------
        
        % initial filter by concenzus
        [fi] = SelectingConsenzus(H5,strMode ,'pareto');%'2d_95p'
        [am,vm,lm,pem,pnm] = ind2sub(size(H5), fi) ;   
        
        % initial filter by concenzus
        [fi] = SelectingConsenzus(H5,strMode ,'pareto');%'2d_95p'
        [am,vm,lm,pem,pnm] = ind2sub(size(H5), fi) ;   
            
        % deviations (v - vectors) for filtered cenarios
        [v_mat] = LockedFualt_Compute_VtVs(t0Set,t1Set,...
                    posE(pem).*nrmFctr,posN(pnm).*nrmFctr,AZs(am),VFs(vm).*nrmFctr,LDs(lm).*nrmFctr);
        v_mat = v_mat .* v_mat; % variances
        
        %             v_mat = v_mat'; % re-ordering for sum
        %             v_mat = sum(v_mat.*v_mat);
        %             i = find(v_mat == min(v_mat));
        min_v = min(v_mat(:));
        i = find(v_mat == min_v); % for optimal consenzus
        [i,~] = ind2sub(size(v_mat),i);
        af = am(i); lf = lm(i);vf = vm(i); pef = pem(i); pnf = pnm(i);
        af = reshape(af,[numel(af),1]); lf = reshape(lf,[numel(af),1]); vf = reshape(vf,[numel(af),1]);

%         fprintf('---min value From v_mat_squared-----\n')
%         
%         disp([  reshape(AZs(af),[numel(af),1])    .* 180./pi,...
%                 reshape(VFs(vf),[numel(vf),1])    .* nrmFctr,...
%                 reshape(LDs(lf),[numel(lf),1])    .* nrmFctr,...
%                 reshape(posE(pef),[numel(pef),1]) .* nrmFctr,...
%                 reshape(posN(pnf),[numel(pnf),1]) .* nrmFctr     ]  );
%}
%%      
%{
        %%% ------------------------------------------------------------------
        %%%  -----------   k-means for selecting best solution/s
        %%% ------------------------------------------------------------------
        ClusterNum = floor(numel(H5)/50000);  %1 scenario for each 10000 options.....
        rng('shuffle');
            % "fixing" azimuth and velocity to the form of "left lock system"
            aaa = AZs(am).*180./pi; vvv = VFs(vm); lll = LDs(lm);
            iii = aaa>90;
            aaa(iii) = aaa(iii)-180; vvv(iii) = vvv(iii)*-1;
            iii = aaa<-90;
            aaa(iii) = aaa(iii)+180; vvv(iii) = vvv(iii)*-1;
            aaa = aaa.*pi./180;
            eee = posE(pem); nnn = posN(pnm);
            % re-sizing the data for kmeas
            if (size(aaa,1) == 1); aaa = aaa';  end;
            if (size(vvv,1) == 1); vvv = vvv';  end;
            if (size(lll,1) == 1); lll = lll';  end;
            if (size(eee,1) == 1); eee = eee';  end;
            if (size(nnn,1) == 1); nnn = nnn';  end;
        rng('shuffle');
        [idx,C] = kmeans([aaa,vvv,lll,eee,nnn],ClusterNum,'Distance','cityblock','Replicates',20);
        %     %'sqeuclidean', 'cityblock', 'cosine','correlation', 'hamming'
        %     if (isShowFigs == true)
        %       figure(gcf), scatter3(C(:,1),C(:,2),C(:,3),'k+');
        %     end
        % for cenario - Centers [order : az, v, ld, posE, posN]
        [Cf] = filterCenters(C,idx, strMode,Az_th, V_th, L_th, pos_th);
        % Cf columns order: [ az, v, ld, posE, posN, opularity] 
        %       Cf is ordered by last column from large to small
        
        [vc_mat] = LockedFualt_Compute_VtVs(t0Set,t1Set,...
            Cf(:,4)'.*nrmFctr,Cf(:,5)'.*nrmFctr,Cf(:,1)',Cf(:,2)'.*nrmFctr,Cf(:,3)'.*nrmFctr);
        vc_mat = vc_mat .* vc_mat;
        min_vc = min(vc_mat(:));
        ic = find(vc_mat == min_vc); % for optimal consenzus
        [ic,~] = ind2sub(size(vc_mat),ic);
        
        if (numel(ic) == 1)
            ac = Cf(ic,1); vc = Cf(ic,2); lc = Cf(ic,3); ec = Cf(ic,4); nc = Cf(ic,5);
        else
            ic = 1; %largest in popularity
            ac = Cf(ic,1); vc = Cf(ic,2); lc = Cf(ic,3); ec = Cf(ic,4); nc = Cf(ic,5); 
        end
        
        fprintf('---min value From KMeans Clustering---\n');
        fprintf(' Az, V, depth, posE, posN\n');
        disp( round(   [Cf(ic,1)*180/pi ,Cf(ic,2:end-1)*nrmFctr,Cf(ic,end)] ) );
        fprintf('---all values From KMeans Clustering---\n');
        disp( round(   [Cf(:,1).*180./pi,Cf(:,2:end-1)*nrmFctr ,Cf(:,end)] )  );
%}        
 %%     
        %%% ------------------------------------------------------------------
        %%%  -----------   best csenarion selection using peaks method
        %%% ------------------------------------------------------------------
        
        % get peaks for the voting matrix
        maxPeaks = ceil(numel(H5)*0.0001);
        Cp_idx = GetVotingPeaks(H5, maxPeaks);
        [acp,vcp,lcp,pecp,pncp] = ind2sub(size(H5), Cp_idx);
        fprintf('--- Peaks From Voting Matrix ---\n');
        switch (strFaultType)
        case 'Locked'
            fprintf('--- Az[deg], V_F[mm/year], LD[m], posE[m], posN[m]---\n');
            if numel(acp) > 100
                fprintf('[displaying only first 100 solutions]\n');
                disp([AZs(acp(1:100))'.*180./pi,  VFs(vcp(1:100))'.*nrmFctr,  LDs(lcp(1:100))'.*nrmFctr,  posE(pecp(1:100))'.*nrmFctr,  posN(pncp(1:100))'.*nrmFctr ]) ;
            else
                disp([AZs(acp)'.*180./pi,  VFs(vcp)'.*nrmFctr,  LDs(lcp)'.*nrmFctr,  posE(pecp)'.*nrmFctr,  posN(pncp)'.*nrmFctr ]) ;
            end
        case 'Slip'
            fprintf('--- Az[deg], V_F[mm/year], posE[m], posN[m]---\n');
            if numel(acp) > 100
                fprintf('[displaying only first 100 solutions]\n');
                disp([AZs(acp(1:100))'.*180./pi,  VFs(vcp(1:100))'.*nrmFctr,  posE(pecp(1:100))'.*nrmFctr,  posN(pncp(1:100))'.*nrmFctr ]) ;
            else
                disp([AZs(acp)'.*180./pi,  VFs(vcp)'.*nrmFctr,  posE(pecp)'.*nrmFctr,  posN(pncp)'.*nrmFctr ]) ;
            end
        otherwise
            error('hough_tests_5dim::wrong Fault Type');
        end %switch
        
        
        % re-selecting by vtv
        switch (strFaultType)
            case 'Locked'
                [v_mat] = LockedFualt_Compute_VtVs(t0Set,t1Set,...
                                posE(pecp).*nrmFctr,posN(pncp).*nrmFctr,AZs(acp),VFs(vcp).*nrmFctr,LDs(lcp).*nrmFctr);
                v_mat = v_mat .* v_mat; % variances

                %             v_mat = v_mat'; % re-ordering for sum
                %             v_mat = sum(v_mat.*v_mat);
                %             i = find(v_mat == min(v_mat));
                min_v = min(v_mat(:));
                i = find(v_mat == min_v); % for optimal consenzus
                [i,~] = ind2sub(size(v_mat),i);
                af1cp = acp(i); lf1cp = lcp(i);vf1cp = vcp(i); pef1cp = pecp(i); pnf1cp = pncp(i);
                af1cp = reshape(af1cp,[numel(af1cp),1]); lf1cp = reshape(lf1cp,[numel(af1cp),1]); vf1cp = reshape(vf1cp,[numel(vf1cp),1]);
                pef1cp = reshape(pef1cp,[numel(pef1cp),1]); pnf1cp = reshape(pnf1cp ,[numel(pnf1cp ),1]);
                fprintf('--- min_VtV Peak FromPeaks of Voting Matrix ---\n');
                fprintf('--- Az[deg], V_F[mm/year], LD[m], posE[m], posN[m]---\n');
                disp([AZs(af1cp)'.*180./pi,  VFs(vf1cp)'.*nrmFctr,  LDs(lf1cp)'.*nrmFctr,  posE(pef1cp)'.*nrmFctr,  posN(pnf1cp)'.*nrmFctr ]) ;
                
                try
                    [vcv] = ComputeVCVfromVotingMatrix5(H5,2,3,af1cp,pef1cp,pnf1cp);
                    vcv.VelHzAcc = vcv.dX .* vStep .* nrmFctr;
                    vcv.LdAcc = vcv.dY .* ldStep .* nrmFctr;
                    disp(vcv(:,end-1:end));
                catch ME
                    fprintf('couldnt compute vcv\n');
                    fprintf(ME.message);
                    fprintf('\ncomputation continoues\n');
                end
                % 95% acceptance
                figure(gcf),yticklabels(string(VFs.* nrmFctr)); figure(gcf),ylabel('Velocities'); figure(gcf),ytickangle(45);
                figure(gcf),xticklabels(string(LDs.* nrmFctr)); figure(gcf),xlabel('Lock Depths');figure(gcf),xtickangle(45);
                figure(gcf),title(['Iteration #',num2str(iteration),' - relative voting over 95%']);
                h = figure(gcf);
                % voting matrix
                figure(h.Number-1),yticklabels(string(VFs.* nrmFctr));figure(h.Number-1),ylabel('Velocities');  figure(h.Number-1),ytickangle(45);
                figure(h.Number-1),xticklabels(string(LDs.* nrmFctr));figure(h.Number-1),xlabel('Lock Depths'); figure(h.Number-1),xtickangle(45);
                figure(h.Number-1),title(['Iteration #',num2str(iteration),' - voting ratios']);
                hFig{h.Number-1} = figure(h.Number-1);
                hFig{h.Number  } = figure(h.Number  );
                
                
                
            case 'Slip'
                [v_mat] = SlipFualt_Compute_VtVs(t0Set,t1Set,...
                                posE(pecp).*nrmFctr,posN(pncp).*nrmFctr,AZs(acp),VFs(vcp).*nrmFctr);
                v_mat = v_mat .* v_mat; % variances

                %             v_mat = v_mat'; % re-ordering for sum
                %             v_mat = sum(v_mat.*v_mat);
                %             i = find(v_mat == min(v_mat));
                min_v = min(v_mat(:));
                i = find(v_mat == min_v); % for optimal consenzus
                [i,~] = ind2sub(size(v_mat),i);
                af1cp = acp(i); vf1cp = vcp(i); pef1cp = pecp(i); pnf1cp = pncp(i);
                af1cp = reshape(af1cp,[numel(af1cp),1]);  vf1cp = reshape(vf1cp,[numel(vf1cp),1]);
                pef1cp = reshape(pef1cp,[numel(pef1cp),1]); pnf1cp = reshape(pnf1cp ,[numel(pnf1cp ),1]);
                fprintf('--- min_VtV Peak FromPeaks of Voting Matrix ---\n');
                fprintf('--- Az[deg], V_F[mm/year], posE[m], posN[m]---\n');
                disp([AZs(af1cp)'.*180./pi,  VFs(vf1cp)'.*nrmFctr,  posE(pef1cp)'.*nrmFctr,  posN(pnf1cp)'.*nrmFctr ]) ;
            otherwise
                error('hough_tests_5dim::wrong Fault Type');
        end %switch
        
 %%
        % --------------------------------------
        %     Prepearing for next iteration
        % --------------------------------------
        % resetting values for next iteration
        switch (strFaultType)
            case 'Locked'
                azStep = azStep / 2;
                ldStep = ldStep / 2;      ldStep = ldStep * nrmFctr;
                vStep = vStep / 2;        vStep = vStep * nrmFctr;
                posStep = posStep / 2;    posStep = posStep * nrmFctr;
            case 'Slip'
                azStep = azStep / 2;
                vStep = vStep / 2;        vStep = vStep * nrmFctr;
                posStep = posStep / 2;    posStep = posStep * nrmFctr;
        otherwise
            error('hough_tests_5dim::wrong Fault Type');
        end %switch
        %------ bring back from normalization (using multiplier, function or look-up-table)
%         % best solution by KMEANS
%         vc = vc.*nrmFctr; lc=lc.*nrmFctr; ec=ec.*nrmFctr; nc=nc.*nrmFctr;
        
        %--------------0------------------0----------------------------- 
        % best solution for next iteration relys on voting peaks method
        %--------------0------------------0----------------------------- 
% %         % max voting values
%         ac = AZs(acp(1));  vc = VFs(vcp(1)).*nrmFctr;  lc=LDs(lcp(1)).*nrmFctr;  ec=posE(pecp(1)).*nrmFctr;  nc=posN(pncp(1)).*nrmFctr;
        % mean peaks voting values
%         ac = mean(AZs(acp(:)));  vc = mean(VFs(vcp(:))).*nrmFctr;  lc=mean(LDs(lcp(:))).*nrmFctr;  ec=mean(posE(pecp(:))).*nrmFctr;  nc=mean(posN(pncp(:))).*nrmFctr;
        % minVtV from peaks voting values
        switch (strFaultType)
           case 'Locked'
               ac = AZs(af1cp)'.*180./pi;  vc = VFs(vf1cp)'.*nrmFctr;  lc = LDs(lf1cp)'.*nrmFctr;  ec = posE(pef1cp)'.*nrmFctr;  nc = posN(pnf1cp)'.*nrmFctr;
               fprintf('selected best solution for next iteration:\n');
               fprintf('--- Az[deg], V_F[mm/year], LD[m], posE[m], posN[m]---\n');
               if (numel(ac) > 100)
                    fprintf('[displaying only first 100 solutions]\n');
                    disp([ac(1:100),vc(1:100),lc(1:100),ec(1:100),nc(1:100)]);
               else
                    disp([ac,vc,lc,ec,nc]);
               end
               
               if (numel(ac) > 1)
                    fprintf('avaraging the selected best-solutions');
                    ac = mean(ac);
                    vc = mean(vc);
                    lc = mean(lc);
                    ec = mean(ec);
                    nc = mean(nc);
                   
                    fprintf('\n\nSelected Avg solution for next iteration:\n');
                    fprintf('--- Az[deg], V_F[mm/year], LD[m], posE[m], posN[m]---\n');
                    disp([ac,vc,lc,ec,nc]);
               end
           case 'Slip'
               ac = AZs(af1cp)'.*180./pi;  vc = VFs(vf1cp)'.*nrmFctr;  ec = posE(pef1cp)'.*nrmFctr;  nc = posN(pnf1cp)'.*nrmFctr;
               fprintf('selected best solution for next iteration:\n');
               fprintf('--- Az[deg], V_F[mm/year], posE[m], posN[m]---\n');
               if (numel(ac) > 100)
                    fprintf('[displaying only first 100 solutions]\n');
                    disp([ac(1:100),vc(1:100),ec(1:100),nc(1:100)]);
               else
                    disp([ac,vc,ec,nc]);
               end
           otherwise
               error('hough_tests_5dim::wrong Fault Type');
        end %switch
        
        fprintf('--------------------------------------\n');
        fprintf('---------- end iteration -------------       iteration # %d \n',iteration);
        fprintf('--------------------------------------\n');
        n = 5;
        switch (strFaultType)
            case 'Locked'
                    AZs = 180./pi.*AZs;                                                           %#ok
                AZs = ac  -n*azStep  : azStep  : ac + n*azStep;       AZs = pi./180.*AZs;
                VFs = vc  -n*vStep   : vStep   : vc + n*vStep ;       VFs = round(VFs,SigDig);
                LDs = lc  -n*ldStep  : ldStep  : lc + n*ldStep;       LDs = LDs(LDs>0);           %#ok
                posE = ec -n*posStep : posStep : ec + n*posStep;      %posE = posE + 5000;
                posN = nc -n*posStep : posStep : nc + n*posStep; 
                
%                 AZs = min(ac)  -n*azStep  : azStep  : max(ac) + n*azStep;       AZs = pi./180.*AZs;
%                 VFs = min(vc)  -n*vStep   : vStep   : max(vc) + n*vStep ;       VFs = round(VFs,SigDig);
%                 LDs = min(lc)  -n*ldStep  : ldStep  : max(lc) + n*ldStep;       LDs = LDs(LDs>0);           %#ok
%                 posE = min(ec) -n*posStep : posStep : max(ec) + n*posStep;      %posE = posE + 5000;
%                 posN = min(nc) -n*posStep : posStep : max(nc) + n*posStep; 
                
                fprintf('new input :\n');
                fprintf('Velocities:  ') ;  disp(VFs);
                fprintf('Azimuts:  ') ;  disp(180./pi.*AZs);
                fprintf('Lock-Depths:  ') ;  disp(LDs);
                fprintf('East-Pos:  ') ;  disp(posE);
                fprintf('North-Pos:  ') ;  disp(posN);
            case 'Slip'
                     AZs = 180./pi.*AZs;                                                          %#ok
                AZs = ac  -n*azStep  : azStep  : ac + n*azStep;       AZs = pi./180.*AZs;
                VFs = vc  -n*vStep   : vStep   : vc + n*vStep ;       VFs = round(VFs,SigDig);
                posE = ec -n*posStep : posStep : ec +n*posStep;    %posE = posE + 5000;
                posN = nc -n*posStep : posStep : nc +n*posStep; 
                
                fprintf('new input :\n');
                fprintf('Velocities:  ') ;  disp(VFs);
                fprintf('Azimuts:  ') ;  disp(180./pi.*AZs);
                fprintf('East-Pos:  ') ;  disp(posE);
                fprintf('North-Pos:  ') ;  disp(posN);
            otherwise
                error('hough_tests_5dim::wrong Fault Type');
        end %switch
        
%         if (iteration > 3) %-- --------------- add true values for checkup
%            AZs = sort([AZs, trueAz*pi/180]) ;
%            VFs = sort([VFs, trueHz_Vel]) ;
%            LDs = sort([LDs, trueLd]) ;
%            posE = sort([posE, truePosE]) ;
%            posN = sort([posN, truePosN]) ;
%         end
        
        % set normalization (using multiplier, function or look-up-table
        switch (strFaultType)
            case 'Locked'
                %Az is in radians anyway
                VFs  = VFs ./nrmFctr;   ldStep = ldStep / nrmFctr;
                LDs  = LDs ./nrmFctr;   vStep = vStep / nrmFctr;
                posE = posE./nrmFctr;   posStep = posStep / nrmFctr;
                posN = posN./nrmFctr;
            case 'Slip'
                %Az is in radians anyway
                VFs  = VFs ./nrmFctr;   vStep = vStep / nrmFctr;
                posE = posE./nrmFctr;   posStep = posStep / nrmFctr;
                posN = posN./nrmFctr;
            otherwise
                error('hough_tests_5dim::wrong Fault Type');
        end %switch
        
        dumb = -1;  %#ok
        continue;
        % continue <<<<<<-------------------------<<<<<<<<<<<<<--------------     continue
        
    else %strcmpi(strMode ,'MultiFault')
        %%%%%%%%% ------------------------------
        %             Multi Fault Mode
        %%%%%%%%% ------------------------------
        error('Houhg_tests_ext5 :: selecting & detecting multiple faults un-handeled yet');
        %%
    end %if strcmpi(strMode ,'___STRING__')
    
end % while -> thersholds for values are too large

%% re-computing position of locking of the fault


%%   
        fprintf('Az diff between known true value center and estimated value:\n');
        % Azimuth(de,dn)
        az_checkup = Azimuth(truePosE - posE(pef1cp)'.*nrmFctr ,truePosN - posN(pnf1cp)'.*nrmFctr);
        disp( (az_checkup - trueAz.*pi./180).*180./pi );
        
        fprintf('********\n');
        fprintf('********\n');
        fprintf('\n');
        fprintf('final selected values are:\n');
        fprintf('Az = %.3f, V = %.3f, Ld = %.3f, posE = %.3f, posN = %.3f\n',ac(1),vc(1),lc(1),ec(1),nc(1));
        
        switch (strFaultType)
            case 'Locked' 
            fprintf('checkup - true values:\n');
            fprintf('Az = %.3f, V = %.3f, Ld = %.3f, posE = %.3f, posN = %.3f\n',trueAz, trueHz_Vel, trueLd, truePosE, truePosN);
        % exiting the current poisitio lock for a fault model
        case 'Slip'
            fprintf('checkup - true values:\n');
            fprintf('Az = %.3f, V = %.3f, posE = %.3f, posN = %.3f\n',trueAz, trueHz_Vel, truePosE, truePosN);
        end %switch
        dumb = -1; %#ok
        dtEnd_FulllCompute = datetime('now'); dtEnd_FulllCompute.Format = 'yyyy-MM-dd hh:mm:ss';
        fprintf('time to full process = %.4f seconds ( %.4f minutes)\n',seconds(dtEnd_FulllCompute-dtStart_FulllCompute),minutes(dtEnd_FulllCompute-dtStart_FulllCompute))
        
        
        diary off; %%%%%%%%%%%%%%%%%% end save of diary 
        fprintf('stopping saving diary...');
        save([mPath2outputFolder,curRun,'\iter_',num2str(iteration+1),'.mat']);
        error('======================================= stopping here ....... not an error, just a break =======================================');
    
   

fprintf('\n\n....... voting process ended! .......\n\n');       return;
catch ME_main
    fprintf(ME_main.identifier);
    fprintf('\n');
    fprintf(ME_main.message);
    fprintf('\n');
    fprintf(ME_main.stack);
    fprintf('\n');
    fprintf(ME_main.cause);
end % main try

diary('off');