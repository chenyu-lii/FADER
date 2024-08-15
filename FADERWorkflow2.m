%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 0. Parameters setup %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;

% FADERDir = '/scratch/tolugboj_lab/Prj4_Nomelt/FADER/';
FADERDir = '/home/lic0a/src/FADER-0.1.0/';
addpath([FADERDir 'Functions/']);

% kd_thresh = 2;
kd_thresh = 1;
man_determine = 1; % mannually check autocorrelation plot
tolerance = 0; % Echo delay time from autocorrelation and homomorphic

saveRF = 1;
saveDir = strcat(FADERDir,'filteredRF/');
% savename = 'NE68';
savename = 'MP0203';
% savename = 'MP02';
% savename = 'BOQS';
% savename = 'BTHS';
% savename = 'LT10';
% smoothopt = 1; % smooth filtered RF
smoothopt = 1;

Hkopt = 1;

% Vp = 6.4; % of crust
% Hwin = [25 55];
Vp = 6.5; % of crust
Hwin = [20 50];
% Hwin = [15 45];
% Hwin = [0 30];
kwin = [1.6 1.9];
resfac = 100;
% pWAc = [0.6 0.3 -0.1]; % weighting factors for H-k stacking
pWAc = [0.6 0.2 -0.2];
% pWAc = [0.9 0.1 0];
wid = 0.1; % width of Gaussian windows in H-k stacking

man_PbS = 1; % mannually confirm PbS arrival for H-k stacking

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 1. Load in RF data %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modify this part to load in your own data

dataopt = 2;

switch dataopt
    
    case 1 % Telewavesim Synthetics
        
        RFmat = load(strcat(FADERDir,'Data/Synthetics/M1_1.0Hz_syn.mat'));
        R = RFmat.rRF;
        t = RFmat.time; y = RFmat.garc; nY = length(y);
        R1 = R(2,:);
        [epiDist, rayP] = raypToEpiDist(y, 1, 1, localBaseDir);
        
    case 2 % Real Data from MTC RF
        
%         network = 'XO';
%         station = 'LT10';
%         network = 'YP';
%         station = 'NE68';
        network = 'KG';
        station = 'MP0203';
%          station = 'MP02';
% MP02: a=2.5 >85% ; MP021 a=1.0 10 events; 
% MP022: a=1.0 MP02, stack >80%, 16 events before stack, 9 after stack; 
% MP023: a=2.5 MP02 stack > 80%;MP024: a=2.5 MP023+MP04; 
% MP025: a=2.5 MP02+MP04; MP026: a=5.0 >85%;
% MP02l1: a=2.5, lp 0.05-1 Hz MP02;
% MP0201: lp 2Hz, a=1 > 80%; MP0202: lp 2Hz, a=2.5 >80%;
% MP0203: lp 2Hz, a=1 >80%, all M>5.8 events; MP0204: a=2.5, same as MP0203
%         station = 'BOQS';
%         station = 'BTHS';
        epiDistRange = [30 95];
        resolutionFactor = 1;
        RFOUTDIR  = [FADERDir 'Data/MTCRF/'];
%         [rrfAmpArray, timeAxisHD, binAxisHD] = ...
%             loadAndPrepRF(network, station, resolutionFactor, epiDistRange, RFOUTDIR);
        [rrfAmpArray, timeAxisHD, binAxisHD] = ...
            loadAndPrepRF2(network, station, resolutionFactor, epiDistRange, RFOUTDIR);
        R = rrfAmpArray; t = timeAxisHD; y = binAxisHD; nY = length(y);
        R1 = R(1,:);
%         R1 = R(7,:);
        [epiDist, rayP] = raypToEpiDist(y, 1, FADERDir);
        
end


%%
% modify the filter process, for each trace do the filter from it's own
% autocorrelation
% Chenyu Li

% No need to edit the rest of this code %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 2. Run Determination test             %
% Step 3. Get tlag & r0 from autocorrelation %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

L=length(y);
N=length(R);
R_flted=zeros(L,N);

% for j=1:2;
for j=1:L;
    
    R1=R(j,:);

    [SedPre,tlaga,r0,tac,ac,dsin] = ...
        DeterminationTest_func(R1,t,kd_thresh,man_determine);

    if SedPre == 0
        continue
        if Hkopt == 0
            fprintf('No reverberation detected. Exiting FADER.\n');
        else
            tlag = 0;
            tPbS = 0;
            fprintf('No reverberation detected. Will run H-k stacking directly.\n');
        end
    else
        fprintf('Reverberation detected. Analyzing ...\n');

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Step 4. Run cepstrum analysis to get tlag and compare %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        tlag=tlaga;

%         if tlaga < 1
%             twin_cstack = [0 2];
%         else
%             twin_cstack = [tlaga-1 tlaga+1];
%         end
% 
%         [tlagc,cstackt,cstackA] = ceps_func(R1,t,twin_cstack);
% 
%         if abs(tlagc - tlaga) < tolerance * max(tlagc, tlaga)
%             tlag = 0.5 * (tlagc + tlaga);
%         else
%             if man_determine   
%                 tlag = tlag_check(tac,ac,dsin,tlaga,cstackt,cstackA,tlagc,tolerance);
%             else
%                 tlag = tlaga;
%             end
%         end

        fprintf('tlag = %3.2f s, r0 = %3.2f.\n',tlag,r0);
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Step 5. Filter the RF using determined parameters %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        fprintf('\nFiltering data ...\n');
        R_flted(j,:) = filterRF_FADER(R1,t,tlag,r0,smoothopt);

    end
    
end

    if saveRF
        fprintf('Saving filtered RF ...\n');
        nname = strcat(saveDir,savename,'.mat');
        save(nname,'R','R_flted','t','y');
    end
    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 6. Optional H-k Stacking %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Hkopt
    
    if SedPre == 1
%         tPbS = tPbS_Confirm(R, t, y, man_PbS);
        tPbS = tPbS_Confirm(R, t, y, 0);
        R2stack = R_flted;
%         R2stack = R_flted(4,:);
    else
        R2stack = R;
    end
    
% H-k for crust
    [stackArray, Hrange, krange, HBest, kBest] = ...
        HkStacking(R2stack, t, rayP, Vp, Hwin, kwin, resfac, pWAc, wid, ...
        SedPre, tlag, tPbS,1);
%     


% check tPbs for each trace individually

%  [stackArray, Hrange, krange, HBest, kBest] = ...
%         HkStacking3(R2stack, t, rayP, Vp, Hwin, kwin, resfac, pWAc, wid, ...
%         SedPre, tlag, tPbS,1);

        
%% H-k for sediment
%     Vpd=4.5;
%     Hwin_d=[0 10];
%     kwin_d=[1.5 5];
%     [stackArray, Hrange, krange, HBest, kBest] = ...
%         HkStacking2(R2stack, t, rayP, Vpd, Hwin_d, kwin_d, resfac, pWAc, wid, ...
%         SedPre, tlag, tPbS,1);
end

%% H-k for no de-reverberation
   R2stack = R;
   SedPre = 0;
   tlag = 0;
   tPbS =0;
[stackArray, Hrange, krange, HBest, kBest] = ...
        HkStacking(R2stack, t, rayP, Vp, Hwin, kwin, resfac, pWAc, wid, ...
        SedPre, tlag, tPbS,1);

fprintf('\nExiting FADER.\n');
