% plot RF wiggles before or after filter

clear all;
FADERDir = '/home/lic0a/src/FADER-0.1.0/';
addpath([FADERDir 'Functions/']);

% plot filtered RF
% data=load('filteredRF/MP0203.mat');
 data=load('filteredRF/MP0502.mat');
% data=load('filteredRF/BTHS.mat');

data_before=data.R;
data_filter=data.R_flted;
t=data.t;

% figure;
% x0=500;
% y0=300;
% width=800;
% height=900;
% set(gcf,'position',[x0,y0,width,height])
% 
% L=length(data_before(:,1));
% for i=1:L;
%     subplot(L,1,i)
%     plot(t,data_before(i,:)/max(data_before(i,:)));
% %     set(gca,'XTick',[])
%     hold on;
%     plot(t,data_filter(i,:)/max(data_filter(i,:)));
%     xlim([-10,50]);
% end

network = 'KG';
%         station = 'MP026';
%          station = 'MP0203';
         station = 'MP0502';
% MP02: a=2.5 >85% ; MP021 >85%, a=1.0 10 events;
% MP022: a=1.0 MP02, stack >80%, 16 events before stack, 9 after stack;  
% MP023: a=2.5 MP02 stack > 80%;MP024: a=2.5 MP023+MP04; 
% MP025: a=2.5 MP02+MP04; MP026: a=5.0 >85%;
% MP02l1: a=2.5, lp 0.05-1 Hz MP02;
% MP0201: lp 2Hz, a=1 > 80%; MP0202: lp 2Hz, a=2.5 >80%;
% MP0203: lp 2Hz, a=1 >80%, all M>5.8 events; MP0204: a=2.5, same as MP0203
% MP0501: MP05 2022, 2023, a=1.0; MP05 2022-2023, a=2.5;
%         station = 'BOQS';
%         station = 'BTHS';

    epiDistRange = [30 95];
    resolutionFactor = 1;
    RFOUTDIR  = [FADERDir 'Data/MTCRF/'];
%     [rrfAmpArray, timeAxisHD, binAxisHD] = ...
%         loadAndPrepRF(network, station, resolutionFactor, epiDistRange, RFOUTDIR);
    [rrfAmpArray, timeAxisHD, binAxisHD] = ...
        loadAndPrepRF2(network, station, resolutionFactor, epiDistRange, RFOUTDIR);
    R = rrfAmpArray; t = timeAxisHD; y = binAxisHD; nY = length(y);
    [epiDist, rayP] = raypToEpiDist(y, 1, FADERDir);
    

   
%% Plot wiggle, modified from RFWigglePlot_any.m
% plot RF waveforms Y asix as number

% if plot filtered
% R=data_filter; 
R=data_before;

% Prepping values to be plotted from RF


tWin = [-5 30];
epiDistRange = [30 95];

% RF variables
nY = length(y);

tStart = tWin(1); tEnd = tWin(2);

%Time window to smooth over where location conversion and reverberation
it = find(t>tStart, 1); %Removing negative time
endt = find(t>tEnd, 1);

% Summary stack
for iY = 1:nY
    R(iY,:) = detrend(R(iY,:));
end
RR = R(:,it:endt);
sumR = sum(RR,  1);

% Plot pure RF traces with out weighting

hh = figure(10);
clf;

set(gcf,'position',[50,50,800,1200]);

tshft = 0;
mm = max(abs(R(:,it:endt)), [] ,'all');

% plot prediected arrival times
Vp=6.5;
% H=27;      % crust from SRF study (Hasan 2007) 
% ratio=1.64;  % ratio from Tang 2016
H=22;
% ratio=1.6;  % a=1
ratio=1.8;  % a=2.5
% ratio=1.63;
% H=30;
Vs=Vp/ratio;

% [tPs, tPps, tPss] = travelTimes(Vp, Vs, H, rayP, 1);
% Estimate PS conversion time from LAB, PlS
[tPs, tPps, tPss, tPls] = travelTimes3(Vp, Vs, H, rayP, 1);


% maxY = max(rayP);
% minY = min(rayP);
% maxY = 90;  % for Y as distance
% minY = 30;
maxY = nY; % for Y as number
minY = 0;

for iY = 1:nY
    
    mm = max(abs(R(iY,it:endt)));
    Rn = R(iY,it:endt) ./ mm;
    Rn = Rn - mean(Rn);
    Tn = t(it:endt); sizeT = length(Tn);
%     Tn = Tn'; % for RFLoad2;
    
    yLev = (nY-iY);  % for Y as number
%     yLev = (rayP(iY)-minY)/(maxY-minY)*nY;
%     yLev = (epiDist(iY)-minY)/(maxY-minY)*nY-1;  % for Y as distance
    yVec = repmat(yLev, 1, sizeT);
    
    
    plot([tPs(iY),tPs(iY)], [yLev, 1+yLev],'k'); hold on;
    plot([tPls(iY),tPls(iY)], [yLev, 1+yLev],'r'); hold on;
    plot([tPps(iY),tPps(iY)], [yLev, 1+yLev],'k'); hold on;
    plot([tPss(iY),tPss(iY)], [yLev, 1+yLev],'k'); hold on
    
    jbfill(Tn, max(Rn+yLev, yLev), yVec, [0 0 1],'k', 1, 1.0);
    jbfill(Tn, min(Rn+yLev, yLev), yVec, [1 0 0],'k', 1, 1.0);
%     jbfill(Tn, min(Rn+yLev, yLev), yVec, [1 1 1],'k', 1, 1.0);

    
end


% Plotting axis
xlim([0 tEnd])
% yticks([0:5:nY]);
yticks(linspace(0, nY, (nY/5)+1))
% set(gca,'yticklabel', floor(linspace(epiDistRange(1),epiDistRange(2),(nY/5)+1)))
set(gca,'yticklabel', floor(linspace(minY,maxY,(nY/5)+1)))
% set(gca,'xticklabel', '')
ylim([-1 nY+1]);  % for Y as number
xlabel('Time (s)','FontSize', 24)
% ylabel('Epicentral distance (deg)','FontSize', 20)

grid on




