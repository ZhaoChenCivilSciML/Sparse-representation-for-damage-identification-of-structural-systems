clear
clc
close all

% tic
%% Define parameters
NoiseLevel = 10;
SensorInd = [3 4 5 6 9 10 15 16 17 18 23 24 25 26]; 
Fs = 1400;
N = 60*Fs;
NoEle = 31;
NoDOF = 28;

% Define stiffness reduction coefficients
rng('default')

ThetaIntact = (rand(1,NoEle)-0.5)*0.2;

ThetaDamaged = zeros(1,NoEle);
ThetaDamaged(1) = 0.2;
ThetaDamaged(15) = 0.15;
ThetaDamaged(27) = 0.15;

ThetaDamagedCombo = ThetaIntact+ThetaDamaged-ThetaIntact.*ThetaDamaged;

%% Compute intact and damaged data
ModeInd1 = 1:3; 
[agIntact,accIntact,~,~,SIntact] = StructuralResponses_eig(N,NoiseLevel,SensorInd,...
    ThetaIntact,Fs,ModeInd1); 
ModeInd2 = 1:3; 
[agDamaged,accDamaged,~,~,SDamaged] = StructuralResponses_eig(N,NoiseLevel,SensorInd,...
    ThetaDamagedCombo,Fs,ModeInd2); 

%% Signal Filter
% Apply this if noisy
% Transfer Function design
% [bb, aa] = butter(5, [0.1/(Fs/2) 600/(Fs/2)], 'bandpass');
% Type 1: filter
% agIntact1 = filter(bb, aa, agIntact);
% accIntact1 = filter(bb, aa, accIntact);
% agDamaged1 = filter(bb, aa, agDamaged);
% accDamaged1 = filter(bb, aa, accDamaged);
% Type 2: filtfilt
% agIntact1 = filtfilt(bb, aa, agIntact');
% accIntact1 = filtfilt(bb, aa, accIntact');
% agDamaged1 = filtfilt(bb, aa, agDamaged');
% accDamaged1 = filtfilt(bb, aa, accDamaged');

% Zero-Pole-Gain design
% [z,p,k] = butter(5, [0.1/(Fs/2) 600/(Fs/2)],'bandpass');
% [sos,g] = zp2sos(z,p,k);
% agIntact1 = sosfilt(sos,agIntact)*g;
% accIntact1 = sosfilt(sos,accIntact)*g;
% agDamaged1 = sosfilt(sos,agDamaged)*g;
% accDamaged1 = sosfilt(sos,accDamaged)*g;

% Dont apply noise filtering but the following.
agIntact1 = agIntact;
accIntact1 = accIntact;
agDamaged1 = agDamaged;
accDamaged1 = accDamaged;


%% Eigen variables identification
% Type 1: OKIDex
% [FreqIntact,ShapeIntact] = OKIDex(agIntact1,accIntact1,Fs);
% [FreqDamaged,ShapeDamaged] = OKIDex(agDamaged1,accDamaged1,Fs);

% Type 2: Directly load
load('intact.mat')
load('damage.mat')

%% Model updating
% ReguFlag = 1: Sparse Bayesian learning (L1); ReguFlag = 2:
% Bayesian learning (L2); ReguFlag = 3: L curve.
disp('Model updating begins')
ReguFlag = 2; 
[ThetaHistoryIntact,~,ThetaIntactCovHistory_approx] = Regu_eig(N,ReguFlag,SensorInd,Fs,zeros(1,NoEle),...
    FreqIntact(ModeInd1),ShapeIntact(:,ModeInd1),ModeInd1,NoDOF);

% Type 2: Compute marginal posterior of individual SRC
R = mvnrnd(-ThetaHistoryIntact(end,:),ThetaIntactCovHistory_approx(:,:,end),1000);
pd = {};
for i = 1:size(-ThetaHistoryIntact,2)
%     figure
%     histfit(R(:,i))
%     title(['No. ',num2str(i),' Element'])
    pd{1,i} = fitdist(R(:,i),'Normal');
    ThetaIntact_std(i) = pd{1,i}.sigma;
end

% Error bar
NoSamples = 1;
CI95_Gaussian = 1.96;
HalfDelta_CI = CI95_Gaussian*ThetaIntact_std/sqrt(NoSamples);

figure('Position',[488,342,560*2,420*0.85])
p = bar([(1:NoEle)',(1:NoEle)'], [-ThetaIntact', -ThetaHistoryIntact(end,:)']);            
p(1).EdgeColor = 'none';
p(2).EdgeColor = 'none';
p(1).FaceColor = 'r';
p(2).FaceColor = [0.07,0.62,1.00];
hold on
% Plot error bars on predicted values
ngroups = length(ThetaIntact);
nbars = 2;
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 2
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    er = errorbar(x, -ThetaHistoryIntact(end,:), HalfDelta_CI,'LineWidth',1);
end
er.Color = [0 0 0];                            
er.LineStyle = 'none'; 
title('Model Updating','Interpreter','latex')
legend([p(:,1),p(:,2)],{'The True','Bayes $\ell_2$'},'Interpreter','latex','Orientation','horizontal','Location','north')
legend('boxoff')
xlabel('Element \#','Interpreter','latex')
ylabel('\boldmath$\theta_{intact}$','Interpreter','latex')
set(gca,'TickLabelInterpreter','latex','FontSize',19)
box off
% ylim([-0.3 0.3])
% yticks([-0.3:0.15:0.3])

%% Damage identification
% ReguFlag = 1: Sparse Bayesian learning (L1); ReguFlag = 2:
% Bayesian learning (L2); ReguFlag = 3: L curve; ReguFlag = 4: STLS.
disp(' ')
disp('Damage identification begins')
ReguFlag = 4; 

ThetaHistoryIntact = zeros(1,NoEle);

[ThetaHistoryDamaged,~,~] = Regu_eig(N,ReguFlag,SensorInd,Fs,ThetaHistoryIntact(end,:),...
    FreqDamaged(ModeInd2),ShapeDamaged(:,ModeInd2),ModeInd2,NoDOF);

% load('D:\Dropbox\Seismic Interferometry, Sparse Bayes & L0 Regu\Plane truss code\Figures for L0 paper\DI_STLS.mat')
% Result visualization

toc

figure('Position',[488,342,560*2,420*0.85])
p = bar([(1:NoEle)',(1:NoEle)'], [-ThetaDamaged', -ThetaHistoryDamaged(end,:)']);        
p(1).EdgeColor = 'none';
p(2).EdgeColor = 'none';
p(1).FaceColor = 'r';
p(2).FaceColor = [0.07,0.62,1.00];
title('Sparse Damage Identification','Interpreter','latex')
legend([p(:,1),p(:,2)],{'The True','STLS'},'Interpreter','latex','Orientation','horizontal','Location','south')
legend('boxoff')
xlabel('Element \#','Interpreter','latex')
ylabel('\boldmath$\theta_{dmg}$','Interpreter','latex')
set(gca,'TickLabelInterpreter','latex','FontSize',19)
box off