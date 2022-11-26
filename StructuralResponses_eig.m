function [ag,acc,vel,dis,S] = StructuralResponses_eig(N,NoiseLevel,SensorInd,...
    Theta,Fs,ModeInd,varargin)
%% Define structural input
% White Noise
rng('default')
ag0 = rand(1,N);
% Filter the Excitation
% Zero-Pole-Gain design
[z,p,k] = butter(5, [0.1/(Fs/2) 600/(Fs/2)],'bandpass');
[sos,g] = zp2sos(z,p,k);
ag0 = sosfilt(sos,ag0)*g;

%% Define material properties
Xi1 = 0.01; % damping ratio
Xi2 = 0.02;

%% Define mass, stiffness and damping matrices.
NoEle = 31;
NoNode = 14;
NoDOF = 2*NoNode;
[KShape,MShape] = TrussShapeMatrix;
M = zeros(NoDOF);
K = zeros(NoDOF);
% K0 = zeros(NoDOF);
for i = 1:NoEle
    M = M+MShape(:,:,i);
    K = K+KShape(:,:,i)*(1-Theta(i));
end
% Add restraints
K = K+KShape(:,:,end);

% K = K(3:27,3:27);
% M = M(3:27,3:27);

[shape, omega0] = eig(K, M);
if nargin == 6
    Omega0 = sqrt(diag(omega0));
    ww1 = Omega0(1); 
    ww2 = Omega0(2);     
elseif nargin == 7
    freq = double(varargin{1});
    ww1 = freq(1)*2*pi;
    ww2 = freq(2)*2*pi;
end
beta = 2*(ww1*Xi1 - ww2*Xi2)/(ww1^2-ww2^2); 
alpha = 2*ww1*Xi1 - beta*ww1^2;
C = alpha*M + beta*K; 
S.M = M;
S.C = C;
S.K = K;
S.alpha = alpha;
S.beta = beta;
S.normalshape = shape(SensorInd,ModeInd);
temp = sqrt(diag(omega0))/2/pi;
S.freq = temp(ModeInd);
S.CNormalShape = shape;
S.Cfreq = sqrt(diag(omega0))/2/pi;

%% Compute acc, vel and dis
RR = zeros(length(SensorInd), size(M,1));
for ii = 1:size(RR, 1)
   RR(ii, SensorInd(ii)) = 1; 
end
L = zeros(size(M,1),1);

% Point Load
L(10) = 1;
L(13) = 1;
A = [zeros(size(M,1)), eye(size(M,1)); -(M\K), -(M\C)];                                     % A matrix 
B = [zeros(size(M,1),1); M\L];                                               % B matrix       
CC = [-M\K, -M\C]; 
D = M\L; 

sys = c2d(ss(A, B, CC, D), 1/Fs);
[acc0, ~, x] = lsim(sys, ag0,0:(1/Fs):(N-1)/Fs);  % The total acceleration is usually measured by an accelerometer.

acc1 = RR*acc0';
vel0 = (x(:,size(M,1)+1:end))'; % The relative velocity
dis0 = (x(:,1:size(M,1)))'; % The relative displacement

%% Add noise
acc = addnoise(acc1, NoiseLevel);
vel = addnoise(vel0, NoiseLevel);
dis = addnoise(dis0, NoiseLevel);

ag = addnoise(ag0, NoiseLevel);
