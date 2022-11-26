function [KShape,MShape] = TrussShapeMatrix
% SShape: Stiffness shape
% MShape: Mass shape
Rho = 2770;
Area = .0025;
Le = [1.52 1.52 2.15 2.15 1.52 1.52 1.52 2.15 2.15 1.52 1.52 1.52 2.15 2.15 ...
    1.52 1.52 1.52 2.15 2.15 1.52 1.52 1.52 2.15 2.15 1.52 1.52 1.52 2.15 ...
    2.15 1.52 1.52]; % element lengths
E = 7e10;
NoDOF = 28;
NoEle = 31;
% Inclination angle. (The angle between a line and the x-axis)
InclAngle = [pi/2 0 -pi/4 pi/4 0 pi/2 0 -pi/4 pi/4 0 pi/2 0 -pi/4 pi/4 0 ...
    pi/2 0 -pi/4 pi/4 0 pi/2 0 -pi/4 pi/4 0 pi/2 0 -pi/4 pi/4 0 pi/2]; 
KShape = zeros(NoDOF,NoDOF,NoEle);
MShape = zeros(NoDOF,NoDOF,NoEle);
for i = 1:NoEle
    T = [cos(InclAngle(i)) sin(InclAngle(i)) 0 0;
        -sin(InclAngle(i)) cos(InclAngle(i)) 0 0;        
        0 0 cos(InclAngle(i)) sin(InclAngle(i));
        0 0 -sin(InclAngle(i)) cos(InclAngle(i))];   
    DOFs = FEMassemble(i);    
    KShape(DOFs,DOFs,i) = E*Area/Le(i)*T'*[1 0 -1 0;0 0 0 0;-1 0 1 0;0 0 0 0]*T;
    MShape(DOFs,DOFs,i) = Area*Rho*Le(i)/6*T'*[2 0 1 0;0 2 0 1;1 0 2 0;0 1 0 2]*T; % consistent mass
%     MShape(DOFs,DOFs,i) = Rho*Area*Le(i)/2*T'*eye(4)*T; % lumped mass
end
KShape([1,2,28],[1,2,28],NoEle+1) = 1e10*eye(3); % stiffness at supports