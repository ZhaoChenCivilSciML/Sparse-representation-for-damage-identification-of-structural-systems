function Phi = SensitivityMatrix_eig(weig,wshape,eigvalFEMA,eigvalMea,shapeFEMA,S,...
    NoEle,NoDOF,NoMode,SensorInd)
KShape = TrussShapeMatrix; % KShape uses the original stiffness values.
for i = 1:NoEle
    for j = 1:NoDOF
        Phi1(j,i) = shapeFEMA(:,j)'*(-KShape(:,:,i))*shapeFEMA(:,j);
        
        % Mode shape derivative: Type 1
%         F = S.K-eigvalFEMA(j)*S.M;
%         dF = -KShape(:,:,i)-Phi1(j,i)*S.M;
%         Phi2(1+(j-1)*NoDOF:j*NoDOF,i) = -inv(F*F+2*S.M*shapeFEMA(:,j)*shapeFEMA(:,j)'*S.M)*F*dF*shapeFEMA(:,j);
        
        % Mode shape derivative: Type 2        
        Temp = zeros(size(shapeFEMA(:,1)));
        for h = 1:NoDOF
            if h ~= j
                a = shapeFEMA(:,h)'*(-KShape(:,:,i))*shapeFEMA(:,j)/...
                    (eigvalFEMA(j)-eigvalFEMA(h));
            else
                a = 0;
            end
            Temp = Temp+a*shapeFEMA(:,h);            
        end
        Phi2(1+(j-1)*NoDOF:j*NoDOF,i) = Temp;
        
    end
end

% Freqs
Phi(1:NoMode,:) = diag(weig./(eigvalMea))*Phi1(1:NoMode,:);

% Shapes
R = zeros(length(SensorInd),NoDOF);
for i = 1:length(SensorInd)
    R(i,SensorInd(i)) = 1;
end
for i = 1:NoMode
    Phi(NoMode+1+(i-1)*length(SensorInd):NoMode+i*length(SensorInd),:) = R*Phi2(1+(i-1)*NoDOF:i*NoDOF,:);
end
Phi(NoMode+1:NoMode+length(SensorInd)*NoMode,:) = diag(wshape)*Phi(NoMode+1:NoMode+length(SensorInd)*NoMode,:);