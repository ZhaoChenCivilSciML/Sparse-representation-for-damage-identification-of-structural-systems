function [ThetaHistory,Lambda, ThetaCovHistory_approx] = Regu_eig(N,ReguFlag,SensorInd,Fs,Theta,FreqMea,shapeMea,...
    ModeInd,NoDOF)
tic
NoEle = length(Theta);
ThetaHistory = zeros(1,NoEle);
knob_opt_history = 0;
Tol = 1e-6;
MaxIter = 20;
MaxIter = 20;
i = 1;
ThetaErr = 1;
DeltaTheta = zeros(1,NoEle);
Ax_b = 0;
ObjFunc = 0;
ObjFuncPt1 = 0;
ObjFuncPt2 = 0;
beta0 = 0;
beta1 = 0;
Eta = 0;
NormDeltaTheta = 0;
ThetaCovHistory_approx = zeros(NoEle);
while ThetaErr(i) > Tol && (i-1) < MaxIter
% while (i-1) < MaxIter
    fprintf(1, 'Iteration = %d; Iterative update = %f\n', [i,ThetaErr(i)]);
    ThetaCompo = Theta+ThetaHistory(i,:)-Theta.*ThetaHistory(i,:);
    [~,~,~,~,S] = StructuralResponses_eig(N,0,SensorInd,ThetaCompo,Fs,ModeInd,FreqMea);
    FreqFEM = S.freq; % unit in Hz
    shapeFEM = S.normalshape;
    eigvalMea = (FreqMea*2*pi).^2;
    eigvalFEM = (FreqFEM*2*pi).^2;
    if size(eigvalMea,1)<size(eigvalMea,2)
        eigvalMea = eigvalMea';
    end
    if size(eigvalFEM,1)<size(eigvalFEM,2)
        eigvalFEM = eigvalFEM';
    end
    %% Residue 
    % Type 1
    [shapeFEM2,eigvalFEM2,shapeMea2] = MACpair(shapeFEM,eigvalFEM,shapeMea);
    
    % Type 2
%     shapeFEM2 = shapeFEM;
%     eigvalFEM2 = eigvalFEM;
    
    % This normalization is feasible when sensors are only deployed on
    % selected DOFs.
    for j = 1:size(shapeMea2,2)
        shapeMea3(:,j) = shapeMea2(:,j)*norm(shapeFEM2(:,j))/norm(shapeMea2(:,j));
    end    
    ResidueEigval = ones(size(eigvalFEM2))-eigvalFEM2./eigvalMea;
    ResidueShape = shapeMea3-shapeFEM2;
    ResidueShape2 = reshape(ResidueShape,[],1);
    
    % Residue weight
    weig = ones(length(ModeInd),1);
    ResidueEigval2 = diag(weig)*ResidueEigval;
    wshape = 1/max(max(shapeMea3))*ones(length(ModeInd)*length(SensorInd),1);
    
%     wshape(end-11:end) = 2;
    
    ResidueShape3 = diag(wshape)*ResidueShape2;
    Residue = [ResidueEigval2;ResidueShape3];
    %% Sensitivity matrix
    eigvalFEMA = (S.Cfreq*2*pi).^2; % FEM eigenvalues at All DOFs
    Phi = SensitivityMatrix_eig(weig,wshape,eigvalFEMA,eigvalMea,S.CNormalShape,S,...
    NoEle,NoDOF,length(ModeInd),SensorInd);
    switch ReguFlag
        case 1
            [DeltaTheta,Lambda(i+1,:)] = SBLCallFn(Phi,Residue);
            ThetaHistory(i+1,:) = ThetaHistory(i,:)+DeltaTheta;
            ThetaErr(i+1) = norm(ThetaHistory(i+1,:)-ThetaHistory(i,:))/norm(ThetaHistory(i+1,:));
        case 2
            [DeltaTheta,beta0(i+1),beta1(i+1),Eta(i+1),Lambda, DeltaTheta_approx_cov] = BayesL2_eig(Phi,Residue,i,DeltaTheta,Ax_b);
%             [DeltaTheta,beta0(i+1),beta1(i+1),Eta(i+1)] = ConsistentBayesL2(Phi,Residue,i,DeltaTheta,Ax_b,ThetaHistory(i,:)');
            ThetaHistory(i+1,:) = ThetaHistory(i,:)+DeltaTheta;
            Ax_b = Phi*DeltaTheta'-Residue;                            
            ObjFuncPt1(i+1) = (norm(Phi*DeltaTheta'-Residue))^2;
            NormDeltaTheta(i+1) = (norm(DeltaTheta))^2;
            ObjFuncPt2(i+1) = Eta(i+1)*NormDeltaTheta(i+1);
            ObjFunc(i+1) = ObjFuncPt1(i+1)+ObjFuncPt2(i+1); % The objective function of Tik Regularization.
            
            ThetaErr(i+1) = norm(ThetaHistory(i+1,:)-ThetaHistory(i,:))/norm(ThetaHistory(i+1,:)); 
            
            ThetaCovHistory_approx(:,:,i+1) = ThetaCovHistory_approx(:,:,i) + DeltaTheta_approx_cov;
            
            % Use the objective function of L2 regularization as a
            % convergence criterion.            
%             ThetaErr(i+1) = abs((ObjFunc(i+1)-ObjFunc(i))/ObjFunc(i+1)); 
        case 3
            addpath('..\Tihkonov regularization\regu')
            [U,s,V] = csvd(Phi);
            reg_para = l_curve(U,s,Residue);
%             reg_para = gcv(U,s,Deltaacc); 
%             reg_para = quasiopt(U,s,Deltaacc);
            DeltaTheta = tikhonov(U,s,V,Residue,reg_para);
            if size(DeltaTheta,1) > 1
                DeltaTheta = DeltaTheta';
            end
            ThetaHistory(i+1,:) = ThetaHistory(i,:)+DeltaTheta;
            % Norm convergence
            ThetaErr(i+1) = norm(ThetaHistory(i+1,:)-ThetaHistory(i,:))/...
                norm(ThetaHistory(i+1,:));
        % Angle convergence
%             A = [Phi;reg_para*eye(size(Phi,2))];
%             r = [Deltaacc;zeros(size(Phi,2),1)];
%             ThetaErr(i+1) = r'*(A*inv(A'*A)*A'*r)/norm(r)/...
%                 norm(A*(inv(A'*A)*A'*r));
            Lambda = 0;
        case 4
            % Tuning the sparsity knob via Bayesian optimization
            knob_bayes = optimizableVariable('knob', [0.01 1]);
            fun = @(knob_bayes)STRidge_BayesOpt(Phi, Residue, knob_bayes);
            bayes_results = bayesopt(fun,knob_bayes,'AcquisitionFunctionName','expected-improvement-plus','UseParallel',true,'PlotFcn',[],'Verbose',0);
            knob_opt = bayes_results.XAtMinObjective.(1)(1);           
            
%             knob_opt = 0.1;

             knob_opt_history(i+1) = knob_opt;
            DeltaTheta = STRidge_NoThresholdSweep(Phi, Residue, knob_opt);
                       
            Lambda = 0;
            ThetaHistory(i+1,:) = ThetaHistory(i,:)+DeltaTheta';
            ThetaErr(i+1) = norm(ThetaHistory(i+1,:)-ThetaHistory(i,:))/norm(ThetaHistory(i+1,:));
        case 5
            [B,FitInfo] = lasso(Phi, Residue, 'CV', 10);
            idxDeltaThetaMinMSE = FitInfo.IndexMinMSE;
            DeltaTheta = B(:, idxDeltaThetaMinMSE); % Initial guess
            Lambda = 0;
            ThetaHistory(i+1,:) = ThetaHistory(i,:)+DeltaTheta';
            ThetaErr(i+1) = norm(ThetaHistory(i+1,:)-ThetaHistory(i,:))/norm(ThetaHistory(i+1,:));
    end    
    i = i+1;
end

% figure
% plot(eigvalFEM2)
% hold on
% plot(eigvalMea)
% legend('eigvalFEM','eigvalMea')
% for i = 1:size(shapeMea3,2)
%     figure
%     plot(shapeFEM2(:,i))
%     hold on
%     plot(shapeMea3(:,i))
%     legend('shapeFEM','shapeMea')
%     title(['Mode Shape',num2str(i)])
% end

% figure
% hold on
% for i = 1:size(ThetaHistory,2)
%     plot(-ThetaHistory(:,i),'LineWidth',2)
% end
% hold off
% % Model updating
% title('Convergence of \boldmath$\theta_{intact}$','Interpreter','latex')
% xlabel('Iteration \#','Interpreter','latex')
% ylabel('\boldmath$\theta_{intact}$','Interpreter','latex')
% set(gca,'TickLabelInterpreter','latex','FontSize',19)
% box off

% Damage identification
% title('Convergence of \boldmath$\theta_{dmg}$','Interpreter','latex')
% xlabel('Iteration \#','Interpreter','latex')
% ylabel('\boldmath$\theta_{dmg}$','Interpreter','latex')
% set(gca,'TickLabelInterpreter','latex','FontSize',19)
% box off
% annotation(gcf,'textbox',...
%     [0.709999999999999,0.29714285714286,0.171428571428571,0.057142857142857],...
%     'String','Ele \#1',...
%     'Interpreter','latex',...
%     'HorizontalAlignment','center',...
%     'FontSize',15,...
%     'FitBoxToText','off',...
%     'EdgeColor','none');
% annotation(gcf,'textbox',...
%     [0.694999999999999,0.423809523809527,0.171428571428571,0.057142857142857],...
%     'String','Ele \#27',...
%     'Interpreter','latex',...
%     'HorizontalAlignment','center',...
%     'FontSize',15,...
%     'FitBoxToText','off',...
%     'EdgeColor','none');
% annotation(gcf,'textbox',...
%     [0.695714285714286,0.481428571428572,0.171428571428571,0.057142857142857],...
%     'String','Ele \#15',...
%     'Interpreter','latex',...
%     'HorizontalAlignment','center',...
%     'FontSize',15,...
%     'FitBoxToText','off',...
%     'EdgeColor','none');

% figure
% plot(ThetaErr)
% set(gca,'yscale','log')
% title('Iterative update')
% 
% if ReguFlag == 2
%     figure
%     plot(beta0)
%     title('beta0')
%     set(gca,'yscale','log')
%     figure
%     plot(beta1)
%     title('beta1')
%     set(gca,'yscale','log')    
%     figure
%     semilogy(ObjFunc)
%     title('Objective Function')
%     set(gca,'yscale','log')
%     figure
%     plot(ObjFuncPt1)
%     title('Objective Function Pt1')
%     figure
%     plot(ObjFuncPt2)
%     title('Objective Function Pt2')
%     set(gca,'yscale','log')
%     figure
%     plot(Eta)
%     title('Eta')
%     set(gca,'yscale','log')
%     figure
%     plot(NormDeltaTheta)
%     title('The squared L2 norm of DeltaTheta')
%     set(gca,'yscale','log')
% end
toc

% if ReguFlag == 4
%     figure
%     plot(knob_opt_history)
%     title('knob_opt_history')
% end