function [thetaNew,Lambda] = SBLCallFn(Phi,deltaX)
    % This function carries out sparse Bayesian learning as well as
    % plotting the results.
    
    %% Tipping's code
    addpath('..\SB2_Release_200')
    LIKELIHOOD = 'Gaussian';
    BASIS = Phi;
    TARGETS = deltaX;    
           
    % updated noise
    OPTIONS = SB2_UserOptions('MONITOR',1,'DIAGNOSTICLEVEL',1);
    [PARAMETER, HYPERPARAMETER, ~] = SparseBayes(LIKELIHOOD, BASIS, TARGETS,OPTIONS);
    
    thetaNew = zeros(1,size(Phi,2));
    thetaNew(PARAMETER.Relevant) = PARAMETER.Value;
    Lambda = zeros(1,size(Phi,2));
    Lambda(PARAMETER.Relevant) = sqrt(1./(HYPERPARAMETER.Alpha));
end