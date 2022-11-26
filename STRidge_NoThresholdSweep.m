function DeltaGamma = STRidge_NoThresholdSweep(Phi, Residue, knob)
% %% The outmost loop controls cross validation
% for i = 1:3      
%     num_train = round(size(Phi, 1)*0.8);
%     train_ind = randperm(size(Phi, 1), num_train);
%     val_ind = setdiff(1:size(Phi, 1), train_ind);
%     Phi_train = Phi(train_ind, :);
%     Residue_train = Residue(train_ind);
%     Phi_val = Phi(val_ind, :);
%     Residue_val = Residue(val_ind);
%     
%     %% Normalize Phi
%     for ell = 1:size(Phi_train, 2)
%         Mreg_train(ell,1) = 1/norm(Phi_train(:, ell));
%         Phi_train_norm(:, ell) = Phi_train(:, ell)*Mreg_train(ell);
%     end
%     
%     DeltaGamma_norm = Phi_train_norm\Residue_train;  % initial guess: Least-squares
%     error_val_best = norm(Phi_val*(DeltaGamma_norm.*Mreg_train)-Residue_val) + 0.001*cond(Phi)*nnz(DeltaGamma_norm);
%     
%     for k=1:10
%         DeltaGamma_norm_new = DeltaGamma_norm*0;
%         biginds = (abs(DeltaGamma_norm.*Mreg_train)>=knob);   % find big coefficients
%         % Regress dynamics onto remaining terms to find sparse DeltaGamma_norm
%         DeltaGamma_norm_new(biginds) = Phi_train_norm(:,biginds)\Residue_train; 
%         
%         DeltaGamma_new = DeltaGamma_norm_new.*Mreg_train;
%         error_val = norm(Phi_val*DeltaGamma_new-Residue_val) + 0.001*cond(Phi)*nnz(DeltaGamma_new);
%         
%         if error_val < error_val_best
%             
%             if nnz(DeltaGamma_new) == 0
%                 break
%             end
%             
%             error_val_best = error_val;
%             DeltaGamma_norm = DeltaGamma_norm_new;
%         else
%             break
%         end
%     end
%     DeltaGamma0(:, i) = DeltaGamma_norm.*Mreg_train;
%     
%     %% W/O normalizing Phi    
% %     DeltaGamma = Phi_train\Residue_train;  % initial guess: Least-squares
% %     error_val_best = norm(Phi_val*DeltaGamma-Residue_val) + cond(Phi)*nnz(DeltaGamma);
% % 
% %     for k = 1:10
% %         DeltaGamma_new = DeltaGamma*0;
% %         biginds = (abs(DeltaGamma)>=knob);   % find big coefficients
% %         % Regress dynamics onto remaining terms to find sparse DeltaGamma_norm
% %         DeltaGamma_new(biginds) = Phi_train(:,biginds)\Residue_train; 
% %         error_val = norm(Phi_val*DeltaGamma_new-Residue_val) + cond(Phi)*nnz(DeltaGamma_new);
% %         
% %         if error_val < error_val_best
% %             
% %             if nnz(DeltaGamma_new) == 0
% %                 break
% %             end
% %             
% %             error_val_best = error_val;
% %             DeltaGamma = DeltaGamma_new;
% %         else
% %             break
% %         end
% %     end
% %     DeltaGamma0(:, i) = DeltaGamma;
%     
% end
% DeltaGamma = mean(DeltaGamma0, 2);

%%% W/O cross-validation
%% Normalize Phi
% for ell = 1:size(Phi, 2)
%     Mreg(ell,1) = 1/norm(Phi(:, ell));
%     Phi_norm(:, ell) = Phi(:, ell)*Mreg(ell);
% end
% 
% DeltaGamma_norm = Phi_norm\Residue;  % initial guess: Least-squares
% error_best = norm(Phi*(DeltaGamma_norm.*Mreg)-Residue) + 0.001*cond(Phi)*nnz(DeltaGamma_norm);
% 
% for k=1:10
%     DeltaGamma_norm_new = DeltaGamma_norm*0;
%     biginds = (abs(DeltaGamma_norm.*Mreg)>=knob);   % find big coefficients
%     % Regress dynamics onto remaining terms to find sparse DeltaGamma_norm
%     DeltaGamma_norm_new(biginds) = Phi_norm(:,biginds)\Residue; 
% 
%     DeltaGamma_new = DeltaGamma_norm_new.*Mreg;
%     error = norm(Phi*DeltaGamma_new-Residue) + 0.001*cond(Phi)*nnz(DeltaGamma_new);
% 
%     if error < error_best
% 
%         if nnz(DeltaGamma_new) == 0
%             break
%         end
% 
%         error_best = error;
%         DeltaGamma_norm = DeltaGamma_norm_new;
%     else
%         break
%     end
% end
% DeltaGamma = DeltaGamma_norm.*Mreg;

%% Lasso
% W/O normalizing Phi(Because Sparse Bayesian Learning doesn't normalize Phi and still has good results)
[B,FitInfo] = lasso(Phi, Residue, 'CV', 10);
idxDeltaGamma1SE = FitInfo.Index1SE;
DeltaGamma = B(:, idxDeltaGamma1SE); % Initial guess
error_best = norm(Phi*DeltaGamma-Residue) + 0.001*cond(Phi)*nnz(DeltaGamma);
for k=1:10
    DeltaGamma_new = DeltaGamma*0;
    biginds = (abs(DeltaGamma)>=knob);   % find big coefficients
    % Regress dynamics onto remaining terms to find sparse DeltaGamma
    DeltaGamma_new(biginds) = Phi(:,biginds)\Residue;     
    
    error = norm(Phi*DeltaGamma_new-Residue) + 0.001*cond(Phi)*nnz(DeltaGamma_new);

    if error < error_best
        if nnz(DeltaGamma_new) == 0 % if all coefficients are non-zero meanwhile the least-square error is smaller the the lasso error
            break
        else
            error_best = error;
            DeltaGamma = DeltaGamma_new;
        end
    else
        break
    end
end