function [x,varargout] = BayesL2_eig(Phi,b,varargin)

NoEle = size(Phi,2);
Nob = size(b,1); % NoN = No. of observations * No. of timesteps
k = 1;
MaxIter = 100;

% Adaptive determination of hyper parameters
if nargin == 5
    if double(varargin{1}) == 1
        alpha0 = 1e-2; beta0 = 1e-5;
        alpha1 = 1e-6; beta1 = 1e-9;
        Eta = NaN(1, MaxIter); Eta(k) = 1e-4;
        varargout{1} = beta0;
        varargout{2} = beta1;
    else
        xx = double(varargin{2});
        Ax_b = double(varargin{3});
        alpha0 = 1e-3; beta0 = 1e-3*0.5*Ax_b'*Ax_b;
        alpha1 = 1e-5; beta1 = 1e-3*0.5*xx*xx';
        Eta = NaN(1, MaxIter); Eta(k) = 1e-2;
        varargout{1} = beta0;
        varargout{2} = beta1;
    end    
end

% Good setup for 5% noise case w 10 sensors
% if nargin == 5
%     if double(varargin{1}) == 1
%         alpha0 = 1e-4; beta0 = 1e-7;
%         alpha1 = 1e-2; beta1 = 1e-5;
%         Eta = NaN(1, MaxIter); Eta(k) = 1e-2;
%         varargout{1} = beta0;
%         varargout{2} = beta1;
%     else
%         xx = double(varargin{2});
%         Ax_b = double(varargin{3});
%         alpha0 = 1e-4; beta0 = 1e-3*0.5*Ax_b'*Ax_b;
%         alpha1 = 1e-2; beta1 = 1e-3*0.5*xx*xx';
%         Eta = NaN(1, MaxIter); Eta(k) = 1e-2;
%         varargout{1} = beta0;
%         varargout{2} = beta1;
%     end    
% end

% Good setup for 0 noise case w 10 sensors
% if nargin == 5
%     if double(varargin{1}) == 1
%         alpha0 = 1e-2; beta0 = 1e-5;
%         alpha1 = 1e-2; beta1 = 1e-5;
%         Eta = NaN(1, MaxIter); Eta(k) = 1e0;
%         varargout{1} = beta0;
%         varargout{2} = beta1;
%     else
%         xx = double(varargin{2});
%         Ax_b = double(varargin{3});
%         alpha0 = 1e-2; beta0 = 1e-3*0.5*Ax_b'*Ax_b;
%         alpha1 = 1e-2; beta1 = 1e-3*0.5*xx*xx';
%         Eta = NaN(1, MaxIter); Eta(k) = 1e0;
%         varargout{1} = beta0;
%         varargout{2} = beta1;
%     end    
% end

II = eye(NoEle);
Tolerance = 1e-3;
EtaErr = 1;
ParamErr = 1;
Param = zeros(NoEle, 1);
Phi2 = Phi';
Phi1 = Phi2*Phi;
while k < MaxIter & EtaErr > Tolerance & ParamErr > Tolerance
%         fprintf(1,'Inverse Analysis \n');
    k = k + 1;
    Param(:, k) = (Phi1 + Eta(k-1)*II)\(Phi2*b);
    Temp = Phi*Param(:, k) - b;
    sigmaSQR = (Temp'*Temp+2*beta0)/(Nob+2*(alpha0+1));
    lambdaSQR = (Param(:, k)'*Param(:, k)+2*beta1)/(NoEle+2*(alpha1+1));
    Eta(k) = sigmaSQR/lambdaSQR;
    EtaErr = abs((Eta(k) - Eta(k-1))/Eta(k));
    ParamErr = norm(Param(:, k) - Param(:, k-1))/norm(Param(:, k));
end
x = Param(:, end)';
hold on
plot(Eta)
set(gca,'yscale','log')
title('Eta in inner loop')
varargout{3} = Eta(k);
varargout{4} = sqrt(lambdaSQR);

NoSamples = 1;
%% Approximate variance by assuming the marginal posterior is a Gaussian distribution
% Refer to KV Yuen's 2011 book
x_col = x';
% Dx = diag(0.05*(x_col));
Dx = diag(ones(size(x_col)))*1e-3;
for i = 1:length(x_col)
    for j = i:length(x_col)
        if i == j
            % Diagonal elements
            J1 = 1/2/lambdaSQR*(x_col+Dx(:,i))'*(x_col+Dx(:,i)); % J(x+dx)
            J2 = 1/2/lambdaSQR*x_col'*x_col;% J(x)
            J3 = 1/2/lambdaSQR*(x_col-Dx(:,i))'*(x_col-Dx(:,i));% J(x-dx)
            for k = 1:NoSamples
                J1 = J1+1/2/sigmaSQR*(Phi(:,:,k)*(x_col+Dx(:,i))-b(:,:,k))'*(Phi(:,:,k)*(x_col+Dx(:,i))-b(:,:,k));        
                J2 = J2+1/2/sigmaSQR*(Phi(:,:,k)*x_col-b(:,:,k))'*(Phi(:,:,k)*x_col-b(:,:,k));
                J3 = J3+1/2/sigmaSQR*(Phi(:,:,k)*(x_col-Dx(:,i))-b(:,:,k))'*(Phi(:,:,k)*(x_col-Dx(:,i))-b(:,:,k));
            end
            H(i,i) = (J1-2*J2+J3)/(0.01*x_col(i))^2;
        else
            % Off-diagonal elements
            J1 = 1/2/lambdaSQR*(x_col+Dx(:,i)+Dx(:,j))'*(x_col+Dx(:,i)+Dx(:,j)); % J(x_col+Dx(:,i)+Dx(:,j))
            J2 = 1/2/lambdaSQR*(x_col+Dx(:,i)-Dx(:,j))'*(x_col+Dx(:,i)-Dx(:,j)); % J(x_col+Dx(:,i)-Dx(:,j))
            J3 = 1/2/lambdaSQR*(x_col-Dx(:,i)+Dx(:,j))'*(x_col-Dx(:,i)+Dx(:,j));% J(x_col-Dx(:,i)+Dx(:,j))
            J4 = 1/2/lambdaSQR*(x_col-Dx(:,i)-Dx(:,j))'*(x_col-Dx(:,i)-Dx(:,j));% J(x_col-Dx(:,i)-Dx(:,j))
            for k = 1:NoSamples
                J1 = J1+1/2/sigmaSQR*(Phi(:,:,k)*(x_col+Dx(:,i)+Dx(:,j))-b(:,:,k))'*(Phi(:,:,k)*(x_col+Dx(:,i)+Dx(:,j))-b(:,:,k));        
                J2 = J2+1/2/sigmaSQR*(Phi(:,:,k)*(x_col+Dx(:,i)-Dx(:,j))-b(:,:,k))'*(Phi(:,:,k)*(x_col+Dx(:,i)-Dx(:,j))-b(:,:,k));
                J3 = J3+1/2/sigmaSQR*(Phi(:,:,k)*(x_col-Dx(:,i)+Dx(:,j))-b(:,:,k))'*(Phi(:,:,k)*(x_col-Dx(:,i)+Dx(:,j))-b(:,:,k));
                J4 = J4+1/2/sigmaSQR*(Phi(:,:,k)*(x_col-Dx(:,i)-Dx(:,j))-b(:,:,k))'*(Phi(:,:,k)*(x_col-Dx(:,i)-Dx(:,j))-b(:,:,k));
            end
            H(i,j) = (J1-J2-J3+J4)/4/(0.01*x_col(i))/(0.01*x_col(j));
        end
    end    
end
H = (H+H') - eye(size(H,1)).*diag(H);
x_cov = inv(H/NoSamples);
varargout{5} = x_cov;