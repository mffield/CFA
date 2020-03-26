function [g, mapping]=CFA(X, opt)
%CFA Performs manifold learning on dataset X 
%
%   [g, mapping]=CFA(X, opt)

% Performs manifold learning on dataset X to reduce the dimensionality.
% The variable no_analyzers determines the number of
% local factor analyzers that is used in the mixture of factor analyzers.
% The variable opt is a structure of settings including the maximum iterations, 
% tolerances and plotting options. Please refer to the article for more details.
% A plotting is included as option within this code for demonstration purposes.

% 
% Written by Matthew Field

% This algorithm is a reduced version of methods published in:
% *M. Field, D. Stirling, Z. Pan and F. Naghdy, "Learning Trajectories for Robot
% Programing by Demonstration Using a Coordinated Mixture of Factor Analyzers,"
% in IEEE Transactions on Cybernetics, vol. 46, no. 3, pp. 706-717, March 2016.
% https://doi.org/10.1109/TCYB.2015.2414277


if ~isfield(opt, 'd')
    d = 2; % set default number of latent dimensions.
else d = opt.d;
end
if ~isfield(opt, 'K')
    c = 20; % set default number of clusters.
else c = opt.K;
end
if ~isfield(opt, 'max_iterations')
    opt.max_iterations = 100; % set default maximum iterations.
end
if ~isfield(opt, 'lle_nn')
    opt.lle_nn = 20; % set default number of nearest neighbours in LLE initialiser.
end
if ~isfield(opt, 'lle_iterations')
    opt.lle_iterations = 10; % set default number of iterations to fit parameters to LLE initial coordinates.
end
if ~isfield(opt, 'verbose')
    opt.verbose = 0; % set flag for additional output.
end
if ~isfield(opt, 'plot')
    opt.plot = 0; % set flag for plotting data and parameters in 3D->2D examples.
end
if ~isfield(opt, 'InitPCA')
    opt.InitPCA = 0; % set the percent variance retained for initial PCA. Set to zero if no initial PCA required.
end

% Make sure data has zero-mean.
    X_mu = mean(X, 1);
    X = X - repmat(X_mu, [size(X, 1) 1]);
    
    [UU1,SS,~]=svd(cov(X)); 
    Xpca = X*UU1(:,1:d);
    Xhat = (UU1(:,1:d)*Xpca')';
    PCA_error = sqrt(sum(sum((X-Xhat).^2)));
    
    X_o = X; 
    dPCA=d;
    if opt.InitPCA
        dPCA = find(cumsum(diag(SS))/sum(diag(SS))>opt.InitPCA,1,'first');
        X = X*UU1(:,1:dPCA);
        
        Xpca = X_o*UU1(:,1:d);
        Xhat = (UU1(:,1:d)*Xpca')';
        PCA_error = sqrt(sum(sum((X_o-Xhat).^2)));
    end
    X = X';
    [D,n] = size(X);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Initialize the parameters of the CFA model
    
    SigmaN = repmat(eye(d), [1 1 n]);
    Pi = (1/c)*ones(1,c);
    Kappa = rand(d, c) * .01;
    SigmaC = repmat(eye(d), [1 1 c]);
    Lambda = zeros(D, d, c);
    Psi = zeros(D, D, c);
    C = zeros(D,D,c);
    R = zeros(d,d,c); 
    Proj = zeros(D,d,c);
    Q=zeros(n,c);
    logdetSigmaN = zeros(n,1)+realmin;
    %%% Initialize with a mixture of Gaussians
    obj=gmdistribution.fit(X',c,'Regularize',0.000001);
    Mu = obj.mu';
    for k=1:c
        [U(:,:,k),EigV,~]=svd(obj.Sigma(:,:,k));
        EigVal(:,k) = diag(EigV);
        sigma2(k) = (1/(D-d))*sum(EigVal(d+1:D,k));
        Psi(:,:,k) = sigma2(k)*eye(D);
        Lambda(:,:,k) = U(:,1:d,k)*sqrt(diag(EigVal(1:d,k))-sigma2(k)*eye(d));
        C(:,:,k)=Lambda(:,:,k)*diag(EigVal(1:d,k))*Lambda(:,:,k)' + Psi(:,:,k);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Initialize the parameters of the CFA model using LLE latent space positions as complete data.
    
    z = lle(X, opt.lle_nn, d+1);
    z = z(1:d,:);
    Z=z;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Setting up graphs of the parameters of the model
    if opt.plot
       [figh, plotstruct] = plotCFAdemo([], [], [], X, EigVal, Mu, U, c, 1);
    end
    
    % Fit the parameters to LLE latent variable coordinates
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Perform the EM algorithm for the optimization in main loop
    tolerance = 5; min_iter=5;
    iter = 1; InitialPhase = 0;
    model.Q = Q; model.SigmaN = SigmaN; model.SigmaC = SigmaC; model.Psi = Psi; model.Pi = Pi; model.Lambda = Lambda;
    model.Mu = Mu; model.Kappa = Kappa; model.Z = Z; model.logdetSigmaN = logdetSigmaN; model.C = C; model.Proj = Proj; model.R = R;
    Fs = zeros(c,1);
    Ind= true(1,c); splitInd=0; c_iter=1;
    while iter <= opt.max_iterations
        
        if iter<opt.lle_iterations
            if ~InitialPhase; disp('*** Start Initial Phase: use Z from LLE ***'); InitialPhase = 1; end
            Z=z;
            [model]=cfa_e_step(X,model,1);
            Fs(Ind,iter) = model.Fs';
        else
            if InitialPhase; disp('*** End Initial Phase ***'); InitialPhase = 0; end
            [model]=cfa_e_step(X,model,1);
            model;
            Fs(:,iter) = model.Fs';
            
        end
        Fh(iter) = sum(Fs(:,iter));
        
        [model]=cfa_m_step(X, model);
        Lambda = model.Lambda; Psi = model.Psi; SigmaC = model.SigmaC; Pi = model.Pi; Kappa = model.Kappa;
        Mu = model.Mu; Q = model.Q; m = model.m; Z = model.Z; C = model.C; Proj = model.Proj; R = model.R;

        %%% Calculate reconstruction error
        x_z = ReconstructX(Z,Q,Lambda,SigmaC,Mu,Kappa);
        
        if opt.InitPCA
            CFA_error = sqrt(sum(sum(((UU1(:,1:dPCA)*x_z)' - X_o).^2)));
        else
            CFA_error = sqrt(sum(sum((x_z - X).^2)));
        end
        disp(['Iteration: ' num2str(iter) ', clusters: ' num2str(c) '| CFA error: ' num2str(CFA_error) ', ratio: ' num2str(CFA_error/PCA_error) ', ~LogLik: ' num2str(sum(Fs(:,iter)))])
        
        
        for k=1:c
            R(:,:,k) = chol(SigmaC(:,:,k));            
            Proj(:,:,k)= Lambda(:,:,k)*R(:,:,k)';
            C(:,:,k)=Lambda(:,:,k)*SigmaC(:,:,k)*Lambda(:,:,k)' + Psi(:,:,k);
%             invC(:,:,k) = invPsi(:,:,k) - invPsi(:,:,k)*Lambda(:,:,k)*invVc(:,:,k)*Lambda(:,:,k)'*invPsi(:,:,k);
%             invC2(:,:,k) = invPsi(:,:,k)*(Psi(:,:,k) - Lambda(:,:,k)*invVc(:,:,k)*Lambda(:,:,k)')*invPsi(:,:,k);
            [U(:,:,k),Sc,~]=svd(C(:,:,k));
            EigK(:,k) = diag(Sc);
        end
        model.Proj = Proj; model.C = C; model.R = R; model.U = U; model.EigK = EigK; 
               
        
        c_iter=c_iter+1;
        % Update number of iterations
        iter = iter + 1;
        if rem(iter, 1) == 0
            if opt.plot
                [figh, plotstruct] = plotCFAdemo(figh, plotstruct, model, X, EigVal, Mu, U, c, 0);
            end
        end        
    end
    
    if opt.plot
        title(plotstruct.sp1,'Finished')
    end
    
    % Transpose to get lowdimensional data representation
    g = Z';
    
    % return the parameters of the mapping
    mapping.Lambda = Lambda;
    mapping.Z = Z;
    mapping.Kappa = Kappa;
    mapping.Mu = Mu;
    mapping.Psi = Psi;
    mapping.Pi = Pi;
    mapping.SigmaC = SigmaC;
    mapping.SigmaN = SigmaN;
    if opt.InitPCA; mapping.InitialPCA = UU1(:,1:dPCA); else mapping.InitialPCA=0;end
    mapping.Q = Q;
    for ii=1:c
        mapping.Sigma(:,:,ii) = Lambda(:,:,ii)*Lambda(:,:,ii)' + Psi(:,:,ii);
    end
    mapping.data_mu = X_mu;
    mapping.Fs=Fs;
%     figure;plot(Fh)
end

function [model]=cfa_e_step(X,model,initialZ)

Lambda = model.Lambda; Psi = model.Psi; Pi = model.Pi; SigmaC = model.SigmaC; SigmaN = model.SigmaN;
Kappa = model.Kappa; Mu = model.Mu; Q = model.Q; logdetSigmaN = model.logdetSigmaN; Z = model.Z;

[D,d,c] = size(Lambda);
n = size(X,2);


% E-step
% ====================================================

% Do some precomputations for speed
invPsi = zeros(size(Psi));
logdetSigmaC = zeros(1, c);
logdetPsi = zeros(1, c);
logPi = zeros(1, c);

for j=1:c
    invPsi(:,:,j) = inv(Psi(:,:,j) + 0.000000000000001*eye(D));
    logdetSigmaC(j) = log(det(SigmaC(:,:,j))+realmin);
    logdetPsi(j) = log(prod(diag(Psi(:,:,j)))+realmin);
    logPi(j) = log(Pi(j)+realmin);
end

% Compute matrices Epsilon, Vc, and m
Eps = zeros(n, c);
S = zeros(n,c);
Vc = zeros(d, d, c); invVc = zeros(d, d, c);
m = zeros(d, n, c);
log2pi = log(2*pi);
const1 = ((D + d)/2)*log2pi;
const2 = (d/2)*log2pi;
I = eye(size(SigmaN(:,:,1)));
for j=1:c
    
    % Precomputations
    Xnc = bsxfun(@minus,X,Mu(:,j));
    Znc = bsxfun(@minus,Z,Kappa(:,j));
    tmpLiPL = Lambda(:,:,j)' * invPsi(:,:,j) * Lambda(:,:,j);
    Vc(:,:,j) = I/(SigmaC(:,:,j) + realmin*eye(d)) + tmpLiPL;
    invVc(:,:,j) = I/(Vc(:,:,j) + realmin*eye(d));
    % Compute Epsilon
    for i=1:n
        normX = (Xnc(:,i) - Lambda(:,:,j)*Znc(:,i));
        S(i,j) = .5*logdetSigmaN(i) - log(Q(i,j)+realmin) + const2;
        Eps(i, j) = (.5 * sum(sum(SigmaC(:,:,j).*(SigmaN(:,:,i) + Znc(:,i) * Znc(:,i)')))) ...
            + (.5 * sum(sum(SigmaN(:,:,i).*tmpLiPL))) ...
            + (.5 * normX' * invPsi(:,:,j) * normX);
        m(:,i,j) = Kappa(:,j) + ((Vc(:,:,j) + realmin*eye(d))\Lambda(:,:,j)')*invPsi(:,:,j)*(Xnc(:,i));
%         m(:,i,j) = Kappa(:,j) + invVc(:,:,j)*Lambda(:,:,j)'*invPsi(:,:,j)*(Xnc(:,i));
    end
    Eps(:,j) = Eps(:,j) - logPi(j) + const1 + 0.5*(logdetSigmaC(j) + logdetPsi(j));
end
Fs = (sum(Q.*(S-Eps)));
F = sum(sum(Q.*(S-Eps)));

% Update estimate of Q. 
% Here we need to use logsumexp trick.
Eps=-Eps;
shiftEps = bsxfun(@minus,Eps,max(Eps,[],2));
expEps = exp(shiftEps);
Q = bsxfun(@rdivide,expEps,(sum(expEps,2)));





I = eye(size(SigmaN(:,:,1)));
if initialZ
    SigmaN = repmat(0.000000000001*eye(d),[1 1 n]);
else
    % Update estimate of SigmaN
    for i=1:n
        tmp = zeros(d, d);
        for j=1:c
            tmp = tmp + Q(i, j) * Vc(:,:,j);
        end
        % Compute covariance matrix
        SigmaN(:,:,i) = I/(tmp + realmin*eye(d));  % code above gave us inv(SigmaN)
        logdetSigmaN(i) = log(det(SigmaN(:,:,i))+realmin);
    end
    
    % Update estimate of Z
    for i=1:n
        tmp = zeros(d, 1);
        for j=1:c
            tmp = tmp + (Q(i, j) * Vc(:,:,j) * m(:,i,j));
        end
        Z(:,i) = SigmaN(:,:,i) * tmp;
    end
%     Zsmooth = smooth(Z',10,'loess'); Z=reshape(Zsmooth,n,d)';

end
model.Q = Q; model.SigmaN = SigmaN; model.Z = Z; model.m = m; model.Fs = Fs;
model.logdetSigmaN = logdetSigmaN; model.Fs = Fs; model.invVc = invVc; model.invPsi = invPsi;

end


function [model]=cfa_m_step(X, model)

Q = model.Q; Z = model.Z; SigmaN = model.SigmaN;
c = size(Q,2);
d = size(Z,1);
[D,n] = size(X);

% M-step
% ====================================================
min_var = 1e-12;                     % minimum STD of Gaussians
% Update estimate of Pi
Pi = sum(Q, 1) ./ n;
for j=1:c
    tmp = 0;
    tmpQ = Q(:,j) ./ (sum(Q(:,j))+realmin);
    
    % Update estimate of Mu and Kappa
    Mu(:,j) = sum(repmat(tmpQ', [D 1]) .* X, 2);
    Kappa(:,j) = sum(repmat(tmpQ', [d 1]) .* Z, 2);
    
    Xnc = bsxfun(@minus,X,Mu(:,j));
    Znc = bsxfun(@minus,Z,Kappa(:,j));
    
    % Update estimate of SigmaC
    for i=1:n
        tmp = tmp + (tmpQ(i) * (SigmaN(:,:,i) + Znc(:,i) * Znc(:,i)'));
    end
    % Enforce some variance
    tmp(1:size(tmp, 1) + 1:end) = max(min_var, tmp(1:size(tmp, 1) + 1:end));
    SigmaC(:,:,j) = tmp;

    % Update estimate of Lambda
    Sc = zeros(D, d);
    for i=1:n
        Sc = Sc + tmpQ(i) * (Xnc(:,i) * Znc(:,i)');
    end
    Lambda(:,:,j) = Sc/(SigmaC(:,:,j) + realmin*eye(d));

    % Update estimate of Psi
    tmpPsi = zeros(D, D); 
%     LL(:,:,j) = Lambda(:,:,j)*Lambda(:,:,j)';
    for i=1:n
        tmpLSL = Lambda(:,:,j) * SigmaN(:,:,i) * Lambda(:,:,j)';
        tmpPsi(1:size(tmpPsi, 2) + 1:end) = tmpPsi(1:size(tmpPsi, 2) + 1:end) + ...
            tmpQ(i) * (((Xnc(:,i) - Lambda(:,:,j) * Znc(:,i)) .^ 2)' + tmpLSL(1:size(tmpLSL, 2) + 1:end));
    end
    Psi(:,:,j) = tmpPsi;
end

model.SigmaC = SigmaC; model.Psi = Psi; model.Pi = Pi; model.Lambda = Lambda;
model.Mu = Mu; model.Kappa = Kappa; 



end
