function logp = gaussLogprob(arg1, arg2, arg3, arg4, arg5, arg6)
% Log pdf of multivariate Gaussian
% L  = gaussLogprob(model, X)
% L = gaussLogprob(mu, Sigma, X)
% X(i,:) is i'th case, can contain NaNs for missing values.
% Sigma can be a vector; this is interpreted as a diagonal matrix.
% *** In the univariate case, Sigma is the variance, not the standard
% deviation! ***
%
% Examples
% L  = gaussLogprob(zeros(3,1), randpd(3), rand(10,3))
% L  = gaussLogprob(zeros(3,1), diag(randpd(3)), rand(10,3))
% L = gaussLogprob(structure(mu, Sigma), X)

% This file is from pmtk3.googlecode.com


switch nargin
    
    case 6,  mu = arg1; Sigma = arg2; W = arg3; sigma2 = arg4; invM = arg5; X = arg6;
    case 3,  mu = arg1; Sigma = arg2; X = arg3;
    case 2, model = arg1; mu = model.mu; Sigma = model.Sigma; X = arg2;
    otherwise
        error('bad num args')
end

mu = mu(:);
if any(isnan(X(:)))
    logp = gaussLogprobMissingData(structure(mu, Sigma), X);
    return;
end
if isscalar(mu)
    X = X(:);
end
[N, d] = size(X);

if d == 1
    X = X(:) - mu(:); % each data case evaluated under a different mu 
else
    X = bsxfun(@minus, X, mu(:)');
    
end
if isvector(Sigma) && (numel(Sigma) > 1) % diagonal case
  sig2 = repmat(Sigma', N, 1);
  tmp  = -(X.^2)./(2*sig2) - 0.5*log(2*pi*sig2);
  logp = sum(tmp, 2);
else
  % Full covariance case
  %   slightly faster version  
  %%%
  if nargin==6      %%% PPCA version
      d=size(W,1);

      invSigma = (1/sigma2)*(eye(d)-W*invM*W');
      logp = -0.5*sum((X*invSigma).*X, 2); % -0.5*X*invSigma*X';
      logZ = (d/2)*log(2*pi) + 0.5*logdet(Sigma);
      logp = logp - logZ;
      %%%
  else
      [R, p]   = chol(Sigma);
      if p~=0
          logp=0;
      else
          logp = -0.5*sum((X/R).^2, 2);
      end
      logZ = 0.5*d*log(2*pi) + sum(log(diag(R)));
      logp = logp - logZ;
  end
  %assert(approxeq(logp, logp2))
end
end

