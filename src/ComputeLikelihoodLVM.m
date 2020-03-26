function Likelihood=ComputeLikelihoodLVM(Z, Q, SigmaC, Kappa)

% Function to compute the log likelihood of a point in the latent variable
% space. Input arguments include the model parameters 
% the set of latent variable points for which to compute the likelihood.

% Written by Matthew Field

[N, c] = size(Q);
p_c = (1/N)*sum(Q,1);
p_z_c=zeros(size(Z,1),c);
for k=1:c
    p_z_c(:,k) = exp(gaussLogprob(Kappa(k,:)',SigmaC(:,:,k),Z));
end
p_z = p_z_c*p_c(:);
Likelihood = p_z;

end