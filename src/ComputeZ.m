function Z=ComputeZ(X,Q,Lambda,SigmaC,Mu,Psi,Kappa,invVc,invPsi)
% Function to compute latent variable, z, from data space and model
% parameters

% Written by Matthew Field

[n,c]=size(Q);
p_c = (1/n)*sum(Q,1);
p_x_c=zeros(1,c); n2=size(X,2);
for k=1:c
    SigmaL = Lambda(:,:,k)*SigmaC(:,:,k)*Lambda(:,:,k)' + Psi(:,:,k);
    p_x_c(:,k) = exp(gaussLogprob(Mu(:,k),SigmaL,X'));
    
    m(:,:,k) = Kappa(:,k) + invVc(:,:,k)*Lambda(:,:,k)'*invPsi(:,:,k)*(X-Mu(:,k));
end
p_x = p_x_c*p_c(:);
p_c_x = bsxfun(@rdivide,bsxfun(@times,p_x_c,p_c),p_x);
Z = sum(bsxfun(@times,m,permute(p_c_x,[3 1 2])),3);

end

