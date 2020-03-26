function X = ReconstructX(Z,Q,Lambda,SigmaC,Mu,Kappa)

% Function to reconstruct vector in data space, X, given vector in latent
% space, Z.

% Written by Matthew Field

[D,c]=size(Mu);
n=size(Q,1); n2=size(Z,2);
p_c = (1/n)*sum(Q,1);
p_x_c=zeros(n2,c); p_z_c=zeros(n2,c); mm=zeros(D,n2,c);

for k=1:c
%     SigmaL = Lambda(:,:,k)*SigmaC(:,:,k)*Lambda(:,:,k)' + Psi(:,:,k);
%     p_x_c(:,k) = exp(gaussLogprob(Mu(:,k),SigmaL,X'));
    
    p_z_c(:,k) = exp(gaussLogprob(Kappa(:,k),SigmaC(:,:,k),Z'));
    [R(:,:,k), p] = chol(SigmaC(:,:,k));
    Proj(:,:,k) = Lambda(:,:,k);
    %  if any(Proj(:,:,k))
    %  Proj(:,:,k) = bsxfun(@rdivide,Proj(:,:,k),sqrt(sum(Proj(:,:,k).^2)));
    %  end
    
    mm(:,:,k)=bsxfun(@plus,Proj(:,:,k)*bsxfun(@minus,Z,Kappa(:,k)),Mu(:,k));
    % p_z_xc(:,k) = gaussLogprob(zeros(d,1),invVc(:,:,k),(Z-m(:,:,k))');
    % p_x_zc(:,k) = gaussLogprob(zeros(D,1),Psi(:,:,k),(X-mm)');

end
% p_x = p_x_c*p_c(:);
% p_c_x = bsxfun(@rdivide,bsxfun(@times,p_x_c,p_c),p_x);
% z_x = sum(bsxfun(@times,m,permute(p_c_x,[3 1 2])),3);


p_z = p_z_c*p_c(:);
p_c_z = bsxfun(@rdivide,bsxfun(@times,p_z_c,p_c),p_z);
X = sum(bsxfun(@times,mm,permute(p_c_z,[3 1 2])),3);

end

