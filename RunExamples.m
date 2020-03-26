%   CFA ALGORITHM
%   v1.0 created from 30/1/2014

% This script is designed to demonstrate the globally coordinated mixture of
% factor analyzers (GCMFA) model or otherwise known as Coordinated Factor Analysis (CFA)
% for dimensionality reduction. It is tested on synthetic data of a 2D plane
% curved in 3D and on human motion capture data. The performance is compared to using a PCA projection.
% An interactive plot is shown to demonstrate mapping from data space to
% latent space and vice versa. Please refer to the article for more details.
% 
% Written by Matthew Field
%
% This algorithm is a reduced version of methods published in:
% *M. Field, D. Stirling, Z. Pan and F. Naghdy, "Learning Trajectories for Robot
% Programing by Demonstration Using a Coordinated Mixture of Factor Analyzers,"
% in IEEE Transactions on Cybernetics, vol. 46, no. 3, pp. 706-717, March 2016.
% https://doi.org/10.1109/TCYB.2015.2414277

%clc; 
clear;
addpath(genpath('.\src'))

% Pick one of the examples
Example = [0    % Simple curve data (1/2 of S)
           1    % S-curve data
           0];  % Motion capture: 69 Euler angles

   randn('seed',3)  
   rand('seed',3) 
if Example(1) 
    
    load data/SwissRoll

    CurveData = X(labels<8,:);
    CurveData_label=labels(labels<8);
     
    opt.K = 4; % set the number of components.
    opt.d = 2; % set the dimensionality of latent space.
    opt.max_iterations = 40; % maximum number of iterations for algorithm
    opt.lle_iterations = 10; opt.lle_nn = 20; % set parameters for local linear embedding
    opt.plot=1; % enable plotting of the parameters.
    opt.labels=CurveData_label;

    [mappedX, mapping]=CFA(CurveData,opt);
    
    % Interact with the mapping function
    mapping.datatype = '3Dcurve'; mapping.type = 'CFA'; mapping.plotPosition=[100 20];
    MappingInteractionFigure(CurveData,mappedX,mapping,CurveData_label,[])
    
elseif Example(2) % Use 2D plane curved into shape of an S.
        
    GenerateSdata; 
    SData = X(:,1:2:end)'; S_label = zeros(1,size(SData,1));
    S_label = round(10*SData(:,1));
    
    opt.K = 12;
    opt.d = 2;
    opt.max_iterations = 40;
    opt.lle_iterations = 10; opt.lle_nn = 10; % set parameters for LLE
    opt.verbose=1; opt.plot=1; opt.labels=S_label; opt.InitPCA = 0;
    tic
    [mappedX, mapping]=CFA(SData,opt);
    toc

    mapping.datatype = '3Dcurve'; mapping.type = 'CFA'; mapping.plotPosition=[100 20];
    MappingInteractionFigure(SData,mappedX,mapping,S_label,[])
        
elseif Example(3)
    
    load data/WalkingData
    D = Data.EulerAngle*180/pi; % convert to degrees.
    D(:,[1 55:57 67:69])=0;     % toe angle and body heading direction set to zero.
    X = D(1:4:end,~ismember(1:size(D,2),[1 55 56 57 67 68 69])); % downsample the data. %%% there was 21
    X_mu = mean(X);             % store the mean.
    [U,S,V]=svd(cov(X));        % compute the eigenvalues/vectors of data covariance (PCA).
    
    Y = (U(:,1)'*(bsxfun(@minus,X,X_mu))')';    %project data along first principal component.
    Ylabel = ceil(Y(:,1)-min(Y(:,1))+1);        %create a color label vector.

    % CFA
    opt.K = 10; opt.d = 2;      % set the number of components of the model and the dimension of the latent space.
    opt.max_iterations = 30;    % set the maximum number of iterations of the algorithm.
    % set parameters for LLE initialization. lle_iterations defines for how many iterations the 
    % algorithm fit the parameters to the latent coordinates Z obtained from LLE.
    % lle_nn defines how many nearest neighbours are used to create
    % neighbourhood graph for LLE (a necessary parameter for LLE, dataset dependent).
    opt.lle_iterations = 10; opt.lle_nn = 25; 
    opt.plot=0; % set the plot flag to zero (it is only for visualizing the parameters in 3D->2D case) 
    
    % InitPCA controls proportion of data variance to maintain in initial
    % PCA projection. This is optional but often improves performance in reconstruction error and computation time.
    % This is because it reduces the number of parameters to fit by eliminating factors below a noise threshold.
    opt.InitPCA = 0.98; % we set to retain 98% of variance but we could modify it to threshold on the minimum eigenvalue - derived from noise floor. 
    [mappedX, mapping]=CFA(X,opt); % run the algorithm
    
    % set identifiers/parameters for graphing
    mapping.type = 'CFA'; mapping.datatype = 'Mocap'; mapping.plotPosition=[100 20];
    
    % Simple PCA to compare the reconstruction with CFA
    PCAmap.Z = X*U(:,1:2); % project along 2 principal components 
    PCAmap.InitialPCA = U(:,1:2); % store the eigenvectors for data reconstruction
    PCAmap.data_mu = X_mu; PCAmap.type = 'PCA'; PCAmap.datatype = 'Mocap'; 
    PCAmap.plotPosition = [100 550];
    
    info.skel = skeltree;
    MappingInteractionFigure(X,PCAmap.Z,PCAmap,Ylabel,info)
    MappingInteractionFigure(X,mappedX,mapping,Ylabel,info)

end
        

