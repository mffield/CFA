function MappingInteractionFigure(X,mappedX,mapping,labels,info)

% This function implements a real-time interactive figure for visualising
% the model. 

% Written by Matthew Field

fig_interact = figure('Color','w','Position',[mapping.plotPosition 450 450]); ax1 = axes; hold on;
handles.h2D = plot3(0,0,0,'-g','LineWidth',2,'Parent',ax1);
figure('Color','w','Position',[[mapping.plotPosition+[500 0]] 450 450]); ax2 = axes; hold on;

title(ax1,[mapping.type ': Latent space  Z']);
title(ax2,[mapping.type ': Data space  X']);

if strcmp(mapping.type,'CFA')
    % c=size(mapping.Kappa,2);
    % colors=repmat([0 0 1],c,1);
    % std_scale=2;  minx = 1.2*min(mappedX,[],2); maxx = 1.2*max(mappedX,[],2);
    % for k=1:c
    %     % The covariance of the Gaussian is: C = Lambda * Lambda' + Psi;
    %     [U(:,:,k),EigV,~]=svd(mapping.SigmaC(:,:,k));
    %     EigVal(:,k) = diag(EigV);
    %     EigK = std_scale*sqrt(EigVal(:,k));
    %     cluster_ellip(k) = hgtransform('Parent',ax1);
    %     [xx,yy,zz]=ellipsoid(0,0,0,EigK(1),EigK(2),0,10);
    %     cluster_surf(k) = surf(xx,yy,zz,'EdgeColor','none','Parent',cluster_ellip(k),'FaceColor',colors(k,:));
    %     alpha(cluster_surf(k),0.2);
    %     set(cluster_ellip(k),'Matrix',[[U(1:2,1:2,k) [0;0];[0 0 0]] [mapping.Kappa(1,k);mapping.Kappa(2,k);0];[0 0 0 1]])
    %     bases(k,1) = quiver3(mapping.Kappa(1,k),mapping.Kappa(2,k),0,EigK(1)*U(1,1,k),EigK(1)*U(2,1,k),0,'LineWidth',2,'Color','g','AutoScale','off','Parent',ax1);
    %     bases(k,2) = quiver3(mapping.Kappa(1,k),mapping.Kappa(2,k),0,EigK(2)*U(1,2,k),EigK(2)*U(2,2,k),0,'LineWidth',2,'Color','r','AutoScale','off','Parent',ax1);
    % end
    % set(sp1,'XLim',[minx(1) maxx(1)],'YLim',[minx(2) maxx(2)],'ZLim',[minx(3) maxx(3)],'View',[25 8])
    % drawnow;
    
    
    No_Meshpoints=150;
    minz = 1.5*min(mapping.Z,[],2);
    maxz = 1.5*max(mapping.Z,[],2);
    [Mx, My]=meshgrid(linspace(minz(1),maxz(1),No_Meshpoints),linspace(minz(2),maxz(2),No_Meshpoints));
    Z = [reshape(Mx,No_Meshpoints^2,1),reshape(My,No_Meshpoints^2,1)];
    
    Lik=ComputeLikelihoodLVM(Z,mapping.Q,mapping.SigmaC,mapping.Kappa');
    Lik(Lik>0.2)=0.2; %0.4
    L = reshape(Lik,No_Meshpoints,No_Meshpoints);
    contourf(Mx,My,L,'EdgeColor','none','Parent',ax1); colormap(ax1,gray)
    set(ax1,'XLim',[minz(1) maxz(1)],'YLim',[minz(2) maxz(2)])
end
if strcmp(mapping.datatype,'Mocap')
    %downsample the data points to plot?
    ds=2; mX = mappedX(1:ds:end,:);
    plot3(mX(:,1),mX(:,2),ones(size(mX,1),1),'rx','Parent',ax1);
    vis = [1 1 1 1 0];
    for nn=1:5
        if  vis(nn); info.vis{nn} = 'on';
        else info.vis{nn} = 'off';
        end
    end
    info = initialise_ax(info, ax2);
    handles.info = info;
end


if strcmp(mapping.datatype,'3Dcurve')
    
    %downsample the data points to plot?
    ds=2; mX = mappedX(1:ds:end,:);
    plot3(mX(:,1),mX(:,2),ones(size(mX,1),1),'xr','Parent',ax1);
    X = X - repmat(mean(X, 1), [size(X, 1) 1]);
    scatter3(X(:,1),X(:,2),X(:,3),20*ones(length(mappedX),1),labels,'filled');
    set(ax2,'View',[30 10])
    handles.h3D = plot3(0,0,0,'-k','LineWidth',2,'Parent',ax2);
end


set(fig_interact,'WindowButtonDownFcn',{@EnableDrawFunction,mapping,handles},'UserData',0);


end


function EnableDrawFunction(h,~,mapping,handles)
    
    if get(h,'UserData')
        disp('Disable Draw')
        set(h,'WindowButtonMotionFcn',[]);
        set(h,'UserData',0)

    else
        disp('Enable Draw')
        set(h,'WindowButtonMotionFcn',{@CallPosition, mapping, handles});
        set(h,'UserData',1)
        set(handles.h2D,'XData',0,'YData',0,'ZData',0)
        if isfield(handles,'h3D');  set(handles.h3D,'XData',0,'YData',0,'ZData',0); end
    end
end


function CallPosition(h,~,mapping,handles)
    
    [z1,z2]=localCheckPointPosition(h,gca);
    if ~isempty(z1)

        Z = [z1;z2];
        Z_ = [get(handles.h2D,'XData'); get(handles.h2D,'YData'); get(handles.h2D,'ZData')];


        if strcmp(mapping.type,'CFA')
            x_pca = ReconstructX(Z(1:2,:),mapping.Q,mapping.Lambda,mapping.SigmaC,mapping.Mu,mapping.Kappa)';
            if mapping.InitialPCA
                x = (x_pca*mapping.InitialPCA');
            else 
                x = x_pca;
            end
        elseif strcmp(mapping.type,'PCA')
            x = (Z(1:2,:)'*mapping.InitialPCA');
        end

        if strcmp(mapping.datatype,'Mocap')
            data.pos = zeros(1,69);
            x = x + mapping.data_mu;
            EulerAngle = ([0 x(1:53) zeros(1,3) x(54:62) zeros(1,3)])*pi/180;
            for jj=1:size(EulerAngle,2)/3
                data.q(4*(jj-1)+1:4*jj) = M2Q(rpy2matrix(EulerAngle([3*(jj-1)+2 3*jj 3*(jj-1)+1])));
            end
            Draw_pose(handles.info, data,1);
            
        elseif strcmp(mapping.datatype,'3Dcurve')
            X_hist = [get(handles.h3D,'XData'); get(handles.h3D,'YData'); get(handles.h3D,'ZData')];
            if size(X_hist,2)==1; X_hist = x'; end
            set(handles.h3D,'XData',[X_hist(1,:) x(1)],'YData',[X_hist(2,:) x(2)],'ZData',[X_hist(3,:) x(3)])
        end

        
        if size(Z_,2)==1; Z_ = Z; end
        if any(Z_(1:2,end)) && any((Z_(1:2,end)'-Z').^2>0.5)
            set(handles.h2D,'XData',[Z_(1,1:end-1) Z(1) Z(1)],'YData',[Z_(2,1:end-1) Z(2) Z(2)],'ZData',ones(1,size(Z_,2)+1))
        else
            set(handles.h2D,'XData',[Z_(1,:) Z(1)],'YData',[Z_(2,:) Z(2)],'ZData',ones(1,size(Z_,2)+1))
        end
        drawnow;
    end
end


function point = localGetNormCursorPoint(figHandle)

point = get(figHandle, 'currentPoint');
figPos = get(figHandle, 'Position');
% Normalise the point of the curstor
point(1) = point(1)/figPos(3);
point(2) = point(2)/figPos(4);

end
function [x, y] = localGetNormAxesPoint(point, axesHandle)

position = get(axesHandle, 'Position');
x = (point(1) - position(1))/position(3);
y = (point(2) - position(2))/position(4);
lim = get(axesHandle, 'XLim');
x = x*(lim(2) - lim(1));
x = x + lim(1);
lim = get(axesHandle, 'YLim');
y = y*(lim(2) - lim(1));
y = y + lim(1);

end

function [x, y] = localCheckPointPosition(fig,ax)

point = localGetNormCursorPoint(fig);
% get the position of the axes
position = get(ax, 'Position');

% Check if the pointer is in the axes
if point(1) > position(1) ...
      & point(1) < position(1) + position(3) ...
      & point(2) > position(2) ...
      & point(2) < position(2) + position(4);
  
  % Rescale the point according to the axes
  [x y] = localGetNormAxesPoint(point, ax);

  % Find the nearest point
else
  x = [];
  y = [];
end
end