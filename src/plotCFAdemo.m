function [figh, plotstruct] = plotCFAdemo(figh, plotstruct, model, X, EigVal, Mu, U, c, init)

std_scale=2;  minx = 1.2*min(X,[],2); maxx = 1.2*max(X,[],2);
colors=[0 0 1];

if init
    if isempty(plotstruct)
        figh = figure; set(figh,'Color','w');  plotstruct.sp1 = axes;% subplot(121);
    else
        cla(plotstruct.sp1)
    end
    plot3(X(1,:),X(2,:),X(3,:),'.k','MarkerSize',10); hold on;  grid on;
    
    for k=1:c
        % The covariance of the Gaussian is: C = Lambda * Lambda' + Psi;
        EigK = std_scale*sqrt(EigVal(:,k));
        plotstruct.cluster_ellip(k) = hgtransform('Parent',plotstruct.sp1);
        [xx,yy,zz]=ellipsoid(0,0,0,EigK(1),EigK(2),EigK(3),10);
        plotstruct.cluster_surf(k) = surf(xx,yy,zz,'EdgeColor','none','Parent',plotstruct.cluster_ellip(k),'FaceColor',colors);
        alpha(plotstruct.cluster_surf(k),0.2);
        set(plotstruct.cluster_ellip(k),'Matrix',[U(1:3,1:3,k) [Mu(1,k);Mu(2,k);Mu(3,k)];[0 0 0 1]])
        plotstruct.bases(k,1) = quiver3(Mu(1,k),Mu(2,k),Mu(3,k),EigK(1)*U(1,1,k),EigK(1)*U(2,1,k),EigK(1)*U(3,1,k),'LineWidth',3,'MaxHeadSize',3,'Color','g','AutoScale','off');
        plotstruct.bases(k,2) = quiver3(Mu(1,k),Mu(2,k),Mu(3,k),EigK(2)*U(1,2,k),EigK(2)*U(2,2,k),EigK(2)*U(3,2,k),'LineWidth',3,'MaxHeadSize',3,'Color','r','AutoScale','off');
    end
    set(plotstruct.sp1,'XLim',[minx(1) maxx(1)],'YLim',[minx(2) maxx(2)],'ZLim',[minx(3) maxx(3)],'View',[25 8])
    drawnow;
else 
    if 1
        cla(plotstruct.sp1);
        plot3(X(1,:),X(2,:),X(3,:),'.k','MarkerSize',10,'Parent',plotstruct.sp1); hold on;
        %                     scatter3(X(1,:),X(2,:),X(3,:),20*ones(length(Z),1),opt.labels,'filled','Parent',sp1);hold on;
        for k=1:c
            % The covariance of the Gaussian is: C = Lambda * Lambda' + Psi;
            plotstruct.cluster_ellip(k) = hgtransform('Parent',plotstruct.sp1);
            [xx,yy,zz]=ellipsoid(0,0,0,std_scale*sqrt(model.EigK(1,k)),std_scale*sqrt(model.EigK(2,k)),std_scale*sqrt(model.EigK(3,k)),10);
            plotstruct.cluster_surf(k) = surf(xx,yy,zz,'EdgeColor','none','Parent',plotstruct.cluster_ellip(k),'FaceColor',colors);
            alpha(plotstruct.cluster_surf(k),0.2);
            set(plotstruct.cluster_ellip(k),'Matrix',[model.U(1:3,1:3,k) [model.Mu(1,k);model.Mu(2,k);model.Mu(3,k)];[0 0 0 1]])
            plotstruct.bases(k,1) = quiver3(model.Mu(1,k),model.Mu(2,k),model.Mu(3,k),std_scale*model.Proj(1,1,k),std_scale*model.Proj(2,1,k),std_scale*model.Proj(3,1,k),'LineWidth',3,'MaxHeadSize',3,'Color','g','AutoScale','off','Parent',plotstruct.sp1);
            plotstruct.bases(k,2) = quiver3(model.Mu(1,k),model.Mu(2,k),model.Mu(3,k),std_scale*model.Proj(1,2,k),std_scale*model.Proj(2,2,k),std_scale*model.Proj(3,2,k),'LineWidth',3,'MaxHeadSize',3,'Color','r','AutoScale','off','Parent',plotstruct.sp1);
        end
        
    end
    for k=1:c

        [xx,yy,zz]=ellipsoid(0,0,0,std_scale*sqrt(model.EigK(1,k)),std_scale*sqrt(model.EigK(2,k)),std_scale*sqrt(model.EigK(3,k)),10);
        set(plotstruct.cluster_surf(k),'XData',xx,'YData',yy,'ZData',zz);
        set(plotstruct.cluster_ellip(k),'Matrix',[model.U(1:3,1:3,k) [model.Mu(1,k);model.Mu(2,k);model.Mu(3,k)];[0 0 0 1]])
        set(plotstruct.bases(k,1),'XData',model.Mu(1,k),'YData',model.Mu(2,k),'ZData',model.Mu(3,k),...
            'UData',std_scale*model.Proj(1,1,k),'VData',std_scale*model.Proj(2,1,k),'WData',std_scale*model.Proj(3,1,k),'AutoScale','off');
        set(plotstruct.bases(k,2),'XData',model.Mu(1,k),'YData',model.Mu(2,k),'ZData',model.Mu(3,k),...
            'UData',std_scale*model.Proj(1,2,k),'VData',std_scale*model.Proj(2,2,k),'WData',std_scale*model.Proj(3,2,k),'AutoScale','off');
    end
    set(plotstruct.sp1,'XLim',[minx(1) maxx(1)],'YLim',[minx(2) maxx(2)],'ZLim',[minx(3) maxx(3)])
    drawnow;
end

