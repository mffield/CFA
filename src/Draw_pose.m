function [totalCOM, COMp, MIDp, FOOTPRINT, JOINTS]=Draw_pose(info, D, plotting)

% Function to plot human motion capture data. This can be used to plot
% joint positions, midpoints between joints, centre of mass position, body centre of mass and
% the span of both footprints. Plots may be positions, lines between joints
% or a 3D model. The skeletal model is based upon Moven (version 3) with 23
% body segments. It is initialised with the function 'initialise_ax' which
% provides the 'info' input argument.

%Input: info - generated with function 'initialise_ax', 
%       D is a vector of Euler angles of length 69 (which is 3x23 body segments)
%       plotting is a flag to indicate whether to update the plot or only calculate joint positions (kinematics mode only).
% Output: totalCOM is a position vector for the overall centre of mass of the body
%       COMp is a list of position vectors of the centre of mass for each body segment
%       MIDp is a list of position vectors for the mid point of every body segment
%       FOOTPRINT is a list of position vectors for four locations on each foot
%       JOINTS is a list of position vectors for every joint in the body.

% Written by Matthew Field

format compact
if plotting
    vis = get(info.segment(1).h, 'Visible');
    set(info.ax, 'XLim', [-1 1])
    set(info.ax, 'YLim', [-1 1])
    skel = info.skel;    
else
    skel = info;
end    
    P = skel.segJoint;
    
%     set(info.ax, 'XLim', [body_frame(1)-1.5 body_frame(1)+1.5])
%     set(info.ax, 'YLim', [body_frame(2)-1.5 body_frame(2)+1.5])

    COMp=zeros(28,3);
%     D.q(1:4)=[0 0 0 1]; 
    FOOTPRINT=repmat([0;0;0;1]',8,1);
    P(1,:)=[0 0 0];
    for ii=1:28; T(:,:,ii)=[eye(3,3) P(ii,:)';0 0 0 1]; RT(:,:,ii)=eye(4,4); end
  
    Tree1 = [1:7 0 8:11 0 12:15 0 16:19 0 20:23 0];
    Tree2 = [0:7 5 9:12 5 14:17 1 19:22 1 24:27];
    Tree3 = [0 1:5 6 7 0 8:11 0 12:15 0 16:19 0 20:23];
    FootIndex = [zeros(1,21) 1 1 0 0 0 1 1];
    CMI = [0 0 0 skel.clen(1:4) skel.clen(7) 0 skel.clen(8:11) 0 skel.clen(12:15) 0 skel.clen(16:19) 0 skel.clen(20:23)];
    MASS = [0 0 0 skel.mass(1:4) skel.mass(7) 0 skel.mass(8:11) 0 skel.mass(12:15) 0 skel.mass(16:19) 0 skel.mass(20:23)];
    fp=1;
    MIDp = zeros(28,3);
    for ii=1:28
        if Tree1(ii)
            [TT2 TX] = Q2M([0 0 0],D.q(1+4*(Tree1(ii)-1):4*Tree1(ii)));
            T(1:3,1:3,ii) = TX(1:3,1:3);
        end
        if ii==1
            RT(:,:,ii) = T(:,:,ii);
        end
        if ii>1
            RT(:,:,ii) = RT(:,:,Tree2(ii))*T(:,:,ii);
            MIDp(ii,:) = (RT(1:3,4,ii)-RT(1:3,4,Tree2(ii)))/2 + RT(1:3,4,Tree2(ii));
            if CMI(ii)
                COMp(ii,:) = (RT(1:3,4,ii)-RT(1:3,4,Tree2(ii)))*CMI(ii) + RT(1:3,4,Tree2(ii));
            end
        end
        if FootIndex(ii)            
            FOOTPRINT(fp,1:3)=(RT(1:3,1:3,ii)*skel.segment(Tree3(ii)).point(3).pos' + RT(1:3,4,Tree2(ii)))'; 
            fp=fp+1;
            if ii==23 || ii==28
                FOOTPRINT(fp,1:3)=RT(1:3,1:3,ii)*skel.segment(Tree3(ii)).point(2).pos' + RT(1:3,4,Tree2(ii)); fp=fp+1;% toe end
                FOOTPRINT(fp,1:3)=RT(1:3,1:3,ii)*skel.segment(Tree3(ii)).point(4).pos' + RT(1:3,4,Tree2(ii)); fp=fp+1;
            end
        end
        if plotting 
            if strcmp(vis,'on') && (Tree3(ii)~=0) && (Tree2(ii)~=0)
                RoT = [RT(1:3,1:3,Tree2(ii)) MIDp(ii,:)'; 0 0 0 1];
                if ii==2
                    RoT = [RT(1:3,1:3,Tree2(ii)) [0 0 0]'; 0 0 0 1];
                end
                set(info.pg(Tree3(ii)),'Matrix',RoT)
            end
        end
    end
    JOINTS = squeeze(RT(1:3,4,:))';
       
    totalCOM = COMp'*MASS'/skel.total_mass;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    %%%% set the standard plots
    if plotting
        set(info.joint_lights, ...
            'XData',JOINTS(:,1),'YData',JOINTS(:,2),'ZData',JOINTS(:,3));
        set(info.mid_points, ...
            'XData',MIDp(:,1),'YData',MIDp(:,2),'ZData',MIDp(:,3));
        con = info.Line_connect;
        for kk=1:26
            set(info.jointLines(kk),'XData',[JOINTS(con(2,kk),1) JOINTS(con(1,kk),1)],...
                'YData',[JOINTS(con(2,kk),2) JOINTS(con(1,kk),2)],...
                'ZData',[JOINTS(con(2,kk),3) JOINTS(con(1,kk),3)]);
        end
        set(info.COM_lights, 'XData',COMp(:,1),'YData',COMp(:,2),'ZData',COMp(:,3));
        set(info.tCOM, 'XData',totalCOM(1),'YData',totalCOM(2),'ZData',totalCOM(3));
    end
end