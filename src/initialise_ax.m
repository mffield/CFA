function info = initialise_ax(info, ax)

% This function is designed to initialise the components of plot for displaying human motion capture data

% Assumes there is a file in same directory 'StructMVN3.mat' which contains
% all of the necessary information to describe the shapes of the 3d model.
% It is based on the Moven motion capture suite (version 3).

% Written by Matthew Field

info.ax = ax;
num_parts=23;
num_COM=14;

info.pg=zeros(1,num_parts);
info.world_frame = hgtransform('Parent',info.ax);
for i=1:num_parts
    info.pg(i) = hgtransform('Parent',info.world_frame);
end

hold(info.ax, 'on')
info.mid_points=plot3(zeros(num_parts,1),zeros(num_parts,1),zeros(num_parts,1),'.',...
    'Color',[1 0.5 0],'MarkerSize', 14, 'Visible',info.vis{3}, 'Parent', info.world_frame);
info.joint_lights=plot3(zeros(num_parts,1),zeros(num_parts,1),zeros(num_parts,1),'.b',...
    'MarkerSize', 14, 'Visible',info.vis{2},'Parent', info.world_frame);
info.COM_lights=plot3(zeros(num_COM,1),zeros(num_COM,1),zeros(num_COM,1),'.g',...
    'MarkerSize', 10, 'Visible',info.vis{4},'Parent', info.world_frame);
info.tCOM = plot3(0,0,0, '.g', 'MarkerSize', 24, 'Visible',info.vis{4}, 'Parent', info.world_frame);
for x=1:num_COM
    info.COM_lines(x) = line(0,0,0,'Color', 'g', 'Marker','.','LineWidth',1,...
        'Visible',info.vis{4},'Parent', info.world_frame);
end

for jj=1:26
    info.jointLines(jj) = line(0,0,0,'Color', [0.5 0.5 0.5],'LineStyle','--','LineWidth',1,'Parent',info.world_frame);
end
info.Line_connect=[2:27; [1:7 5 9:12 5 14:17 1 19:22 1 24:26]];

info.foot_points=plot3(zeros(8,1),zeros(8,1),zeros(8,1),'.m',...
    'MarkerSize', 14, 'Visible',info.vis{4}, 'Parent', info.world_frame);
info.foot_convhull = trisurf(1,0,0,0,'Parent',info.ax,'EdgeColor','none',...
    'Visible',info.vis{4});
info.foot_shadow = plot(info.ax,0,0,'-m','LineWidth',1,'Visible',info.vis{4});

info.mapping_arrows = quiver3(0,0,0,0,0,0,'AutoScale','off','Color','k','LineStyle','-','LineWidth',3,'Parent',info.world_frame);


set(gcf,'Renderer','opengl');
light('Position', [-0.5 0.5 1],'Parent', info.ax)
lighting flat

set(info.ax, 'XLim', [-0.8 0.8], 'YLim', [-0.8 0.8], 'ZLim', [-1 1], 'View', [45, 30]);
set(get(info.ax,'XLabel'),'String','x_1','Interpreter','tex')
set(get(info.ax,'YLabel'),'String','x_2','Interpreter','tex')
set(get(info.ax,'ZLabel'),'String','x_3','Interpreter','tex')
info.num_skel=1;
hold(info.ax, 'on')

p = mfilename('fullpath');
pc_path = p(1:find(p=='\',1,'last'));
mac_path = p(1:find(p=='/',1,'last'));

if isempty(pc_path);  load([mac_path '/StructMVN3.mat']); 
elseif isempty(mac_path);    load([pc_path '/StructMVN3.mat']);
else disp('Error: Failed to interpret path of current directory.')
end
% 'mass distribution' + shapes + tree

info.COM_male = COM_male; info.COM_female = COM_female;
Def_skel = Struct;
info.Def_segment = Def_skel.segment;
info.segment = Def_skel.segment;
segment = info.segment;


colors = [repmat([0.33 0.33 0.34],7,1);repmat([0.8 0.1 0.1],4,1);...
    repmat([0.1 0.8 0.1],4,1);repmat([0.8 0.1 0.1],4,1);repmat([0.1 0.8 0.1],4,1)];
for ii=1:num_parts
    info.segment(ii).h = patch('Vertices',segment(ii).vert,'Faces',segment(ii).faces,...
        'Visible',info.vis{1},'FaceColor',colors(ii,:),...
        'EdgeColor','none','Parent',info.pg(ii));
    alpha(info.segment(ii).h,0.6);
end
info.Def_segment = info.segment;
     
    

    
end