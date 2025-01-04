function f = draw3Dstreamlines(StreamLines, Xmax, Xmin, Ymax, Ymin, Zmax, Zmin, ScalarField, MinAlpha, Seed, CLabel, OpenColorBar, AdjustIntegralConstant)
%% Adjusting Scalar Fields （调整标量场）
if min(ScalarField, [], 'all') < 0 && AdjustIntegralConstant
    ScalarField = ScalarField + abs(min(ScalarField, [], 'all'));
end
%% Draw streamlines （绘制流线）
N_x = size(ScalarField, 2); N_y = size(ScalarField, 1); N_z = size(ScalarField, 3);
% MinAlpha = 0.5;
MaxAlpha = 1.0;
f = figure; 
% f = 1;
hold on
axis tight
axis equal
grid off
box on
ax = gca;
ax.BoxStyle = 'full';
N_Streamline = size(StreamLines, 1);
Width = cell(N_Streamline, 1);
XLength = Xmax - Xmin;
YLength = Ymax - Ymin;
ZLength = Zmax - Zmin;
Scaler = 0.1 * mean([XLength, YLength, ZLength]) / 20;
for n = 1:N_Streamline
    Width{n} = Scaler.*ones(size(StreamLines{n}, 1), 1);
end
xlim([Xmin, Xmax]);
ylim([Ymin, Ymax]);
zlim([Zmin, Zmax]);
h = streamtube(StreamLines(:, 1), Width, [1 20]);
MeanCData = zeros(N_Streamline, 1);
MaxD = -realmax;
MinD = realmax;
for n = 1:N_Streamline
    XData = h(n).XData;
    XData(XData < Xmin) = Xmin;
    XData(XData > Xmax) = Xmax;
    YData = h(n).YData;
    YData(YData < Ymin) = Ymin;
    YData(YData > Ymax) = Ymax;
    ZData = h(n).ZData;
    ZData(ZData < Zmin) = Zmin;
    ZData(ZData > Zmax) = Zmax;
    sz = size(XData);
    PData = [XData(:), YData(:), ZData(:)];
    [IndexX, IndexY, IndexZ] = convertAbsoluteCoordinatesToGrid3D(PData, N_x, ...
       Xmax, Xmin, N_y, Ymax, Ymin, N_z, Zmax, Zmin);
    Index = sub2ind([N_y, N_x, N_z], IndexY, IndexX, IndexZ);
    CData = ScalarField(Index);
    h(n).CData = reshape(CData, sz);
    if ~isempty(CData)
        MeanCData(n) = mean(CData, 'all','omitnan');
        % MeanCData(n) = max(CData, [], 'all');
        MaxD = max(MaxD, max(CData, [], 'all'));
        MinD = min(MinD, min(CData, [], 'all'));
    end
end
if MinAlpha < MaxAlpha && MaxAlpha <= 1.0 && MinAlpha >= 0.0
    AlphaData = zeros(N_Streamline, 1);
    CDataRange = range(MeanCData);
    CDataMin = min(MeanCData, [], 'all');
    AlphaRange = MaxAlpha - MinAlpha;
    for n = 1:N_Streamline
        AlphaData(n) = min(AlphaRange / CDataRange * (MeanCData(n) - CDataMin) + MinAlpha, 1.0);
    end
    for n = 1:N_Streamline
        h(n).FaceAlpha = AlphaData(n);
    end
else
    for n = 1:N_Streamline
        h(n).FaceAlpha = MaxAlpha;
    end
end
if Seed
    for n = 1:N_Streamline
        scatter3(StreamLines{n, 2}(1), StreamLines{n, 2}(2), StreamLines{n, 2}(3), 128, 'red', 'filled', 'hexagram');
    end
end

set(gca,'fontsize',16,'FontWeight','bold');%设置坐标轴字体大小
set(gca,'linewidth',0.5,'FontWeight','bold'); %坐标线粗0.5磅

% xlabel('x');
% ylabel('y');
% zlabel('z');
% set(get(gca,'XLabel'),'FontSize',16,'FontWeight','bold');%
% set(get(gca,'YLabel'),'FontSize',16,'FontWeight','bold');
% set(get(gca,'ZLabel'),'FontSize',16,'FontWeight','bold');
% set(get(gca,'TITLE'),'FontSize',16,'FontWeight','bold');

set(gca, 'xtick', [], 'xticklabel', [], 'XLabel', [])
set(gca, 'ytick', [], 'yticklabel', [], 'YLabel', [])
set(gca, 'ztick', [], 'zticklabel', [], 'ZLabel', [])

MaxD = max(ScalarField, [], 'all') + eps;
MinD = min(ScalarField, [], 'all') - eps;
clim([MinD, MaxD]);
Ticks = clim;
Ticks = linspace(Ticks(1), Ticks(2), 6);

TickLabels = cell(size(Ticks));
for i = 1:size(Ticks, 2)
    % TickLabels{i} = num2str(Ticks(i)); 
    TickLabels{i} = erase(sprintf('%3.2g', Ticks(i)), '+'); 
end
load('cool2warm.mat', 'map');
colormap(map);
if(OpenColorBar)
    c = colorbar('Ticks', Ticks);
    c.Label.String = "";
    c.TickLabels = TickLabels;
end
% if isempty(CLabel)
%     c.Label.String = 'User defined scalar field';
% else
%     c.Label.String = CLabel;
% end
view(3)
shading interp
%% Adjust the position of the colorbar （调整colorbar的位置）
if (sqrt(XLength^2 + YLength^2) < ZLength || 1) && OpenColorBar
    ax.Units = "centimeters";
    AxisPosition = tightPosition(ax, IncludeLabels=true);
    camera_position = campos;
    camera_target = camtarget;
    camera_up = camup;
    camera_direction = camera_target - camera_position;
    camera_direction = camera_direction / norm(camera_direction);
    camera_right = cross(camera_direction, camera_up);
    camera_right = camera_right / norm(camera_right);
    camera_up = cross(camera_right, camera_direction);
    DownEdge = [Xmin, Ymin, Zmin, 1]';
    UpEdge = [Xmax, Ymax, Zmax, 1]';
    LeftEdge = [Xmin, Ymax, Zmin, 1]';
    RightEdge = [Xmax, Ymin, Zmin, 1]';
    View = [camera_right, -dot(camera_right, camera_position);
        camera_up, -dot(camera_up, camera_position);
        -camera_direction, dot(camera_direction, camera_position);
        0, 0, 0, 1];
    DownEdge = View*DownEdge;
    UpEdge = View*UpEdge;
    LeftEdge = View*LeftEdge;
    RightEdge = View*RightEdge;
    FigWidth = RightEdge(1) - LeftEdge(1);
    FigHeight = UpEdge(2) - DownEdge(2);
    % AxisPosition(3) = min(sqrt((Xmax - Xmin)^2 + (Ymax - Ymin)^2)/(Zmax - Zmin)*AxisPosition(4), AxisPosition(3));
    if FigWidth < FigHeight || 1
        AxisPosition(3) = FigWidth/FigHeight*AxisPosition(4);
        c.Units = "centimeters";
        c.AxisLocationMode = "manual";
        ColorbarPosition = c.Position;
        ColorbarPosition(1) = AxisPosition(1) + AxisPosition(3) + 0.2;
        ColorbarPosition(3) = 0.55;
        ColorbarPosition(2) = AxisPosition(2);
        ColorbarPosition(4) = AxisPosition(4);
        c.Position = ColorbarPosition;
        ax.Position = AxisPosition;
    else
        % AxisPosition(4) = FigHeight/FigWidth*AxisPosition(3);
    end
end
hold off
end