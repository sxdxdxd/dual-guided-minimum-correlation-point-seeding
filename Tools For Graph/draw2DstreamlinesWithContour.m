function f = draw2DstreamlinesWithContour(StreamLines, x, y, ScalarField, Seed, CLabel, OpenColorBar, AdjustIntegralConstant)
if min(ScalarField, [], 'all') < 0 && AdjustIntegralConstant
    ScalarField = ScalarField + abs(min(ScalarField, [], 'all'));
end
f = figure;
hold on
axis tight
axis equal
grid off
box on
ax = gca;
ax.BoxStyle = 'full';
contourf(x, y, ScalarField, 'LineStyle', 'none');
if ~isempty(StreamLines)
    lineobj = streamline(StreamLines(:, 1));
    [lineobj(:).Color] = deal('k');
end
if strcmp(Seed ,'On')
    Seeds = cell2mat(StreamLines(:, 2));
    scatter(Seeds(:, 1), Seeds(:, 2), 'red', 'filled', 'hexagram');
end
load('cool2warm.mat', 'map');
colormap(map);
% xlabel('x');
% ylabel('y');
set(gca, 'fontsize', 16, 'FontWeight', 'bold');%设置坐标轴字体大小
set(gca, 'linewidth', 1.5, 'FontWeight', 'bold'); 
% set(get(gca, 'XLabel'), 'FontSize', 16, 'FontWeight','bold');
% set(get(gca, 'YLabel'), 'FontSize', 16, 'FontWeight','bold');
% set(get(gca, 'TITLE'), 'FontSize', 16, 'FontWeight','bold');
set(gca, 'xtick', [], 'xticklabel', [], 'XLabel', [])
set(gca, 'ytick', [], 'yticklabel', [], 'YLabel', [])
MaxD = max(ScalarField, [],'all');
MinD = min(ScalarField, [],'all');
% MinD = 0;
clim([MinD, MaxD]);
Ticks = clim;
Ticks = linspace(Ticks(1), Ticks(2), 6);
TickLabels = cell(size(Ticks));
for i = 1:size(Ticks, 2)
    % TickLabels{i} = num2str(Ticks(i)); 
    TickLabels{i} = erase(sprintf('%3.2g', Ticks(i)), '+'); 
end
if(OpenColorBar)
    c = colorbar('Ticks', Ticks);
end
% if isempty(CLabel)
%     c.Label.String = 'User defined scalar field';
% else
%     c.Label.String = CLabel;
% end
c.Label.String = "";
c.TickLabels = TickLabels;
shading interp
hold off
end