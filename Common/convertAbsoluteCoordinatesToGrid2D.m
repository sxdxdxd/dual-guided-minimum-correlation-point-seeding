function [IndexPassedX, IndexPassedY] = convertAbsoluteCoordinatesToGrid2D(PointsPassedByStreamLine, N_x, Xmax, Xmin, N_y, Ymax, Ymin)
%% 去除超出边界范围的点
%% Remove points that are beyond the boundary
PointsPassedByStreamLine(PointsPassedByStreamLine(:, 1) < Xmin, :) = [];
PointsPassedByStreamLine(PointsPassedByStreamLine(:, 1) > Xmax, :) = [];
PointsPassedByStreamLine(PointsPassedByStreamLine(:, 2) < Ymin, :) = [];
PointsPassedByStreamLine(PointsPassedByStreamLine(:, 2) > Ymax, :) = [];
%% 变绝对坐标为网格坐标
%% Convert absolute coordinates to grid coordinates
IndexPassedX = floor((PointsPassedByStreamLine(:, 1) - Xmin) ./ (Xmax - Xmin) .* (N_x - 1) + 1);
IndexPassedY = floor((PointsPassedByStreamLine(:, 2) - Ymin) ./ (Ymax - Ymin) .* (N_y - 1) + 1);
end