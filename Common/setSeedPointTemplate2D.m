function L = setSeedPointTemplate2D(x, y, ExtremePointIdList, Xmin, Xmax, Ymin, Ymax)
%% 初始化参数
%% Initialization parameters
delta_x = x(1, 2) - x(1, 1);
delta_y = y(2) - y(1);
N = size(ExtremePointIdList, 1);
SeedPointList = zeros(9 * N, 2);
%% 在标量场的极大值点处放置种子点模板
%% Place seed templates at the maximum point of the scalar field
for i = 1:N
    [X, Y] = SeedPointTemplate2D(x(ExtremePointIdList(i)), y(ExtremePointIdList(i)), delta_x, delta_y);
    SeedPointList(9 * (i - 1) + 1:9 * i, 1) = X;
    SeedPointList(9 * (i - 1) + 1:9 * i, 2) = Y;
end
%% 去除超出边界范围的种子点
%% Remove seeds that are beyond the boundary
SeedPointList(SeedPointList(:, 1) < Xmin, :) = [];
SeedPointList(SeedPointList(:, 1) > Xmax, :) = [];
SeedPointList(SeedPointList(:, 2) < Ymin, :) = [];
SeedPointList(SeedPointList(:, 2) > Ymax, :) = [];
L = SeedPointList;
end
%% function 1
function [X, Y] = SeedPointTemplate2D(x, y, delta_x, delta_y)
x0 = [0, -1, 1, -2, 0, 2, -1, 1, 0];
y0 = [2, 1, 1, 0, 0, 0, -1, -1, -2];
X = x0 .* delta_x + x;
Y = y0 .* delta_y + y;
end