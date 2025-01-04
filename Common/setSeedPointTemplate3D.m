function L = setSeedPointTemplate3D(x, y, z, ExtremePointIdList, Xmin, Xmax, Ymin, Ymax, Zmin, Zmax)
%% 初始化参数
%% Initialization parameters
delta_x = x(1, 2, 1) - x(1, 1, 1);
delta_y = y(2, 1, 1) - y(1, 1, 1);
delta_z = z(1, 1, 2) - z(1, 1, 1);
N = size(ExtremePointIdList, 1);
SeedPointList = zeros(27 * N, 3);
%% 在标量场的极大值点处放置种子点模板
%% Place seed templates at the maximum point of the scalar field
for i = 1:N
    j = ExtremePointIdList(i);
    [X, Y, Z] = SeedPointTemplate3D(x(j), y(j), z(j), delta_x, delta_y, delta_z);
    SeedPointList(size(X, 2) * (i - 1) + 1:size(X, 2) * i, 1) = X;
    SeedPointList(size(X, 2) * (i - 1) + 1:size(X, 2) * i, 2) = Y;
    SeedPointList(size(X, 2) * (i - 1) + 1:size(X, 2) * i, 3) = Z;
end
%% 去除超出边界范围的种子点
%% Remove seeds that are beyond the boundary
SeedPointList(SeedPointList(:, 1) < Xmin, :) = [];
SeedPointList(SeedPointList(:, 1) > Xmax, :) = [];
SeedPointList(SeedPointList(:, 2) < Ymin, :) = [];
SeedPointList(SeedPointList(:, 2) > Ymax, :) = [];
SeedPointList(SeedPointList(:, 3) < Zmin, :) = [];
SeedPointList(SeedPointList(:, 3) > Zmax, :) = [];
L = SeedPointList;
end
%% function 1
function [X, Y, Z] = SeedPointTemplate3D(x, y, z, delta_x, delta_y, delta_z)
%27 octahedral
x_positive = [0, -3,  0,  3,  0, -2,  2,  2, -2];
y_positive = [0,  0, -3,  0,  3, -2, -2,  2,  2];
z_positive = [6,  3,  3,  3,  3,  2,  2,  2,  2];
x0 = [0, -3,  3, -6,  0,  6, -3,  3,  0];
y0 = [6,  3,  3,  0,  0,  0, -3, -3, -6];
z0 = [0,  0,  0,  0,  0,  0,  0,  0,  0];

%27 square
% x_positive = [0, -3,  0,  3,  0, 3,  3,  -3, -3];
% y_positive = [0,  0, -3,  0,  3, 3, -3,  3,  -3];
% z_positive = [3,  3,  3,  3,  3,  3,  3,  3,  3];
% x0 = [0, -3,  0,  3,  0, 3,  3,  -3, -3];
% y0 = [0,  0, -3,  0,  3, 3, -3,  3,  -3];
% z0 = [0,  0,  0,  0,  0,  0,  0,  0,  0];

xx = [x_positive, x0, -x_positive];
yy = [y_positive, y0, -y_positive];
zz = [z_positive, z0, -z_positive];
X = xx .* delta_x + x;
Y = yy .* delta_y + y;
Z = zz .* delta_z + z;
end