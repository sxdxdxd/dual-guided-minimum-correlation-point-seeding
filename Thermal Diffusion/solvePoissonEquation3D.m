function T = solvePoissonEquation3D(x, y, z, phi)
    %% 每个方向格子的数量
    %% The number of grids in each direction
    Nx = size(x, 2) - 1;
    Ny = size(x, 1) - 1;
    Nz = size(x, 3) - 1;
    %% 每个方向格子的边长
    %% The length of the grid in each direction
    dx = x(1, 2, 1) - x(1, 1, 1);
    dy = y(2, 1, 1) - y(1, 1, 1);
    dz = z(1, 1, 2) - z(1, 1, 1);
    %% 格心坐标
    %% Grid center coordinates
    X = x + dx/2;
    X = X(1:end - 1, 1:end - 1, 1:end - 1);
    Y = y + dy/2;
    Y = Y(1:end - 1, 1:end - 1, 1:end - 1);
    Z = z + dz/2;
    Z = Z(1:end - 1, 1:end - 1, 1:end - 1);
    %% 由格点处源函数值插值求格心处源函数值
    %% Find the source function value at the grid center from the source function value at the grid point
    Phi = interp3(x, y, z, phi, X, Y, Z);
    %% 在格心坐标上求解poisson方程
    %% Solve the Poisson equation in grid center coordinates
    T_GridCenter = solvePoisson3D(Phi, Nx, Ny, Nz, dx, dy, dz);
    %% 由格心处解函数值插值求格点处解函数值(由于不希望在边界上放置种子点，干脆不做外插)
    %% Find the solution function value at the grid point from the solution function value at the grid center
    %(Since we don't want to place seeds on the boundary, we simply don't do extrapolation)
    T_GridPoint = interp3(X, Y, Z, T_GridCenter, x(2:end - 1, 2:end - 1, 2:end - 1), y(2:end - 1, 2:end - 1, 2:end - 1), z(2:end - 1, 2:end - 1, 2:end - 1));
    T = gather(T_GridPoint);
end
function T = solvePoisson3D(Phi, Nx, Ny, Nz, dx, dy, dz)
%% 计算特征值
%% Calculate eigenvalues
kx = 0:Nx - 1;
ky = 0:Ny - 1;
kz = 0:Nz - 1;
lambda_x = 2 * (cos(pi * kx / Nx) - 1) / dx^2;
lambda_y = 2 * (cos(pi * ky / Ny) - 1) / dy^2;
lambda_z = 2 * (cos(pi * kz / Nz) - 1) / dz^2;
% lambda_x = (16.*cos(pi * kx / Nx) - cos(2*pi * kx / Nx) - 15) / 6 / dx^2;%Another possible eigenvalue
% lambda_y = (16.*cos(pi * ky / Ny) - cos(2*pi * ky / Ny) - 15) / 6 / dy^2;%Another possible eigenvalue
% lambda_z = (16.*cos(pi * kz / Nz) - cos(2*pi * kz / Nz) - 15) / 6 / dz^2;%Another possible eigenvalue
[Lambda_x, Lambda_y, Lambda_z] = meshgrid(lambda_x, lambda_y, lambda_z);
%% 对源场做DCT-II变换
%% Perform DCT-II transform on the source field
Phi_hat = dctn(Phi);
%% 除以特征值
%% Divide by the eigenvalue
Lambda = Lambda_x + Lambda_y + Lambda_z;
Phi_hat(1, 1, 1) = 0;
Lambda(1, 1, 1) = 1;
T_hat = Phi_hat ./ Lambda;
%% 对温度场做DCT-II逆变换
%% Perform DCT-II inverse transform on the temperature field
T = idctn(T_hat);
end