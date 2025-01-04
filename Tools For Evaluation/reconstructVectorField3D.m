function [U, V, W] = reconstructVectorField3D(Mask, u, v, w, OvercorrectionFactor, Mu, DiffusionThreshold, MaxIterations)
%% 基于最优化方法的流线扩散
%% Streamline diffusion based on optimization method
u_Marked = Mask .* u;
v_Marked = Mask .* v;
w_Marked = Mask .* w;
[U, V, W] = OptimizationFunc(u_Marked, v_Marked, w_Marked, OvercorrectionFactor, Mu, DiffusionThreshold, MaxIterations);
end
function [U_r, V_r, W_r] = OptimizationFunc(u_Marked, v_Marked, w_Marked, OvercorrectionFactor, Mu, DiffusionThreshold, MaxIterations)
N_x = size(u_Marked, 2);
N_y = size(u_Marked, 1);
N_z = size(w_Marked, 3);
%% 向量场归一化
%% Vector field normalization
Length = sqrt(u_Marked.^2 + v_Marked.^2 + w_Marked.^2);
Length(Length == 0) = 1;
U_Temp = u_Marked ./ Length;
V_Temp = v_Marked ./ Length;
W_Temp = w_Marked ./ Length;
%% 计算系数
%% Calculation coefficient
OffsetWeight_w = U_Temp.^2 + V_Temp.^2 + W_Temp.^2;
OffsetWeight_x = OffsetWeight_w .* U_Temp;
OffsetWeight_y = OffsetWeight_w .* V_Temp;
OffsetWeight_z = OffsetWeight_w .* W_Temp;
%% 迭代优化
%% Iterative optimization
PrevDiffusionError = intmax;
[X, Y, Z] = meshgrid(1:N_x, 1:N_y, 1:N_z);
Ix = gpuArray(X);
Iy = gpuArray(Y);
Iz = gpuArray(Z);
for iIter = 0:MaxIterations
    [U_R, V_R, W_R] = arrayfun(@updateGrid, Ix, Iy, Iz, N_x, N_y, N_z);
    DiffusionError = sum((U_Temp - U_R).^2 + (V_Temp - V_R).^2 + (W_Temp - W_R).^2, 'all');
    DiffusionError = sqrt(DiffusionError / numel(U_R));
    if iIter > max([N_x, N_y, N_z])
        ErrorRate = DiffusionError / PrevDiffusionError;
        if DiffusionThreshold <= ErrorRate && ErrorRate <= 1.0
            break;
        end
    end
    if( iIter > 0 && PrevDiffusionError < DiffusionError)
        break;
    end
    PrevDiffusionError  = DiffusionError;
    U_Temp = U_R;
    V_Temp = V_R;
    W_Temp = W_R;
end
U_r = gather(U_R);
V_r = gather(V_R);
W_r = gather(W_R);
    function [u_r, v_r, w_r] = updateGrid(Ix, Iy, Iz, Nx, Ny, Nz)
        XLeftIndex = max(1, Ix - 1); XRightIndex = min(Ix + 1, Nx);
        YBottomIndex = max(1, Iy - 1); YTopIndex = min(Iy + 1, Ny);
        ZFrontIndex = max(1, Iz - 1); ZBackIndex = min(Iz + 1, Nz);
        u_r = (1 - OffsetWeight_w(Iy, Ix, Iz)) .* U_Temp(Iy, Ix, Iz) + ...
            Mu .* (U_Temp(Iy, XLeftIndex, Iz) + U_Temp(Iy, XRightIndex, Iz) + ...
            U_Temp(YBottomIndex, Ix, Iz) + U_Temp(YTopIndex, Ix, Iz) + ...
            U_Temp(Iy, Ix, ZFrontIndex) + U_Temp(Iy, Ix, ZBackIndex) - 6*U_Temp(Iy, Ix, Iz)) +...
            OvercorrectionFactor * OffsetWeight_x(Iy, Ix, Iz);
        v_r = (1 - OffsetWeight_w(Iy, Ix, Iz)) .* V_Temp(Iy, Ix, Iz) + ...
            Mu .* (V_Temp(Iy, XLeftIndex, Iz) + V_Temp(Iy, XRightIndex, Iz) + ...
            V_Temp(YBottomIndex, Ix, Iz) + V_Temp(YTopIndex, Ix, Iz) + ...
            V_Temp(Iy, Ix, ZFrontIndex) + V_Temp(Iy, Ix, ZBackIndex) - 6*V_Temp(Iy, Ix, Iz)) +...
            OvercorrectionFactor * OffsetWeight_y(Iy, Ix, Iz);
        w_r = (1 - OffsetWeight_w(Iy, Ix, Iz)) .* W_Temp(Iy, Ix, Iz) + ...
            Mu .* (W_Temp(Iy, XLeftIndex, Iz) + W_Temp(Iy, XRightIndex, Iz) + ...
            W_Temp(YBottomIndex, Ix, Iz) + W_Temp(YTopIndex, Ix, Iz) + ...
            W_Temp(Iy, Ix, ZFrontIndex) + W_Temp(Iy, Ix, ZBackIndex) - 6*W_Temp(Iy, Ix, Iz)) +...
            OvercorrectionFactor * OffsetWeight_z(Iy, Ix, Iz);
    end
end