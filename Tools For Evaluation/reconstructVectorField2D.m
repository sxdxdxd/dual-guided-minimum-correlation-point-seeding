function [U, V] = reconstructVectorField2D(Mask, u, v, OvercorrectionFactor, Mu, DiffusionThreshold, MaxIterations)
%% 基于最优化方法的流线扩散
%% Streamline diffusion based on optimization method
u_Marked = Mask .* u;
v_Marked = Mask .* v;
[U, V] = OptimizationFunc(u_Marked, v_Marked, OvercorrectionFactor, Mu, DiffusionThreshold, MaxIterations);
end
function [U_r, V_r] = OptimizationFunc(u_Marked, v_Marked, OvercorrectionFactor, Mu, DiffusionThreshold, MaxIterations)
N_x = size(u_Marked, 2);
N_y = size(u_Marked, 1);
%% 向量场归一化
%% Vector field normalization
Length = sqrt(u_Marked.^2 + v_Marked.^2);
Length(Length == 0) = 1;
U_Temp = u_Marked ./ Length;
V_Temp = v_Marked ./ Length;
%% 计算系数
%% Calculation coefficient
OffsetWeight_w = U_Temp.^2 + V_Temp.^2;
OffsetWeight_x = OffsetWeight_w .* U_Temp;
OffsetWeight_y = OffsetWeight_w .* V_Temp;
%% 迭代优化
%% Iterative optimization
PrevDiffusionError = intmax;
[X, Y] = meshgrid(1:N_x, 1:N_y);
Ix = gpuArray(X);
Iy = gpuArray(Y);
for iIter = 0:MaxIterations
    [U_R, V_R] = arrayfun(@updateGrid, Ix, Iy, N_x, N_y);
    DiffusionError = sum((U_Temp - U_R).^2 + (V_Temp - V_R).^2, 'all');
    DiffusionError = sqrt(DiffusionError / numel(U_R));
    if iIter > max([N_x, N_y])
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
end
U_r = gather(U_R);
V_r = gather(V_R);
    function [u_r, v_r] = updateGrid(Ix, Iy, Nx, Ny)
        XLeftIndex = max(1, Ix - 1); XRightIndex = min(Ix + 1, Nx);
        YBottomIndex = max(1, Iy - 1); YTopIndex = min(Iy + 1, Ny);
        u_r = (1 - OffsetWeight_w(Iy, Ix)) * U_Temp(Iy, Ix)...
            + Mu * (U_Temp(Iy, XLeftIndex) + U_Temp(Iy, XRightIndex) + U_Temp(YBottomIndex, Ix) + U_Temp(YTopIndex, Ix) - 4*U_Temp(Iy, Ix))...
            + OvercorrectionFactor * OffsetWeight_x(Iy, Ix);
        v_r = (1 - OffsetWeight_w(Iy, Ix)) * V_Temp(Iy, Ix)...
            + Mu * (V_Temp(Iy, XLeftIndex) + V_Temp(Iy, XRightIndex) + V_Temp(YBottomIndex, Ix) + V_Temp(YTopIndex, Ix) - 4*V_Temp(Iy, Ix))...
            + OvercorrectionFactor * OffsetWeight_y(Iy, Ix);
    end
end