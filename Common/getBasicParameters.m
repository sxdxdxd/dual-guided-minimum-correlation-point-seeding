function Options_Common = getBasicParameters(X, V, N_Seeds)
Options_Common.Dims = ndims(X) - 1;
Options_Common.N_x = size(X, 2);%The number of grids along the x-axis %沿x轴方向的网格数量
Options_Common.N_y = size(X, 1);%The number of grids along the y-axis %沿y轴方向的网格数量
Options_Common.alpha = 0.5;% Importance of scalar fields %标量场的重要性
if ndims(X) == 4
    Stepsize = 0.1;%Integral step size when generating streamlines (compared to one cell) %生成流线时的积分步长(与一个元胞相比)
    Options_Common.N_z = size(X, 3);%The number of grids along the z-axis %沿z轴方向的网格数量
    Max_Vertices = ceil((Options_Common.N_x + Options_Common.N_y + Options_Common.N_z) / 3 / Stepsize) * 10;%单次计算流线的最大积分步数
    %The maximum number of integration steps for a single streamline
else
    Stepsize = 0.5;%生成流线时的积分步长(与一个元胞相比)
    Max_Vertices = ceil((Options_Common.N_x + Options_Common.N_y) / 2 / Stepsize) * 10;
end
%% 生成流线时所用的参数
%% Parameters used when generating streamlines
Options_Common.N_Seeds = N_Seeds;% Total number of seed points, i.e. the number of streamlines generated in the end %种子点总数,即最终生成流线的数量
Options_Common.Stepsize = Stepsize;
Options_Common.Max_Vertices =  Max_Vertices;
%% 将流线重建为流场时所用的参数
%% Parameters used to reconstruct streamlines into flow fields
Options_Common.Mu = 0.1;%This is the weight coefficient of the smoothing term %此为平滑项的权重系数
Options_Common.DiffusionThreshold = 0.995;
Options_Common.OptimizeMaxIterations = 9 * max(size(X));%This is the maximum number of iterations when streamlining diffusion %此为流线扩散时的最大迭代次数
Options_Common.OvercorrectionFactor = 1;
% if ndims(X) == 4
%     Options_Common.OvercorrectionFactor = 1.94;
% else
%     Options_Common.OvercorrectionFactor = 1;
% end
end