%% addpath
addpath(genpath('Common'));
addpath(genpath('Dataset'));
addpath(genpath('Tools For Graph'));
addpath(genpath('Thermal Diffusion'));
addpath(genpath('Tools For Evaluation'));
%% 清理
%% clear
clc;
clear;
close all;
if ~isempty(gcp('nocreate'))
    parfevalOnAll(@gpuDevice, 0, []);
else
    reset(gpuDevice);
end
%% 准备并行资源
%% Prepare parallel resources
GPU = gpuDevice();
if isempty(gcp('nocreate'))
    parpool('Threads');
    % parpool('Processes');
end
%% 读取数据集
%% load DataSets
%% 3D
% load back_step_flow.mat; X = cat(4, x, y, z); V = cat(4, u, v, w); ScalarField = max(-u, 0); alpha = 1; N_Seeds = 10;
% load cascade_impactor.mat; X = cat(4, x, y, z); V = cat(4, u, v, w); ScalarField = calculateOmega3D(x, y, z, u, v, w); alpha = 0; N_Seeds = 80;
% load cavity_flow.mat; X = cat(4, x, y, z); V = cat(4, u, v, w); ScalarField = calculateOmega3D(x, y, z, u, v, w); alpha = 0.0; N_Seeds = 70;
% load incense_stick.mat; X = cat(4, x, y, z); V = cat(4, u, v, w); ScalarField = calculateOmega3D(x, y, z, u, v, w); alpha = 0.0; N_Seeds = 50;
% load tangaroa.mat; X = cat(4, x, y, z); V = cat(4, u, v, w); ScalarField = calculateOmega3D(x, y, z, u, v, w); alpha = 0.0; N_Seeds = 50;
load wall_jet_mixer.mat; X = cat(4, x, y, z); V = cat(4, u, v, w); ScalarField = k.^2 ./ ep; alpha = 0.25; N_Seeds = 20;
% load water_evaporate.mat; X = cat(4, x, y, z); V = cat(4, u, v, w); ScalarField = Humidity; alpha = 0.25; N_Seeds = 30;
% load water_purification_reactor.mat; X = cat(4, x, y, z); V = cat(4, u, v, w); ScalarField = k.^2 ./ ep; N_Seeds = 30; alpha = 0.5;
%% 2D
% load darcy.mat; X = cat(3, x, y); V = cat(3, u, v); ScalarField = c; alpha = 0.3; N_Seeds = 10;
% load piped_cylinder.mat; X = cat(3, x, y); V = cat(3, u, v); ScalarField = calculateOmega2D(x, y, u, v); alpha = 0; N_Seeds = 40;
% load tesla_microvalve.mat; X = cat(3, x, y); V = cat(3, u, v); ScalarField = p; alpha = 0.9; N_Seeds = 20; %tesla microvalve:forward flow
% load tesla_microvalve.mat; X = cat(3, x, y); V = cat(3, u2, v2); ScalarField = p2; alpha = 0.9; N_Seeds = 20;%tesla microvalve:reverse flow
%% 设置参数
%% Set parameters
ScalarField(isnan(ScalarField)) = 0;
ScalarField(isinf(ScalarField)) = 0;
Options_Common = getBasicParameters(X, V, N_Seeds);
Options_Common.alpha = alpha;% Phi的重要性
IsThreeDimensionalData = ndims(ScalarField) == 3;
if IsThreeDimensionalData
    x = X(:, :, :, 1); y = X(:, :, :, 2); z = X(:, :, :, 3);
    u = V(:, :, :, 1); v = V(:, :, :, 2); w = V(:, :, :, 3);
    OmegaField = calculateOmega3D(x, y, z, u, v, w);
else
    x = X(:, :, 1); y = X(:, :, 2);
    u = V(:, :, 1); v = V(:, :, 2);
    OmegaField = calculateOmega2D(x, y, u, v);
end
Options_Ours = getParametersForOurMethod(Options_Common);
%% 生成流线
%% Generate streamlines
tic;
StreamLines_Ours = generateStreamlineByOurMethod(X, V, ScalarField, OmegaField, Options_Ours);%dual-guided minimum correlation point seeding Method
Time = toc;
disp(['使用双引导的最小相关点播种法生成了', num2str(size(StreamLines_Ours, 1)), '根流线，共用时', num2str(Time,'%.4f'), '秒,', ...
    '平均每根流线用时', num2str(Time / size(StreamLines_Ours, 1),'%.4f'), '秒']);
disp(['Using the double-guided minimum correlation point seeding method, ', ...
    num2str(size(StreamLines_Ours, 1)), ' streamlines were generated, which took a total of ', num2str(Time,'%.4f'), ' seconds,', ...
'with an average time of ', num2str(Time / size(StreamLines_Ours, 1),'%.4f'), ' seconds per streamline']);
%% 评估流线质量
%% Evaluate streamline quality
ExpressionData = calculateExpressivenessForStreamLines(X, V, ScalarField, cell2mat(StreamLines_Ours(:, 1)), Options_Common);
disp(['Expressiveness: ', num2str(ExpressionData.Expressiveness, '%.5f')]);
disp(['ConditionalEntropy: ', num2str(ExpressionData.ConditionalEntropy, '%.5f')]);
disp(['PSNR of Flow Field: ', num2str(ExpressionData.PSNR_V, '%.5f')]);
disp(['PSNR of Scalar Field ', num2str(ExpressionData.PSNR_Phi, '%.5f')]);
disp(['AAD (Average Angle Differenc): ', num2str(ExpressionData.AAD, '%.5f')]);
disp(['ARE (Average Relative Error): ', num2str(ExpressionData.ARE, '%.5f')]);
disp(' ');
%% 绘制流线
%% Draw streamlines
close all
DrawSeed = 0;
MinAlpha = 0.2;
CLabel = '';
ColorField = ScalarField;
OpenColorBar = 0;
AdjustIntegralConstant = 1;
f = drawStreamlines(StreamLines_Ours, X, ColorField, MinAlpha, DrawSeed, CLabel, OpenColorBar, AdjustIntegralConstant);