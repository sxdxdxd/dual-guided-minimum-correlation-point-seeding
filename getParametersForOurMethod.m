function Options_Ours = getParametersForOurMethod(Options_Common)
%% OurMethod主要参数
Options_Ours.N_Seeds = Options_Common.N_Seeds;%种子点总数,即最终生成流线的数量
Options_Ours.alpha = Options_Common.alpha;%alpha
Options_Ours.N_threshold = 0;%放置种子点Log（N）阶段与N阶段的分界点(剩余种子数)
if Options_Common.Dims == 3
    Options_Ours.UseTemplateInit = true;
    Options_Ours.UseTemplateIter = true;
    Options_Ours.InitialRatio = 0.1;% beta 即N_Seeds中最少有多少比例的流线初始流线，需小于1
    Options_Ours.Proportion = 0.1;%gamma Log（N）阶段时，单次放置种子点占剩余种子点的比例
else
    Options_Ours.UseTemplateInit = true;
    Options_Ours.UseTemplateIter = true;
    Options_Ours.InitialRatio = 0;
    Options_Ours.Proportion = 0;
end
Options_Ours.UseMask_01 = false;
Options_Ours.UseMask_0N = false;
Options_Ours.UseMask_MN = false;
Options_Ours.UseMaskGaussian = true;
Options_Ours.KernelSize = 5;
Options_Ours.sigma = 1;
%% 处理标量场中负数的方法
Options_Ours.handlingNegative = 'AddConst';
%'AddConst'方法在标量场中给每个点加一个常数，使得场中最小的值不再为负；
%'Abs'方法在标量场中给每个点取绝对值
%% 在标量场中搜寻极大值点所用的参数
Options_Ours.MaximumAspectRatio = sqrt(2);%使用KDTree方法划分场时（LeafNode之前的层）每个块的最大长宽比，需大于等于1；
%% 生成流线时所用的参数
Options_Ours.Stepsize = Options_Common.Stepsize;
Options_Ours.Max_Vertices = Options_Common.Max_Vertices;%单次计算流线的最大积分步数，即单段流线的最大点数
end