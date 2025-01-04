function Binfield = calculateBinField2D(u, v, N_Bins)
assert(N_Bins <= intmax('int8'));
%% 生成角度场
%% Generate angle field
Theta = atan2(v, u);
N_Y = size(Theta, 1);
N_X = size(Theta, 2);
%% 计算Bin场
%% Calculate Bin field
Bin = linspace(-pi - 1e-15, pi + 1e-15, N_Bins + 1);%-pi + 2*pi.*NN./60
Bin = Bin(2:end);
Bin_Field = zeros(N_Y, N_X, 'int8');
parfor Iy = 1:N_Y
    for Ix = 1:N_X
        Binindex = lower_bound(Bin, Theta(Iy, Ix));
        Bin_Field(Iy, Ix) = Binindex;
    end
end
Binfield = Bin_Field;
end
