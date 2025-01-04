function Binfield = calculateBinField3D(U, V, W, BinMap)
%% 生成角度场
%% Generate angle field
assert(sum(~isinf(BinMap(:, 1:end - 1, 1)), "all") <= intmax('int16'));
[Theta, Phi, ~] = cart2sph(U, V, W);
Theta = Theta + pi;
Phi = Phi + pi / 2;
N_Y = size(Theta, 1);
N_X = size(Theta, 2);
N_Z = size(Theta, 3);
%% 计算Bin场
%% Calculate Bin field
Bin_Field = zeros(N_Y, N_X, N_Z, 'int16');
for Iz = 1:N_Z
    for Iy = 1:N_Y
        for Ix = 1:N_X
            Binindex = locateVectorInBins(Theta(Iy, Ix, Iz), Phi(Iy, Ix, Iz), BinMap);
            Bin_Field(Iy, Ix, Iz) = Binindex;
        end
    end
end
Binfield = Bin_Field;
end
%% function 定位向量在直方图中的位置
%% locates the position of the vector in the histogram
function BinIndex = locateVectorInBins(theta, phi, BinMap)
Index_Phi = lower_bound(BinMap(:, end, 1), phi);
Index_Theta = lower_bound(BinMap(Index_Phi, 1:end - 1, 1), theta);
if BinMap(Index_Phi, Index_Theta, 2) > 360
    Index_Theta = Index_Theta - 1;
end
BinIndex = BinMap(Index_Phi, Index_Theta, 2);
end