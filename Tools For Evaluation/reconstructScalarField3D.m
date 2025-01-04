function ScalarField_R = reconstructScalarField3D(X, Mask, ScalarField)
Index = (Mask > 0 );
x = X(:, :, :, 1); y = X(:, :, :, 2); z = X(:, :, :, 3);
Index([1, end], [1, end], [1, end]) = 1;
IndexNegative = ~Index;
Mask_Scalar = ScalarField(Index);
Mask_X = x(Index);
Mask_Y = y(Index);
Mask_Z = z(Index);
S_r = griddata(Mask_X, Mask_Y, Mask_Z, double(Mask_Scalar), x(IndexNegative), y(IndexNegative), z(IndexNegative));
S_R = zeros(size(x));
S_R(Index) = Mask_Scalar;
S_R(IndexNegative) = S_r;
Smax = max(ScalarField, [], 'all');
Smin = min(ScalarField, [], 'all');
S_R(S_R > Smax) = Smax;
S_R(S_R < Smin) = Smin;
ScalarField_R = S_R;
end