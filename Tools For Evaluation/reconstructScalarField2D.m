function ScalarField_R = reconstructScalarField2D(X, Mask, ScalarField)
Index = (Mask > 0 );
x = X(:, :, 1); y = X(:, :, 2);
Index([1, end], [1, end]) = 1;
IndexNegative = ~Index;
Mask_Scalar = ScalarField(Index);
Mask_X = x(Index);
Mask_Y = y(Index);
S_r = griddata(Mask_X, Mask_Y, double(Mask_Scalar), x(IndexNegative), y(IndexNegative));
S_R = zeros(size(x));
S_R(Index) = Mask_Scalar;
S_R(IndexNegative) = S_r;
Smax = max(ScalarField, [], 'all');
Smin = min(ScalarField, [], 'all');
S_R(S_R > Smax) = Smax;
S_R(S_R < Smin) = Smin;
ScalarField_R = S_R;
end