function f = espritfreq(X,d)

% Obtain single vector to construct X
x = X(1, :);

m = 5;
X=zeros([m length(x)-m]);
for i=1:length(x)-m
    X(:, i) = transpose(x(i:i+(m-1)));
end

% Redefine X 
Y = X(2:end, :);
X = X(1:end-1, :);

Z = [X;Y];
M = size(Z, 1)/2;
[U, E, V] = svd(Z);
%[U_hat, E, V] = svd(Z,'econ')
U_hat = U(:, 1:d);
U_x = U_hat(1:M, :);
U_y = U_hat(M+1:end, :);


phi = eig(inv(U_x'*U_x)*U_x'*U_y);
%phi = eig(pinv(U_x)*U_y);
angles = angle(phi);

f = angles./(2*pi);
%phi = 2*pi*sin(theta)
