function f = ESPRIT_freq(X,d)
% ESPRIT algorithm for estimating the frequencies
[M,N] = size(X);

m = floor(N/3);
%m = length(X)/2;

% Modify X

% 1
% X1 = X(:, 1:m);
% X2 = X(:, m+1:end);
% X_mod = [X1; flipud(X2)];

% 2
% x_m = X(1,:);
% x_m = x_m';
% 
% for r=1:m
%    for c=1:m
%       X_mod(r,c) = x_m(r+c-1); 
%    end
% end

% 3
Z = [];

for i=1:1
    x = X(i,:);
    %X_mod = hankel(X(i,1:N-m),X(i,m:N));
    hankel = zeros([m N-m]);
    for j = 1:(N-m)
       hankel(:,j) = transpose(x(j:j+m-1)); 
    end

%         m = 5;
%     X=zeros([m length(x)-m]);
%     for i=1:length(x)-m
%         X(:, i) = transpose(x(i:i+(m-1)));
%     end  
    
    A = hankel(2:end, :);
    B = hankel(1:end-1, :);
    
    Z = [B;A];
%     disp(size(X))
%     disp(size(hankel))
%     Z = [X;hankel];
%     disp(size(Z))
end

% 4
% X_mod = zeros(m,m);
% for n=1:N-m+1
%    X_mod = X_mod+(X(2,n:n+m-1))*(X(2,n:n+m-1))';
% end

% 5
% Z = zeros(m,N-m+1);
% for k = 1:N-m+1
%     Z(:,k) = X(1:m,k:k+m-1);
% end

%6


% SVD of X
[U,~,~] = svd(Z);
M = size(Z, 1)/2;
% Signal subspace with first d elements
%Ux = U(1:floor(M/2)-1,1:d);
%Uy = U(2:floor(M/2),1:d);

U_hat = U(:, 1:d);
Ux = U_hat(1:M, :);
Uy = U_hat(M+1:end, :);

D = pinv(Ux)*Uy;

% SVD of D to find DOA and  T
[~, eigenval] = eig(D);
f = angle(sort(diag(eigenval)))./(2*pi);

%f = mod(f,1);
end

