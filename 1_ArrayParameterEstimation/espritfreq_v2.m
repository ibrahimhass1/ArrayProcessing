%% Parameters:
%
% X: MxN received signal matrix, where each row is the signal of one of the
% antennas
% d: Number of sources

%% Output
%
% dx1 vector with the estimates of the normalized frequecies for each 
% the d sources

function f = espritfreq(X,d)
    
    % Obtain single vector to construct X
    x = X
    
    % Construct data matrix with 2 times more columns than rows 
    m = floor(size(X, 2)/3);
    m
    M = size(X, 1);
    X= zeros([m (size(x,2)-m)*M]);
    size(X)
    % Construct X (Report eq. 26)
    for i=1:size(x,2)-m
        (i-1)*M+1
        (i+1)*M
        X(:, (i-1)*M+1:i*M) = transpose(x(:,i:i+m-1));
    end
    X
    % Construct X and Y from data (Report eq. 14, 15)
    X = X(1:end-1, :); 
    Y = X(2:end, :);

    % Construct Z (Report eq. 17)
    Z = [X;Y];
    
    % Take svd
    [U, E, V] = svd(Z); 

    % Take part of U that spans the columns space of Z (Report eq. 18)
    U_z = U(:, 1:d);

    % Obtains Ux and Uy using Uz
    M = size(Z, 1)/2;
    U_x = U_z(1:M, :);
    U_y = U_z(M+1:end, :);

    % Find phi's using evd (Report eq. 22)
    phi = eig(inv(U_x'*U_x)*U_x'*U_y);
    angles = angle(phi);
    
    % Using the estimates, calculate the frequencies (Report eq. 29)
    f = angles./(2*pi);
end
