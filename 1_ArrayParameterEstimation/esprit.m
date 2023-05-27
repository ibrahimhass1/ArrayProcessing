function theta = esprit(X,d)

    % Redefine X 
    Y = X(2:end, :);
    X = X(1:end-1, :);
    
    Z = [X;Y];
    M = size(Z, 1)/2;
    [U, E, V] = svd(Z);
    U_hat = U(:, 1:d);
    U_x = U_hat(1:M, :);
    U_y = U_hat(M+1:end, :);
    
    
    phi = eig(inv(U_x'*U_x)*U_x'*U_y);
    
    angles = angle(phi);
    
    theta = 180/pi*asin(angles./pi);

end