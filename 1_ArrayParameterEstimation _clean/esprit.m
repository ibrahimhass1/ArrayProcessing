%% Parameters:
%
% X: MxN received signal matrix, where each row is the signal of one of the
% antennas
% d: Number of sources

%% Output
%
% dx1 vector with the estimates of the AoA in degrees for each the d sources

function theta = esprit(X,d)

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
    
    % Calcluate AoA using the estimates of phi (Report eq. 24)
    angles = angle(phi);
    theta = 180/pi*asin(angles./pi);
end