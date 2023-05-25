%% Q1
function [X,A,S] = genDataRK(M,N,Delta,theta,f,SNR)

% Check if theta's are within the specified range
if or((min(theta) < -90), (max(theta) >= 90))
  error('Source direction is out of range');
end

% Check if normalized frequencies are within the specified range
if or((min(f) < 0), (max(f) >= 1))
  error('Source frequency is out of range');
end

% Obtain number of sources
d = length(theta);

% Initialize source matrix, dimensions dxN
S = zeros([N d]);
% Construct the source matrix for N samples,
S =  transpose(exp(1i*2*pi*f*(0:N-1)));

% Calculate corresponding phase shift, note the the sin function expects
% the angle to be in radians
phi = 2*pi*Delta*sin(deg2rad(theta));

% Initialize array response matrix, dimensions Mxd
A = zeros([M d]);
% Construct the array resonse matrix for M receivers
A = exp(1i*(0:M-1)'.*phi');

% Use Khatri-Roa to construct X from A and S
X = zeros([M*N d]);
for i = 1:d
    % Perform column-wise Kronecker product, which is the Khatri-Roa
    % product
    X(:, i) = kron(A(:, i), S(:, i));
end

X = X*ones([d 1]);

% Unvec the stacked output vectors
X =transpose(reshape(X, N, M));


% Estimate received signal power for each antenna
powerPerSource = mean(abs(X).^2, 1);
% Calculate desired noise level to obatain the given SNR per antenna
scaledNoisePower = powerPerSource./10^(SNR/10);

% Construct complex-valued additive noise matrix 
Noise = zeros([M N]);
for i=1:M
    Noise(i, :) = wgn(1, N, scaledNoisePower(i), 'linear', 'complex');
end

% The measurement matrix is now given by the source signal scaled by the 
% arra response plus the noise matrix
X = X + Noise;

