function [theta,f] = joint(X,d,m)

% Obtain data parameters
M = size(X, 1);
N = size(X, 2);

% Apply smoothing in time with a factor m
for column_id = 1: N-m+1
    new_column = zeros([m*M 1]);
    i =1;

    for sample_index = column_id:column_id+m-1
        new_column((i-1)*M+1: i*M) = X(:, sample_index);
        i = i+1;
    end

    Z(:, column_id) = new_column;
end

% Compute trucated SVD of Z
[Ux,Sx,Vx] = svd(Z,"econ");

Ux = Ux(:,1:d);
Sx = Sx(1:d,1:d);
Vx = Vx(:,1:d);

% Creating selection matrices which make use of the shift invariance in
% both F and A, phi and theta respectively
Jxphi = kron([eye(m-1),zeros(m-1,1)],eye(M));
Jyphi = kron([zeros(m-1,1),eye(m-1)],eye(M));

Jxtheta = kron(eye(m),[eye(M-1),zeros(M-1,1)]);
Jytheta = kron(eye(m),[zeros(M-1,1),eye(M-1)]);

% These selection matrices allow for estimation of Phi, by taking the
% submatrices which are the first and last M*(m-1) rows of Ux, and for
% estimation of Theta the first and last M-1 rows of all m blocks are
% stacked.

Uxphi = Jxphi*Ux;
Uyphi = Jyphi*Ux;

Uphipseud = inv(Uxphi'*Uxphi)*Uxphi';

TphiT = Uphipseud*Uyphi;

Uxtheta = Jxtheta*Ux;
Uytheta = Jytheta*Ux;

Uthetapseud = inv(Uxtheta'*Uxtheta)*Uxtheta';

TthetaT = Uthetapseud*Uytheta;

% Perform joint diagonalization, see section 9.1.3 of reader
% first d x d matrix in estimates is related to the frequencies
% second d x d matrix in estimates is related to the angle
[~, estimates] = joint_diag([TphiT,TthetaT],1e-8);

phiEst = diag(estimates(:,1:2));
thetaEst = diag(estimates(:,3:4));

f = sort(angle(phiEst)/(2*pi));
theta = sort(asind(angle(thetaEst)/pi));




















% % Creating selection matrices which make use of the shift invariance in
% % both F and A, phi and theta respectively
% % Also see page 178 of reader, there the same is done but for joint
% % angle delay estimation
% 
% Jxphi = kron([eye(N-1),zeros(N-1,1)],eye(M));
% Jyphi = kron([zeros(N-1,1),eye(N-1)],eye(M));
% 
% Jxtheta = kron(eye(N),[eye(M-1),zeros(M-1,1)]);
% Jytheta = kron(eye(N),[zeros(M-1,1),eye(M-1)]);
% 
% % These selection matrices allow for estimation of Phi, by taking the
% % submatrices which are the first and last M*(m-1) rows of Ux, and for
% % estimation of Theta the first and last M-1 rows of all m blocks are
% % stacked.
% % Also see page 178 of reader (eq 9.12, 9.13), there the same is done but for joint
% % angle delay estimation
% size(Jxphi)
% size(Ux)
% Uxphi = Jxphi*Ux;
% Uyphi = Jyphi*Ux;
% 
% Uphipseud = inv(Uxphi'*Uxphi)*Uxphi';
% 
% TphiT = Uphipseud*Uyphi;
% 
% Uxtheta = Jxtheta*Ux;
% Uytheta = Jytheta*Ux;
% 
% Uthetapseud = inv(Uxtheta'*Uxtheta)*Uxtheta';
% 
% TthetaT = Uthetapseud*Uytheta;
% 
% % Perform joint diagonalization, see section 9.1.3 of reader
% % first d x d matrix in estimates is related to the frequencies
% % second d x d matrix in estimates is related to the angle
% [~, estimates] = joint_diag([TphiT,TthetaT],1e-8);

% phiEst = diag(estimates(:,1:2));
% thetaEst = diag(estimates(:,3:4));
% 
% P = 100;
% 
% f = angle(phiEst)/(2*pi)*P;
% theta = asind(angle(thetaEst)/pi);
