function [MVDR_estimates] = MVDR_estimates(X,a)

    
    Rx = X.*X;
    MVDR_estimates = inv(Rx).*a/(a'.*Rx^(-1).*a)
