
function X = gendata_conv(s,P,N,sigma)

X = zeros([2*P N-1]);

for column=1:N-1
    for row=1:2*P
        t = column-1 + (row-1)/P;

        x = 0;
        for i=1:N
            x = x + h(t-(i-1))*s(i);
        end 
        X(row, column) = x;

    end
end

X = X + wgn(2*P, N-1, sigma^2, 'linear', 'complex');

