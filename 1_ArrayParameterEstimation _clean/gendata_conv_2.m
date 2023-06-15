function x = gendata_conv_2(s,P,N,sigma)

    x = zeros([N 1]);
    
    for i=1:N
        t = i-1/P;
    
        x(i) = 0;
        for j=1:N
            x(i) = x(i) + h(t-(j-1))*s(j);
        end 
        x(i) = x(i);
    end
    
   % x = x + sigma^2*((randn([N 1])+ randn([N 1])*1i)./sqrt(2));
    x = x+ wgn(1, N, sigma^2, 'linear', 'complex');
end