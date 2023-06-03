% function x = gendata_conv(s,P,N,sigma)
% 
% x = zeros([N 1]);
% 
% for i=1:N
%    for k=1:N
%             x = x + h(t-(i-1))*s(i);
%         end 
% 
% 
% for column=1:N-1
%     for row=1:2*P
%         t = column-1 + (row-1)/P;
% 
%         x = 0;
%         for i=1:N
%             x = x + h(t-(i-1))*s(i);
%         end
%         x = x + sigma^2*((randn()+ randn()*1i)./sqrt(2));
%         X(row, column) = x;
% 
%     end
% end



function X = gendata_conv(s,P,N,sigma)

X = zeros([2*P N-1]);

for column=1:N-1
    for row=1:2*P
        t = column-1 + (row-1)/P;

        x = 0;
        for i=1:N
            x = x + h(t-(i-1))*s(i);
        end 
        x = x + sigma^2*((randn()+ randn()*1i)./sqrt(2));
        X(row, column) = x;

    end
end

