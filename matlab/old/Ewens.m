%% EWENS sampling
data_individual

figure(2); clf;
for ij=1:length(DNA);

counts=DNA{ij};

k=max(counts);

%make clone size abundnaces
f=0;
for j=1:k;
    f(j)=sum(counts==j);
end

th=0.01:0.01:0.99;
for i=1:length(th);
    
thprod=prod((th(i)+0:(k-1)));

%probEw = factorial(k)/thprod*prod(th.^f./(factorial(f).*(1:k).^f))

%FACf=sqrt(2*pi*f).*(f/exp(1)).^f;

%do it all in logs using Stirlings for large factorials

lFf=log(factorial(f));
lFf(lFf>100)=f(lFf>100).*log(f(lFf>100))-f(lFf>100);
lEw(i) = log(factorial(k)) - log(thprod) + ...
        sum(f.*log(th(i)) - lFf + f.*log((1:k)));
    
end

hold on
plot(th,lEw/max(lEw))
xlabel('Ewens \theta')
ylabel('normalized log probability, $\log(p(f;\theta))/\max{\log(p(f;\theta))}$','Interpreter','Latex')

end
