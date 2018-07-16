%check zeta approximation!

al=linspace(1.001,2,15);
R=round(logspace(2,7,10));
for i=1:length(al)
    for j=1:length(R)
        r=1:R(j);
        
        p_r=r.^-al(i); %powerlaw
        S=sum(p_r);
        Z=zeta(al(i));
        I=(1-R(j)^(1-al(i)))/(al(i)-1);
        
        errI(i,j)=abs(S-I)/S;
        errZ(i,j)=abs(S-Z)/S;
        
    end
end

%%

subplot(221)
contourf(R,al,log10(errZ))
set(gca,'XScale','log')
colorbar

subplot(222)
plot(al,log10(errZ))

subplot(223)
contourf(R,al,log10(errI))
set(gca,'XScale','log')
colorbar

subplot(224)
plot(al,log10(errZ))
hold on
plot(al,log10(errI))

%%
figure(1); clf

semilogy(al,errZ*100,'Color',[0 1 0.9])
hold on
%legend('Zeta approx','Integral approx')
semilogy(al,errI*100,'Color',[0.9 0 1])
ylim([0.01,100])
xlabel('power law exponent, \alpha')
ylabel('percent error')
grid on

