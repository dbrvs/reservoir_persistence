%test the residuals of fitting to CDF

%make the theoretical "true" distribution
R=1e4; %richness
L0=1e7; %reservoir size, note all L0 are possible for a given alpha and R, see a few lines later
r=1:R; %ranks

%sample the data of real model
ns=100;
al=1.5;
a_r=r.^(-al); 
pa=a_r/sum(a_r);
[r,a,cpa]=sampling(pa,ns);

%try different models and see if it returns the right answer
for i=1:10
    al=i/5;
    a_r=r.^(-al); 
    pa=a_r/sum(a_r);
    
    [r1,a1,cpa1]=sampling(pa,ns);

end


[r1,a1,cpa1]=sampling(pa,ns);
[r2,a2,cpa2]=sampling(pa,ns);

%%
clf
subplot(131)
semilogx(r1,cpa1)
hold on
semilogx(r2,cpa2)

subplot(132)
scatter(r1,cpa1-cpa2)
set(gca,'XScale','log')

subplot(133)
scatter(r1,(cpa1-cpa2).^2)
set(gca,'XScale','log')
set(gca,'YScale','log')
ylim([1e-6,1])
legend(['rms=' num2str(sqrt(mean((cpa1-cpa2).^2)))])

