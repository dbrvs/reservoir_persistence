figure(1)
clf;

N=5000;

al=0

for i=1:5

    R=10^i;

    [r,a,cpa]=sampling(ones([R,1])/R,N);

    subplot(121)
    hold on
    loglog(r,a/sum(a),'-o')
    ylim([1e-2,1])
    set(gca,'Xscale','Log')
    set(gca,'Yscale','Log')
    title('uniform')
    xlabel('rank')
    ylabel('proportional abundance')
    
    pa=(1:R).^(-al);
    pa=pa/sum(pa);
    [r,a,cpa]=sampling(pa,100);

    subplot(122)
    hold on
    loglog(r,a/sum(a),'-o')
    set(gca,'Xscale','Log')
    set(gca,'Yscale','Log')
    title('power law \alpha=-1')
    ylim([1e-2,1])
    xlabel('rank')
    
end

legend('R=10^1','10^2','10^3','10^4','10^5')