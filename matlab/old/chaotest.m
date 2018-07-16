%Chao test

clear
data_individual

figure(1); clf;
for ij=1:length(DNA);

    times=tDNA(ij);
    counts=DNA{ij};
    Robs=length(counts);

    f1=sum(counts==1);
    f2=sum(counts==2);
    n=max(counts);

    nc=(n-1)/n;

    RChao1=Robs+nc*f1*(f1-1)/(2*(f2+1));

    fc=f1/f2;

    varRchao1=f2*(nc^2/4*fc^4 + nc^2*fc^3 + nc/2*fc^2);

    CIc=exp(1.96*sqrt(1+varRchao1/(RChao1-Robs)^2));

    CI=[Robs+(RChao1-Robs)/CIc;Robs+(RChao1-Robs)*CIc];

    f3=sum(counts==3);
    f4=sum(counts==4)+1;

    RiChao1=RChao1+(n-3)/(4*n)*f3/f4*max([0,f1-(n-3)/(2*(n-1))*f2*f3/4]);
    
    scatter(times,Robs,'o')
    hold on
    scatter(times,RChao1,'x')
    plot([times,times],CI)
    %scatter(times,RiChao1,'D')
    set(gca,'YScale','log')
    xlabel('time (years on ART)')
    ylabel('richness, R')
    legend('obs','Chao1','95%CI')

end

