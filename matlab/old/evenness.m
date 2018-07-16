%Evenness

clear
data_individual

figure(3); clf;
for ij=1:length(DNA);

    times=tDNA(ij);
    counts=DNA{ij};
    Robs=length(counts);

    pj=counts/sum(counts);
    S=-sum(pj.*log(pj)); %shannon entropy
    E=S/log(Robs); %eveness
     
    GS=1-sum(pj.^2); %Gini Simpson
    
    subplot(131)
    scatter(times,S,100,'o')
    hold on
    %set(gca,'YScale','log')
    xlabel('time (years on ART)')
    ylabel('Shannon entropy, S')
    %legend('obs','Chao1','95%CI')

    subplot(132)
    scatter(times,E,100,'o')
    hold on
    %set(gca,'YScale','log')
    xlabel('time (years on ART)')
    ylabel('evenness, E')
    %legend('obs','Chao1','95%CI')

    subplot(133)
    scatter(times,GS,100,'o')
    hold on
    %set(gca,'YScale','log')
    xlabel('time (years on ART)')
    ylabel('Gini-Simpson index, 1-\lambda')
    %legend('obs','Chao1','95%CI')

end

