%% plot best fits
%DBR 2018

%inputs:
%R - richness
%L0 - reservoir size (number cells)
%al0 - powerlaw parameter
%smat - for contour plot, all scores in al x R matrix
%models - for sorting, list of all models [score al R] x number models
%ns - number of samples from real data
%dr - ranks from real data
%da - real abundance data
%dcpa - cumulative proportaional abundance frmo real data
%fn - file name to print figure

%outputs:
%best_al - the best powerlaw parameter

function best_al=plotter(R,L0,al,smat,models,ns,dr,da,dcpa,fn)
    
    %calculate Chao estimate
    n1=length(da(da==1));
    n2=length(da(da==2));
    Robs=length(da);
    Rchao=round(Robs+n1^2/(2*n2));

    % sort the models to get best fits
    [~,bf_idx] = sort(models(:,1)); %sort best fits
    num_bfs=5;
    
    %loop over num_bfs fits, and plot
    for i=1:num_bfs

        j=bf_idx(i); %use the sorted models
        r=1:models(j,3); %make the ranks up to R for the model

        bm=r.^(-models(j,2)); bm=bm/sum(bm); %compute the 'best' model

        num_replicates=10; %how many time to sample for best models
        maxd=10000; %max data length for replicate array filling
        rep_data=zeros([1,maxd]); %replicate data arary
        %loop over replicates
        for j=1:num_replicates
            mabund=mnrnd(ns,bm,1); %model abundance
            dm=mabund(mabund>0); %don't know which ones were missed
            rep_data(1:length(dm))=rep_data(1:length(dm))+dm; %replicate data
        end
        rep_data=rep_data/sum(rep_data); %normalize

        %sampled data fit
        subplot(221); hold on
        pl=plot(1:maxd,100*cumsum(rep_data),'LineWidth',1); pl.Color(4) = 0.7;

        %extrapolate cumulative abundance
        subplot(223); hold on
        m_cpa=cumsum(bm)*100;
        pl=plot(r,m_cpa,'LineWidth',1); pl.Color(4) = 0.7;

        %extrapolated rank abundance
        subplot(224); hold on
        aex=round(bm*L0);
        pl=loglog(r,aex,'LineWidth',1); pl.Color(4) = 0.7;
    end

    %sampled data and best fit cpas
    subplot(221)
    xlabel('sampled rank')
    ylabel({'sampled cumulative'; 'proportional abundance'})
    scatter(dr,100*dcpa,10,'k')
    xlim([1,dr(end)+1])
    ylim([-5,110])

    %richness and model param heatmap
    subplot(222)
    contourf(log10(R),al,log10(smat),30,'LineColor','none')
    hold on
    ylabel('model parameter, \alpha')
    xlabel('model richness, log_{10}R')
    colorbar;
    caxis([min(log10(models(:,1))),max(log10(models(:,1)))])
    %shade out below R chao
    if Rchao>min(models(:,3))
        %shade parts of plot
        Rshade=linspace(log10(R(1)),log10(Rchao),100);
        h=area(Rshade,ones([length(Rshade),1])*al(1),al(end));
        set(h(1),'FaceColor','k');
    end
    %shade out above R max
    R0=L0; %assume max R could be richness, i.e. flat
    alr=linspace(0,5,100);
    psi=0.5*(1-R0.^(1-alr))./(L0*(alr-1));
    mrshade=log10(exp(log(psi)./-alr));%approximate log10 max richness
    H=area(fliplr(mrshade),fliplr(alr),alr(end));
    set(H(1),'FaceColor','k');
    xlim([log10(R(1)),log10(R(end))])
    ylim([al(1),al(end)])
    %title(['best fit \alpha = ' num2str(best_al) ],'Fontsize',8,'FontWeight','Normal')

    %extrapolated cpa
    subplot(223)
    xlabel('extrapolated rank')
    ylabel({'extrapolated cumulative';'proportional abundance'})
    grid on
    ylim([0,100])
    set(gca,'XScale','log')
    set(gca,'XTick',logspace(1,7,7))
    xlim([1,1e7])
    ylim([-5,110])

    %extrapolated rank abundance
    subplot(224)
    xlabel('extrapolated rank')
    ylabel({'extrapolated abundance'})
    set(gca,'YScale','log')
    set(gca,'XScale','log')
    set(gca,'XTick',logspace(1,7,7))
    xlim([1,R(end)])
    ylim([1,L0])
    set(gca,'YTick',logspace(0,ceil(log10(L0)),ceil(log10(L0))+1))
    
    %print the figure
    w=7;
    h=6;
    u='inches';
    pp=0.01;

    set(gcf,'Units',u);
    screenpos = get(gcf,'Position');

    set(gcf,...
      'Position',[screenpos(1:2) w h],...
      'PaperUnits',u,...
      'PaperPosition',[pp*w pp*h w h],...
      'PaperSize',[w*(1+2*pp) h*(1+2*pp)]);

    %tightfig;
    print(fn,'-dpng','-r600') 

    %best_params=[models(bf_idx(1),2),log10(models(bf_idx(1),3))]; %the best fit

    %calculate best parameter
    best_al=round(mean(models(bf_idx(1:10),2)),2);

end