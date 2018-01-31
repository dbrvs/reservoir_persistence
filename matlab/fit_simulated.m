%%
clf; clear
%example of finding best fit from simulated undersampled data

%make the theoretical "true" distribution
true_R=1e4; %richness
true_L0=1e7; %reservoir size, note all L0 are possible for a given alpha and R, see a few lines later
true_r=1:true_R; %ranks
true_al=1.5; %the pwl parameter we want to estimate
true_a=true_r.^(-true_al); 
true_pa=true_a/sum(true_a);
true_a=round(true_pa*true_L0);
true_L0=sum(true_a); %recalculate the reservoir size with rounding

%sample the data
num_samples=100;
sim_data=mnrnd(num_samples,true_pa,1);
sim_r=1:length(sim_data);
sim_pa=sim_data/sum(sim_data);
sim_cpa=cumsum(sim_pa); %cumulative proportional abundance

%now show what happens when you don't know about missing data
collapsed_data=-sort(-sim_data(sim_data>0)); %don't know which ones were missed
collapsed_r=1:length(collapsed_data);
collapsed_pa=collapsed_data/sum(collapsed_data);
collapsed_cpa=cumsum(collapsed_pa); %cumulative proportional abundance

%% plot true data and do the logloglinear pa fit and plot

figure(1); clf
subplot(121)
hold on
loglog(true_r,true_pa,'-g','LineWidth',3)
plot(sim_r,sim_pa,'ro')
plot(collapsed_r,collapsed_pa,'ko','MarkerSize',5)
hold off
xlabel('rank')
ylabel('proportional abundance')
set(gca,'XScale','log')
set(gca,'YScale','log')
set(gca,'XTick',[1,1e1,1e2,1e3,1e4,1e5])
ylim([1e-4,1])
set(gca,'YTick',[1e-4 1e-3 1e-2 1e-1 1])
legend(['true \alpha = ' num2str(true_al)],'simulated','collapsed')

naive_fit=polyfit(log10(collapsed_r),log10(collapsed_pa),1); %fit to sampled data
naive_R=1e3; %naive richness from Chao?
naive_r=1:naive_R;
naive_a=naive_r.^(naive_fit(1));
naive_pa=naive_a/sum(naive_a);

subplot(122)
hold on
loglog(true_r,true_pa,'-g','LineWidth',3)
loglog(naive_r,naive_pa,'-k','LineWidth',3)
loglog(collapsed_r,collapsed_pa,'ko','MarkerSize',5)
xlabel('rank')
ylabel('proportional abundance')
set(gca,'XScale','log')
set(gca,'YScale','log')
set(gca,'XTick',[1,1e1,1e2,1e3,1e4,1e5])
ylim([1e-4,1])
set(gca,'YTick',[1e-4 1e-3 1e-2 1e-1 1])
legend(['true \alpha = ' num2str(true_al)],['naive \alpha = ' num2str(round(-naive_fit(1),2))])

%print the figure
w=6;
h=3;
u='inches';
p = 0.01;

set(gcf,'Units',u);
screenpos = get(gcf,'Position');

set(gcf,...
  'Position',[screenpos(1:2) w h],...
  'PaperUnits',u,...
  'PaperPosition',[p*w p*h w h],...
  'PaperSize',[w*(1+2*p) h*(1+2*p)]);

print -dpng example_naive.png -r600

%% setup for real fit

% general model parameters
ins=1;%score index

num_param=50; num_R=50; 

num_fits=num_param*num_R;

score_mat=zeros([num_param,num_R]); %initialize score
models=zeros([num_param*num_R,3]); %initialize model list

R=logspace(1,5,num_R); %richness range

al=linspace(0.5,3,num_param); %alpha range

%% the fitting loop

figure(2)
clf
hold on

tic
for j=1:num_R
    for i=1:num_param
        r=1:R(j); %ranks
        f_r=r.^(-al(i)); %pwl1
        loglog(r,round(f_r/sum(f_r)*true_L0))
        mscore=calcscore(f_r,collapsed_pa,num_samples);
        score_mat(i,j)=mscore;
        models(ins,:)=[mscore al(i) R(j)];
        ins=ins+1;
    end
end
toc
%%
xlabel('model rank')
ylabel('model abundance')
set(gca,'XTick',logspace(0,5,6))
set(gca,'YTick',logspace(0,7,8))
xlim([1,2e5])
ylim([1,1e7])
    %print the figure
    w=3;
    h=3;
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
    print('example_models','-dpng','-r600') 

%% plot best fits
figure(3)
clf
filename='example_bestfit_shaded';
best_models=plotter(R,true_L0,al,score_mat,models,...
                            num_samples,collapsed_r,collapsed_data,collapsed_cpa,filename);

            