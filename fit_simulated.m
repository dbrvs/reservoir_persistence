%%
clf; clear
%example of finding best fit from simulated undersampled data

%make the theoretical "true" distribution
true_R=1e5; %richness
true_L=1e7; %reservoir size, note all L0 are possible for a given alpha and R, see a few lines later
true_r=1:true_R; %ranks
true_al=1.5; %the pwl parameter we want to estimate

true_a=true_r.^(-true_al); %do powerlaw
true_pa=true_a/sum(true_a); %normalize it
true_a=round(true_pa*true_L); %make it discrete
true_R=sum(true_a>0); %recalculate richness
true_L=sum(true_a); %recalculate the reservoir size with rounding

%true_a=true_a(true_a>0); %
%true_r=1:true_R;
%true_pa=true_pa(1:length(true_r));

%plot the true distribution
figure(1)
loglog(true_r,true_a,'Color',[0.8 0 1],'LineWidth',3)

xlabel('simulated rank')
ylabel('simulated abundance')
set(gca,'YTick',logspace(0,7,8))
ylim([1e0,1e7])
set(gca,'XTick',logspace(0,5,6))
xlim([1e0,1e5])
title(['true \alpha = ' num2str(true_al)...
    ', true R = ' num2str(true_R)...
    ', true L = ' num2str(true_L)])


%%
%sample the data
num_samples=1000;
sim_data=mnrnd(num_samples,true_pa,1);
sim_r=1:length(sim_data);
sim_pa=sim_data/sum(sim_data);
sim_cpa=cumsum(sim_pa); %cumulative proportional abundance

%now show what happens when you don't know about missing data
collapsed_data=-sort(-sim_data(sim_data>0)); %don't know which ones were missed
sim_R=length(collapsed_data);
collapsed_r=1:sim_R;
collapsed_pa=collapsed_data/sum(collapsed_data);
collapsed_cpa=cumsum(collapsed_pa); %cumulative proportional abundance

%plot the true distribution
figure(2)
clf
subplot(121)
loglog(sim_r,sim_data,'x','MarkerSize',10,'Color',[0.1 0.6 0])
xlabel('true rank')
ylabel('sampled abundance')
set(gca,'XTick',[1,1e1,1e2,1e3,1e4,1e5])
ylim([0.8,100])
xlim([0.8,1e3])

subplot(122)
loglog(collapsed_r,collapsed_data,'ko','MarkerSize',10)
xlabel('observed rank')
ylabel('observed abundance')
set(gca,'XScale','log')
set(gca,'YScale','log')
set(gca,'XTick',[1,1e1,1e2,1e3,1e4,1e5])
ylim([0.8,100])
xlim([0.8,1e3])
%set(gca,'YTick',[1e-4 1e-3 1e-2 1e-1 1])
%legend(['true \alpha = ' num2str(true_al)],'simulated','collapsed')


%% setup for real fit

% general model parameters
ins=1;%score index

num_al=50; num_R=50; 

num_fits=num_al*num_R;

score_mat=zeros([num_al,num_R]); %initialize score
models=zeros([num_al*num_R,3]); %initialize model list

al=linspace(0,5,num_al); %alpha range

%% the fitting loop
if num_al*num_R<500
    figure(2)
    clf
    hold on
end

tic
for i=1:num_al

    %calculate Chao estimate
    n1=length(sim_data(sim_data==1));
    n2=length(sim_data(sim_data==2));
    Rchao=log10(round(sim_R+n1*(n1-1)/(2*(n2+1))));

    %calculate approximate max richness
    psi=0.5*(1-true_L^(1-al(i)))/(true_L*(al(i)-1));
    maxR=min([7,log10(exp(log(psi)./-al(i)))]);%approximate log10 max richness
    R(i,:)=logspace(Rchao,maxR,num_R); %richness range
    
    for j=1:num_R
 
        r=1:R(i,j); %ranks
        f_r=r.^(-al(i)); %pwl1
        if num_al*num_R<500
            loglog(r,round(f_r/sum(f_r)*true_L))
        end
        mscore=calcscore(f_r,collapsed_pa,num_samples);
        score_mat(i,j)=mscore.avg;
        %score_var(i,j)=mscore.std;
        models(ins,:)=[mscore.avg al(i) R(i,j)];
        ins=ins+1;
    end
toc
end

if num_al*num_R<500
    xlabel('model rank')
    ylabel('model abundance')
    set(gca,'XTick',logspace(0,7,8),'XScale','log',...
            'YTick',logspace(0,7,8),'YScale','log')
    xlim([0.8,2e7])
    ylim([0.5,1e7])
end

%% plot best fits
figure(3)
clf
%R=logspace(1,6,num_R);
plotter(R',true_L,al,score_mat,models,num_samples,collapsed_r,collapsed_data,collapsed_cpa,'simulated \alpha=1.5');

            
%% compare various modeling techniques

naive_fit=polyfit(log10(collapsed_r),log10(collapsed_pa),1); %fit to sampled data

naive_R=1e4; %naive richness from Chao?
naive_r=1:naive_R;
naive_a=naive_r.^(naive_fit(1));
naive_pa=naive_a/sum(naive_a);

[CZMalpha, CZMxmin, CZML]=plfit(collapsed_data,'finite','xmin',1);
CZM_a=naive_r.^(-CZMalpha);
CZM_pa=CZM_a/sum(CZM_a);

my_pa=naive_r.^(-best_models);
my_pa=my_pa/sum(my_pa);

figure(5)
clf
hold on
loglog(true_r,true_pa,'Color',[0.8 0 1],'LineWidth',3)
loglog(collapsed_r,collapsed_pa,'ko','MarkerSize',10)
loglog(naive_r,naive_pa,'--r','LineWidth',1)
loglog(naive_r,CZM_pa,'--g','LineWidth',1)
loglog(naive_r,my_pa,'-k','LineWidth',1)
xlabel('rank')
ylabel('proportional abundance')
set(gca,'XScale','log')
set(gca,'YScale','log')
set(gca,'XTick',[1,1e1,1e2,1e3,1e4,1e5])
ylim([5e-3,1])
xlim([0.8,1e3])
set(gca,'YTick',[1e-4 1e-3 1e-2 1e-1 1])
legend(['true \alpha = ' num2str(true_al)],'observed samples',...
        ['naive \alpha = ' num2str(round(-naive_fit(1),2))],...
        ['CZM \alpha = ' num2str(round(CZMalpha,2))],...
        ['our model \alpha= ' num2str(round(best_models,2))])

%% compare model score and model variance (sensitivity to sampling)

figure(4)
clf
hold on
for i=1:num_R
    plot(score_mat(i,:),score_var(i,:))
    
end