% logseries comparison to power law?
%%





clf

r=1:1000;

%log-series
th=0.9;
ls = -th.^r./r/log(1-th);

%power-law
al=1.1;
pl = r.^(-al);
pl = pl/sum(pl);

%log-normal
mu=2; sig=5;
ln = lognpdf(r,mu,sig);

%true distributions
figure(1);
semilogx(r,ls)
hold on
plot(r,pl)
plot(r,ln)
legend('log-series','power-law','lognormal')
ylim([-0.1,1])
xlabel('rank')
ylabel('proportional abundance')
%%
%simulated distributions
num_samples=100;
sim_ls=mnrnd(num_samples,ls,1);
sim_pl=mnrnd(num_samples,pl,1);

collapsed_pla=-sort(-sim_pl(sim_pl>0)); %don't know which ones were missed
collapsed_plr=1:length(collapsed_pla);
collapsed_pla=collapsed_pla/sum(collapsed_pla);

collapsed_lsa=-sort(-sim_ls(sim_ls>0)); %don't know which ones were missed
collapsed_lsr=1:length(collapsed_lsa);
collapsed_lsa=collapsed_lsa/sum(collapsed_lsa);

subplot(122)
loglog(collapsed_lsr,collapsed_lsa)
hold on
plot(collapsed_plr,collapsed_pla)

%% KS test

[h,p] = kstest2(collapsed_lsa,collapsed_lsr,'Alpha',0.05)



%% KS tests for uniformity or not?

clear

%check a totally uniform distribution
R=20;
for i=1:2
    ao=mnrnd(R,ones(R,1)/R,1);
    ao=-sort(-ao(ao>0));
    acdf{i}=cumsum(ao/sum(ao));
end

%h = kstest2(x1,x2) returns a test decision for the null hypothesis that 
%the data in vectors x1 and x2 are from the same continuous 
%distribution, using the two-sample Kolmogorov-Smirnov test. 
%The alternative hypothesis is that x1 and x2 are from 
%different continuous distributions. The result h is 1 
%if the test rejects the null hypothesis at the specified level,
%and 0 otherwise.

%ask, are they different?
%if different h=1 if same h=0
[h,p]=kstest2(acdf{1},acdf{2},'Alpha',0.05)


%% now use real data
clear
data_individual() %loads up the data
data_array={WB_1,WB_7,WB_12,WL_1,WL_4,WL_12,WR_2,...
            WR_8,WR_12,M1_0,M1_5,M1_11,M2_13,M3_0,M3_7,M4_12,M5_14};

fn={'WB-1','WB-7','WB-12','WL-1','WL-4','WL-12','WR-2',...
    'WR-8','WR-12','M1-0','M1-5','M1-11','M2-13','M3-0',...
    'M3-7','M4-12','M5-14'};

R=1000;
p_r=ones(R,1)/R; %uniform reservoir

for i=1:length(data_array)
    Robs(i)=length(data_array{i});
    Nobs(i)=sum(data_array{i});
    ad=data_array{i};
    
    %random sample 100 times
    for j=1:100
        ao=mnrnd(Nobs(i),p_r,1);
        ao=ao(ao>0);

        [ksh(j),ksp(j)] = kstest2(ad,ao,'Alpha',0.05);
    end

    KSH(i)=sum(ksh);
 end


figure(1); clf

subplot(121)
hold on
for i=1:length(data_array)
    if KSH(i)<50;
        mrk='x';
    else
        mrk='s';
    end
    scatter(Robs(i),KSH(i),mrk)
    %scatter(Nd(i),log10(ksp(i)),mrk)
end
xlabel('observed richness, R^{obs}')
set(gca,'XScale','log')
ylabel('proportion of sims with significant KS test (%)')
legend(fn)

%%
data_array={S01,S02,S03,S04,S05,S06,S07,S08,S09,S10,S11,S12}; %only good ones with N>20 [2,5,8,9,10]
fn={'S01','S02','S03','S04','S05','S06',...
    'S07','S08','S09','S10','S11','S12'};

for i=1:length(data_array)
    Robs(i)=sum(data_array{i});
    ad=data_array{i};
    
    %random sample 100 times
    for j=1:100
        ao=mnrnd(Robs(i),ones(R,1)/R,1);
        ao=ao(ao>0);

        [ksh(j),ksp(j)] = kstest2(ad,ao,'Alpha',0.05);
    end

    KSH(i)=sum(ksh);
end

subplot(122)
hold on
for i=1:length(data_array)
    if KSH(i)<50;
        mrk='x';
    else
        mrk='s';
    end
    scatter(Robs(i),KSH(i),mrk)
    %scatter(Nd(i),log10(ksp(i)),mrk)
end
set(gca,'XScale','log')
xlabel('Sample size')
ylabel('proportion of significant (p<0.05) KS tests (%)')
legend(fn)