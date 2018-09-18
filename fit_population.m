%script to fit population level rank abundances
%DBR 2018

%The rank abundances were created by combining all counts on population level
%and then transforming to rank abundance

clc; clear;

data_individual() %loads up the data

%%
dflag=1; %0=DNA 1=Hosmane

if dflag==0
    %DNA Reservoir parameters
    fn='total-DNA';
    true_L0=1e9;
    tot_ra=-sort(-[WB_1,WB_6,WB_12,WL_1,WL_4,WL_11,WR_1,...
     WR_8,WR_12,M1_0,M1_4,M1_11,M2_12,M3_0,M3_7,M4_12,M5_14]);
elseif dflag==1
    %Replication competent Reservoir parameters
    fn='total-repcomp';
    true_L0=1e7;
    tot_ra=-sort(-[S01,S02,S03,S04,S05,S06,S07,S08,S09,S10,S11,S12]);
end

%calculate a few things with data
data_pa  = tot_ra/sum(tot_ra);
data_cpa = cumsum(data_pa);
data_r   = 1:length(tot_ra);
num_samples=sum(tot_ra);

num_al=50; num_R=50; num_fits=num_al*num_R; % general model parameters
al=linspace(0,3,num_al); %alpha range
R=logspace(3,6,num_R); %richness range

%% the fitting loop
score_mat=zeros([num_al,num_R]); %initialize score
models=zeros([num_al*num_R,3]); %initialize model list

tic
ins=1;%score index
%loop over power law parameters
for i=1:num_al
    %loop over possible R
    for j=1:num_R
        r=1:R(j); %ranks
        f_r=r.^(-al(i)); %pwl1
        mscore=calcscore(f_r,data_pa,num_samples);
        score_mat(i,j)=mscore.avg;
        %score_var(i,j)=mscore.std;
        models(ins,:)=[mscore.avg al(i) R(j)];
        ins=ins+1;
    end
end
toc


%% plot best fits

plotter(R,true_L0,al,score_mat,models,...
                            num_samples,data_r,tot_ra,data_cpa,fn);

