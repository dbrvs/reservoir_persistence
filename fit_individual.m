% data to fit

clc; clear all; close all

tic
data_individual() %loads up the data

pflag=1; %1 is to make plots
num_al=50; num_R=50; %number of models

for dflag=0 %0=DNA 1=Hosmane

if dflag==0
    %DNA Reservoir parameters
    L0=1e9;
    data_array=DNA;
    fn={'WB-1','WB-7','WB-12','WL-1','WL-4','WL-12','WR-2',...
        'WR-8','WR-12','M1-0','M1-5','M1-11','M2-13','M3-0',...
        'M3-7','M4-12','M5-14'};
elseif dflag==1
    %Replication competent Reservoir parameters
    L0=1e7;
    data_array=IU; %only good ones with N>20 [2,5,8,9,10]
    fn={'S03','S06','S09','S10','S11'};
end

best_params=zeros(length(data_array),2);

%% loop over data sets
for dd=9 %pick a single one
%parfor dd=1:length(data_array)

    data=-sort(-data_array{dd}); %make sure correctly ranked

    %calculate a few things with data
    data_pa  = data/sum(data);
    data_cpa = cumsum(data_pa);
    data_r   = 1:length(data);
    num_samples = sum(data); %number of experimental samplings

    % general model parameters
    num_fits=num_al*num_R;
    R=logspace(3,7,num_R); %richness range
    al=linspace(0.1,2,num_al); %alpha range
    
    %fitting loop
    ins=1; %score index
    score_mat=zeros([num_al,num_R]); %initialize score
    models=zeros([num_fits,3]); %initialize model list
    tic
    for i=1:num_al

        %calculate Chao estimate
        %n1=sum(data==1);
        %n2=sum(data==2);
        %Robs=length(data);
        %Rchao=log10(Robs+n1*(n1-1)/(2*(n2+1)));

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
    %disp(fn{dd})

    %% plot best fits
    if pflag==1
        figure()
        clf
        plotter(R,L0,al,score_mat,models,...
                                num_samples,data_r,data,data_cpa,fn{dd});
    end    

    sortedmodels=sort(models);
    best_params(dd,:)=sortedmodels(1,2:3);
    

end

%disp(dflag)
%disp(best_params)

end

toc
