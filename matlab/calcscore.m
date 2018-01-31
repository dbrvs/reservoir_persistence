function rms_score=fit_calcscore(f_r,data_pa,ns)

    pa=f_r/sum(f_r); %normalize the function

    num_replicates=5; %how many replicate simulations/fits

    rms_score=0; %initialize average for replicates
    
    for i=1:num_replicates
        
        sample_abund=mnrnd(ns,pa,1); %model abundance
        sample_collapse=sample_abund(sample_abund>0); %don't know which ones were missed

        sample_pa=-sort(-sample_collapse)/sum(sample_collapse); %make sure sorted, and normalize

        %make sure lengths are ok
        sub_mat=zeros(2,max(length(data_pa),length(sample_pa)));
        sub_mat(1,1:length(data_pa))=data_pa;
        sub_mat(2,1:length(sample_pa))=sample_pa;

        %calculate proportional cumulative abundances
        data_cpa=cumsum(sub_mat(1,:));
        sample_cpa=cumsum(sub_mat(2,:));

        rms_score=rms_score+sqrt(sum((data_cpa-sample_cpa).^2)); %rms score from cpas
        %mscore=max(abs(data_cpa-sample_cpa)); %KS score from cpas
    end
    
    rms_score=rms_score/num_replicates; %take average over replicates

end