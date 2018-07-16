
%%
clear
WB_1=[5 3 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];

data=-sort(-WB_1);
Lobs=sum(data);
Robs=length(data);


figure(1); clf;
hold on

al=linspace(0,5,50);
R=round(logspace(2,5,50));
ksmat=zeros([length(al),length(R)]);
ksp=zeros([length(al),length(R)]);

in1=1;
for i=al
    in2=1;
    for j=R
        r=1:j;
        clear p_r d
        p_r=r.^-i; %powerlaw
        true_pa=p_r/sum(p_r);

        replicates=100;
        pp=0;
        for k=1:replicates
            sim_data=mnrnd(Lobs,true_pa,1);
            collapsed_data=-sort(-sim_data(sim_data>0)); %don't know which ones were missed

            [h,p(k)]=kstest2(data,collapsed_data,'Alpha',0.01);
           
         
          
            if h==0 %they are from the same distribution
                pp=pp+1; %how many times out of replicates was the test significant
                %disp('same')
                ksmat(in1,in2)=pp;            
            end
        end
        
        ksp(in1,in2)=mean(p);

        %sim_data=mnrnd(Lobs,true_pa,1);
        %collapsed_data=-sort(-sim_data(sim_data>0)); %don't know which ones were missed
        %[h,ksp(in1,in2)]=kstest2(data,collapsed_data,'Alpha',0.01); %probability of same

        %semilogx(R,l1,'ko');
        %plot(R,l2,'rs');
        in2=in2+1;
    end
in1=in1+1;
disp(in1)
end

%%
 subplot(121)
 contourf(R,al,ksmat)
 set(gca,'XScale','log')
 colorbar
 
 subplot(122)
contourf(R,al,log10(ksp))
set(gca,'XScale','log')
colorbar
