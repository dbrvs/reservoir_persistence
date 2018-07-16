clear;

%data from Hosmane paper
S=[[3,0,0,0,0,0,0,0,1,0,0];...
[0,0,0,0,0,1,0,0,0,0,1];...
[8,3,3,0,0,0,0,1,0,0,0];...
[10,0,0,0,0,0,0,0,0,0,0];...
[4,3,0,1,0,0,0,0,0,0,0];...
[19,1,2,0,0,0,0,0,0,0,0];...
[5,0,0,0,0,0,0,0,0,0,0];...
[4,0,0,0,0,0,0,0,0,0,0];...
[8,1,0,1,0,0,1,0,1,0,0];...
[8,3,0,0,0,1,0,0,0,0,0];...
[11,0,0,0,0,0,0,1,0,0,0];...
[4,0,0,1,0,0,0,0,0,0,0]];

pop_S_nN=sum(S);

pop_S_ra=[];
%make rank abundance
for i=1:length(pop_S_nN)
    pop_S_ra=[pop_S_ra, ones([1,pop_S_nN(i)])*i];    
end
pop_S_ra=fliplr(pop_S_ra); %resort big->small


% data from all DNA
d1=[52  4  1  0  1];
d2=[44  3  2  1  1];
d3=[35  6  1  1];
d4=[32  1  0  0  0  0  2];
d5=[29  1  1  0  0  0  0  0  0  0  0  0  0  0  0  1];
d6=[29  2  2  0  0  0  0  0  1];
d7=[75  0  1];
d8=[53  2  0  0  1];
d9=[46  2  1  2];
d10=[196  13   0   2   1   0   2   0   1   0   0   0   0   0   1];
d11=[49  4  0  1  0  0  1];
d12=[590  49  13   7   1   5   3   0   1   0   0   0   0   0   1   0   0   0 0   0   0   0   0   0   1   1   0   0   0   0   0   0   0   0   0   0 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0 0   0   0   0   0   0   0   1];
d13=[30  4];
d14=[40  4  0  0  0  0  1];
d15=[168  10   0   2   1   0   1   0   0   0   0   0   0   0   1];
d16=[140   9   0   2   1   0   1   0   0   0   0   0   0   0   1];
d17=[93  7  0  2  1  0  1  0  0  0  0  0  0  0  1];

dtot={d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,d12,d13,d14,d15,d16,d17};

dmat=zeros([17,length(d12)]);
for i=1:length(dtot)
    d=dtot{i};
    dmat(i,1:length(d))=d;
end

pop_D_nN=sum(dmat);

pop_D_ra=[];
for i=1:length(pop_D_nN)
    pop_D_ra=[pop_D_ra, ones([1,pop_D_nN(i)])*i];
end

pop_D_ra=fliplr(pop_D_ra); %resort big->small

% plot population rank abundances
subplot(121)
stairs(1:length(pop_D_ra),pop_D_ra,'-ko','MarkerSize',3)
set(gca,'XScale','log')
xlabel('population rank (DNA)')
ylabel('population abundance (DNA)')
ylim([0,65])
xlim([0.8,3e3])
xticks([1,10,100,1e3])

subplot(122)
plot(1:length(pop_S_ra),pop_S_ra,'-ro','MarkerSize',3)
set(gca,'XScale','log')
xlabel('population rank (Hosmane)')
ylabel('population abundance (Hosmane)')
ylim([0,12])
xlim([0.8,200])
xticks([1,10,100])

    %print the figure
    w=6;
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

    print('population_ra','-dpng','-r600') 

    
%% check for normality

for i=1:length(S)
    [h(i),p(i)]=kstest(S{i})
