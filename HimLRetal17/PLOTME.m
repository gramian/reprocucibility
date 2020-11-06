%%% summary: GAMM / PAMM 2017 companion code
%%% project: Fast Low-Rank Empirical Cross Gramians
%%% authors: Christian Himpe ( 0000-0003-2194-6754 )
%%% license: 2-Clause BSD (2017)
%$

load('runme.mat')

figure;
ttt = twx + tsv;
h = barh([fliplr(ttt);fliplr(tsv);fliplr(twx)]);
colormap(prism);
set(gca,'YTickLabel',{'Total','POD','WX'});
set([gca; findall(gca,'Type','text')],'FontSize',18);
l = legend('Distributed','Full','location','east');
set(l,'FontSize',18);
pbaspect([2,1,1]);
xlabel('[s]');
print('-dpdf','timing.pdf');

figure('Name',mfilename,'NumberTitle','off');
semilogy(1:m,l1(1:m),'r','linewidth',5); hold on;
semilogy(1:m,l2(1:m),'g--','linewidth',5);
semilogy(1:m,l8(1:m),'b','linewidth',5); hold off;
xlim([1,m]);
ylim([1e-9,1]);
xlabel('Reduced Order');
ylabel('Output Error');
pbaspect([1,0.5,0.5]);
set(gca,'YTick',10.^[-10:-2:0]);
set([gca; findall(gca,'Type','text')],'FontSize',18);
l = legend('L1 Error ','L2 Error ','L8 Error ','location','west');
set(l,'FontSize',18);
print('-dpdf','error.pdf');
