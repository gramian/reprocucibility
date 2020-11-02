
load('mathmod.mat');

if(N>100) n  = 100; else n = N; end

h = figure;

subplot(10,1,[1,7]);
semilogy(1:n,L2ERROR(1,1:n)./ORIGL2,'r','linewidth',2);
hold on;
semilogy(1:n,L8ERROR(1,1:n)./ORIGL8,'g--','linewidth',2);
semilogy(1:n,L2ERROR(2,1:n)./ORIGL2,'b','linewidth',2);
semilogy(1:n,L8ERROR(2,1:n)./ORIGL8,'m--','linewidth',2);
hold off;
legend('Balanced Truncation (L2)',...
       'Balanced Truncation (L\infty)',...
       'Cross Gramian (L2)',...
       'Cross Gramian (L\infty)');
xlabel('Reduced States');
ylabel('Relative Output Error');
xlim([1,n]);
ylim([1e-17,1]);

subplot(10,1,[9,10])
barh(fliplr(OFFLINE));
set(gca,'YTickLabel',{'WX';'BT'});
xlabel('Offline Time [s]');

print(h,'-depsc','mm.eps'); 
