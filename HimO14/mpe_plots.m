
load('mpe.mat');


h = figure;
subplot(1,2,1);
semilogy(1:N-1,L2ERROR(1,1:N-1)/ORIGLL2,'r','LineWidth',2);
hold on;
semilogy(1:N-1,L2ERROR(2,1:N-1)/ORIGLL2,'k','LineWidth',2);
semilogy(1:N-1,L2ERROR(3,1:N-1)/ORIGLL2,'b','LineWidth',2);
hold off;
xlim([1 N-1]); ylim([1e-12 1]); 
legend('Balanced POD','Balanced Truncation','Direct Truncation');
set([gca; findall(gca, 'Type','text')], 'FontSize', 16);

subplot(1,2,2); 
semilogy(1:N-1,L2ERROR(10,1:N-1)/ORIGNL2,'r','LineWidth',2);
hold on;
semilogy(1:N-1,L2ERROR(11,1:N-1)/ORIGNL2,'k','LineWidth',2);
semilogy(1:N-1,L2ERROR(12,1:N-1)/ORIGNL2,'b','LineWidth',2);
hold off;
xlim([1 N-1]); ylim([1e-12 1]); 
legend('Balanced POD','Balanced Truncation','Direct Truncation'); 
set([gca; findall(gca, 'Type','text')], 'FontSize', 16);


h = figure;
subplot(1,2,1);
semilogy(1:N-1,L2ERROR(4,1:N-1)/ORIGLL2,'r','LineWidth',2);
hold on;
semilogy(1:N-1,L2ERROR(5,1:N-1)/ORIGLL2,'k','LineWidth',2);
semilogy(1:N-1,L2ERROR(6,1:N-1)/ORIGLL2,'b','LineWidth',2);
hold off;
xlim([1 N-1]); ylim([1e-7 1]); 
legend('Sensitivity Gramian','Identifiability Gramian','Cross-Identifiability Gramian');
set([gca; findall(gca, 'Type','text')], 'FontSize', 16);

subplot(1,2,2); 
semilogy(1:N-1,L2ERROR(13,1:N-1)/ORIGNL2,'r','LineWidth',2);
hold on;
semilogy(1:N-1,L2ERROR(14,1:N-1)/ORIGNL2,'k','LineWidth',2);
semilogy(1:N-1,L2ERROR(15,1:N-1)/ORIGNL2,'b','LineWidth',2);
hold off;
xlim([1 N-1]); ylim([1e-7 1]); 
legend('Sensitivity Gramian','Identifiability Gramian','Cross-Identifiability Gramian');
set([gca; findall(gca, 'Type','text')], 'FontSize', 16);


h = figure;
subplot(1,2,1);
semilogy(1:N-1,L2ERROR(7,1:N-1)/ORIGLL2,'r','LineWidth',2);
hold on;
semilogy(1:N-1,L2ERROR(8,1:N-1)/ORIGLL2,'k','LineWidth',2);
semilogy(1:N-1,L2ERROR(9,1:N-1)/ORIGLL2,'b','LineWidth',2);
hold off;
xlim([1 N-1]); ylim([1e-7 1]); 
legend('Controllability-Based','Observability-Based','Cross-Gramian-Based');
set([gca; findall(gca, 'Type','text')], 'FontSize', 16);

subplot(1,2,2); 
semilogy(1:N-1,L2ERROR(16,1:N-1)/ORIGNL2,'r','LineWidth',2);
hold on;
semilogy(1:N-1,L2ERROR(17,1:N-1)/ORIGNL2,'k','LineWidth',2);
semilogy(1:N-1,L2ERROR(18,1:N-1)/ORIGNL2,'b','LineWidth',2);
hold off;
xlim([1 N-1]); ylim([1e-7 1]); 
legend('Controllability-Based','Observability-Based','Cross-Gramian-Based');
set([gca; findall(gca, 'Type','text')], 'FontSize', 16);


h = figure;
subplot(6,1,1); hold on; barh(1,OFFLINE(3),'b'); barh(2,OFFLINE(2),'k'); barh(3,OFFLINE(1),'r'); hold off;
ylabel({'Linear','State','Reduction'}); xlim([0 60]);
set([gca; findall(gca, 'Type','text')], 'FontSize', 16); set(gca,'XTickLabel',{},'Ytick',[1,2,3],'YTicklabel',{'WX','WC+WO','WY'});

subplot(6,1,2); hold on; barh(1,OFFLINE(12),'b'); barh(2,OFFLINE(11),'k'); barh(3,OFFLINE(10),'r'); hold off;
ylabel({'Noninear','State','Reduction'}); xlim([0 60]);
set([gca; findall(gca, 'Type','text')], 'FontSize', 16); set(gca,'XTickLabel',{},'Ytick',[1,2,3],'YTicklabel',{'WX','WC+WO','WY'});

subplot(6,1,3); hold on; barh(1,OFFLINE(6),'b'); barh(2,OFFLINE(5),'k'); barh(3,OFFLINE(4),'r'); hold off;
ylabel({'Linear','Parameter','Reduction'}); xlim([0 60]);
set([gca; findall(gca, 'Type','text')], 'FontSize', 16); set(gca,'XTickLabel',{},'Ytick',[1,2,3],'YTicklabel',{'         WJ','WI','WS'});

subplot(6,1,4); hold on; barh(1,OFFLINE(15),'b'); barh(2,OFFLINE(14),'k'); barh(3,OFFLINE(13),'r'); hold off;
ylabel({'Nonlinear','Parameter','Reduction'}); xlim([0 60]);
set([gca; findall(gca, 'Type','text')], 'FontSize', 16); set(gca,'XTickLabel',{},'Ytick',[1,2,3],'YTicklabel',{'         WJ','WI','WS'});

subplot(6,1,5); hold on; barh(1,OFFLINE(9),'b'); barh(2,OFFLINE(8),'k'); barh(3,OFFLINE(7),'r'); hold off; 
ylabel({'Linear','Combined','Reduction'}); xlim([0 60]);
set([gca; findall(gca, 'Type','text')], 'FontSize', 16); set(gca,'XTickLabel',{},'Ytick',[1,2,3],'YTicklabel',{'         WJ',' WC+WI','WS+WO'});

subplot(6,1,6); hold on; barh(1,OFFLINE(18),'b'); barh(2,OFFLINE(17),'k'); barh(3,OFFLINE(16),'r'); hold off;
ylabel({'Nonlinear','Combined','Reduction'}); xlim([0 60]); xlabel('t[s]'); 
set([gca; findall(gca, 'Type','text')], 'FontSize', 16); set(gca,'Ytick',[1,2,3],'YTicklabel',{'         WJ','WC+WI','WS+WO'});


h = figure;
subplot(1,2,1);
loglog(1:M-1,L2ERROR(19,1:M-1)/ORIGBL2,'r','LineWidth',2);
hold on;
loglog(1:M-1,L2ERROR(20,1:M-1)/ORIGBL2,'b','LineWidth',2);
hold off;
xlim([1 M-1]); ylim([1e-16 1]); 
legend('Balanced Truncation','Direct Truncation');
set([gca; findall(gca, 'Type','text')], 'FontSize', 16);


cmap = autumn(100); cmap = flipud(cmap); tmp = linspace(0,1,100)';
for n = 1:3, cmap(:,n) = interp1( 10.^tmp, cmap(:,n), 1+9*tmp, 'linear'); end

ticks = [1,N/2,N]; 

h = figure;
subplot(1,3,1);
h = surf(l2error(:,:,7)./ORIGLL2); shading interp; view(150,30); set(h,'EdgeColor',[0 0 0]);
set(gca,'xtick',linspace(0,sqrt(N),3),'ytick',linspace(0,sqrt(N),3));
colormap(cmap); set(gca,'zscale','log','XTickLabelMode','manual','YTickLabelMode','manual','Xticklabel',ticks,'Yticklabel',ticks);
set(gca,'ZMinorGrid','off');
xlim([1 sqrt(N)]); ylim([1 sqrt(N)]); zlim([1e-7 1]);
xlabel('Parameters'); ylabel('States'); title('WS+WO');
%set([gca; findall(gca, 'Type','text')], 'FontSize', 16);

subplot(1,3,2);
h = surf(l2error(:,:,8)./ORIGLL2); shading interp; view(150,30); set(h,'EdgeColor',[0 0 0]);
set(gca,'xtick',linspace(0,sqrt(N),3),'ytick',linspace(0,sqrt(N),3));
colormap(cmap); set(gca,'zscale','log','XTickLabelMode','manual','YTickLabelMode','manual','Xticklabel',ticks,'Yticklabel',ticks);
set(gca,'ZMinorGrid','off');
xlim([1 sqrt(N)]); ylim([1 sqrt(N)]); zlim([1e-7 1]);
xlabel('Parameters'); ylabel('States'); title('WC+WI');
%set([gca; findall(gca, 'Type','text')], 'FontSize', 16);

subplot(1,3,3);
h = surf(l2error(:,:,9)./ORIGLL2,'EdgeColor','k'); shading interp; view(150,30); set(h,'EdgeColor',[0 0 0]);
set(gca,'xtick',linspace(0,sqrt(N),3),'ytick',linspace(0,sqrt(N),3));
colormap(cmap); set(gca,'zscale','log','XTickLabelMode','manual','YTickLabelMode','manual','Xticklabel',ticks,'Yticklabel',ticks);
set(gca,'ZMinorGrid','off');
xlim([1 sqrt(N)]); ylim([1 sqrt(N)]); zlim([1e-7 1]);
xlabel('Parameters'); ylabel('States'); title('WJ');
%set([gca; findall(gca, 'Type','text')], 'FontSize', 16);


h = figure;
subplot(1,3,1);
h = surf(l2error(:,:,16)./ORIGNL2,'EdgeColor','k'); shading interp; view(150,30); set(h,'EdgeColor',[0 0 0]);
set(gca,'xtick',linspace(0,sqrt(N),3),'ytick',linspace(0,sqrt(N),3));
colormap(cmap); set(gca,'zscale','log','XTickLabelMode','manual','YTickLabelMode','manual','Xticklabel',ticks,'Yticklabel',ticks);
set(gca,'ZMinorGrid','off');
xlim([1 sqrt(N)]); ylim([1 sqrt(N)]); zlim([1e-7 1]);
xlabel('Parameters'); ylabel('States'); title('WS+WO');
%set([gca; findall(gca, 'Type','text')], 'FontSize', 16);

subplot(1,3,2);
h = surf(l2error(:,:,17)./ORIGNL2,'EdgeColor','k'); shading interp; view(150,30); set(h,'EdgeColor',[0 0 0]);
set(gca,'xtick',linspace(0,sqrt(N),3),'ytick',linspace(0,sqrt(N),3));
colormap(cmap); set(gca,'zscale','log','XTickLabelMode','manual','YTickLabelMode','manual','Xticklabel',ticks,'Yticklabel',ticks);
set(gca,'ZMinorGrid','off');
xlim([1 sqrt(N)]); ylim([1 sqrt(N)]); zlim([1e-7 1]);
xlabel('Parameters'); ylabel('States'); title('WC+WI');
%set([gca; findall(gca, 'Type','text')], 'FontSize', 16);

subplot(1,3,3);
h = surf(l2error(:,:,18)./ORIGNL2,'EdgeColor','k'); shading interp; view(150,30); set(h,'EdgeColor',[0 0 0]);
set(gca,'xtick',linspace(0,sqrt(N),3),'ytick',linspace(0,sqrt(N),3));
colormap(cmap); set(gca,'zscale','log','XTickLabelMode','manual','YTickLabelMode','manual','Xticklabel',ticks,'Yticklabel',ticks);
set(gca,'ZMinorGrid','off');
xlim([1 sqrt(N)]); ylim([1 sqrt(N)]); zlim([1e-7 1]);
xlabel('Parameters'); ylabel('States'); title('WJ');
%set([gca; findall(gca, 'Type','text')], 'FontSize', 16);

