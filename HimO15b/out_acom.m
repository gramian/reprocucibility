% out_acom (Version 1.2)
% by Christian Himpe, 2015 ( http://wwwmath.uni-muenster.de/u/himpe )
% released under BSD 2-Clause License ( opensource.org/licenses/BSD-2-Clause )
%*

load('acom.mat');

m = 2;

figure;
subplot(1,3,[1]);
%z = squeeze(Z(1,2,m:M)); semilogy(m:M,z,'k','Linewidth',2); hold on;
z = squeeze(Z(2,1,m:M)); semilogy(m:M,z,'r','Linewidth',2); hold on;
z = squeeze(Z(3,1,m:M)); semilogy(m:M,z,'g','Linewidth',2);
z = squeeze(Z(4,1,m:M)); semilogy(m:M,z,'b','Linewidth',2);
z = squeeze(Z(5,1,m:M)); semilogy(m:M,z,'m','Linewidth',2); hold off;
xlim([m,M]);
ylim([1e-1,1e3]);
xlabel('Iteration')
ylabel('Offline Time [s]');
set(gca,'Xtick',[m:M]);
set(gca,'XtickLabel',[m-1:M-1]);

subplot(1,3,[2]);
z = squeeze(Z(1,2,m:M)); semilogy(m:M,z,'k','Linewidth',2); hold on;
z = squeeze(Z(2,2,m:M)); semilogy(m:M,z,'r','Linewidth',2);
z = squeeze(Z(3,2,m:M)); semilogy(m:M,z,'g','Linewidth',2);
z = squeeze(Z(4,2,m:M)); semilogy(m:M,z,'b','Linewidth',2);
z = squeeze(Z(5,2,m:M)); semilogy(m:M,z,'m','Linewidth',2); hold off;
xlim([m,M]);
ylim([1e-1,1e3]);
xlabel('Iteration')
ylabel('Online Time [s]');
set(gca,'Xtick',[m:M]);
set(gca,'XtickLabel',[m-1:M-1]);
title('Generic Forward Model');

subplot(1,3,[3]);
z = squeeze(Z(1,2,m:M)); semilogy(m:M,z,'k','Linewidth',2); hold on;
z = squeeze(Z(2,1,m:M)+Z(2,2,m:M)); semilogy(m:M,z,'r','Linewidth',2);
z = squeeze(Z(3,1,m:M)+Z(3,2,m:M)); semilogy(m:M,z,'g','Linewidth',2);
z = squeeze(Z(4,1,m:M)+Z(4,2,m:M)); semilogy(m:M,z,'b','Linewidth',2);
z = squeeze(Z(5,1,m:M)+Z(5,2,m:M)); semilogy(m:M,z,'m','Linewidth',2); hold off;
xlim([m,M]);
ylim([1e-1,1e3]);
xlabel('Iteration')
ylabel('Total Time [s]');
set(gca,'Xtick',[m:M]);
set(gca,'XtickLabel',[m-1:M-1]);

print -depsc generic1.eps;

figure;
z = squeeze(Z(1,3,m:M)); semilogy(m:M,z,'k','Linewidth',2); hold on;
z = squeeze(Z(2,3,m:M)); semilogy(m:M,z,'r','Linewidth',2);
z = squeeze(Z(3,3,m:M)); semilogy(m:M,z,'g','Linewidth',2);
z = squeeze(Z(4,3,m:M)); semilogy(m:M,z,'b','Linewidth',2);
z = squeeze(Z(5,3,m:M)); semilogy(m:M,z,'m','Linewidth',2); hold off;
xlim([m,M]);
ylim([1e-2,1]);
xlabel('Iteration')
ylabel('Relative L^2 Output Error');
set(gca,'Xtick',[m:M]);
set(gca,'XtickLabel',[m-1:M-1]);
legend('Full-Order Inversion','Original Algorithm','Data-Misfit Enhanced','Monte-Carlo Enhanced','Monte-Carlo & Data-Misfit Enhanced','Location','NorthEast');

print -depsc generic2.eps;

% State, Param, Online Speed Up vs Full, Offline Speed Up vs Original , Total Speed Up vs Full

T = [ (m:M)' , (m:M)' , squeeze(Z(1,2,1))./squeeze(Z(5,2,m:M)) , squeeze(Z(2,1,m:M))./squeeze(Z(5,1,m:M)) , squeeze(Z(1,2,1))./squeeze(Z(5,1,m:M)+Z(5,2,m:M)) ]




figure;
subplot(1,3,[1]);
%z = squeeze(Z(6,2,m:M)); semilogy(m:M,z,'k','Linewidth',2); hold on;
z = squeeze(Z(7,1,m:M)); semilogy(m:M,z,'r','Linewidth',2); hold on;
z = squeeze(Z(8,1,m:M)); semilogy(m:M,z,'g','Linewidth',2);
z = squeeze(Z(9,1,m:M)); semilogy(m:M,z,'b','Linewidth',2);
z = squeeze(Z(10,1,m:M)); semilogy(m:M,z,'m','Linewidth',2); hold off;
xlim([m,M]);
ylim([1e-1,1e3]);
xlabel('Iteration')
ylabel('Offline Time [s]');
set(gca,'Xtick',[m:M]);
set(gca,'XtickLabel',[m-1:M-1]);

subplot(1,3,[2]);
z = squeeze(Z(6,2,m:M)); semilogy(m:M,z,'k','Linewidth',2); hold on;
z = squeeze(Z(7,2,m:M)); semilogy(m:M,z,'r','Linewidth',2);
z = squeeze(Z(8,2,m:M)); semilogy(m:M,z,'g','Linewidth',2);
z = squeeze(Z(9,2,m:M)); semilogy(m:M,z,'b','Linewidth',2);
z = squeeze(Z(10,2,m:M)); semilogy(m:M,z,'m','Linewidth',2); hold off;
xlim([m,M]);
ylim([1e-1,1e3]);
xlabel('Iteration')
ylabel('Online Time [s]');
set(gca,'Xtick',[m:M]);
set(gca,'XtickLabel',[m-1:M-1]);
title('fMRI Connectivity Model');

subplot(1,3,[3]);
z = squeeze(Z(6,2,m:M)); semilogy(m:M,z,'k','Linewidth',2); hold on;
z = squeeze(Z(7,1,m:M)+Z(7,2,m:M)); semilogy(m:M,z,'r','Linewidth',2);
z = squeeze(Z(8,1,m:M)+Z(8,2,m:M)); semilogy(m:M,z,'g','Linewidth',2);
z = squeeze(Z(9,1,m:M)+Z(9,2,m:M)); semilogy(m:M,z,'b','Linewidth',2);
z = squeeze(Z(10,1,m:M)+Z(10,2,m:M)); semilogy(m:M,z,'m','Linewidth',2); hold off;
xlim([m,M]);
ylim([1e-1,1e3]);
xlabel('Iteration')
ylabel('Total Time [s]');
set(gca,'Xtick',[m:M]);
set(gca,'XtickLabel',[m-1:M-1]);

print -depsc fmri1.eps;

figure;
z = squeeze(Z(6,3,m:M)); semilogy(m:M,z,'k','Linewidth',2); hold on;
z = squeeze(Z(7,3,m:M)); semilogy(m:M,z,'r','Linewidth',2);
z = squeeze(Z(8,3,m:M)); semilogy(m:M,z,'g','Linewidth',2);
z = squeeze(Z(9,3,m:M)); semilogy(m:M,z,'b','Linewidth',2);
z = squeeze(Z(10,3,m:M)); semilogy(m:M,z,'m','Linewidth',2); hold off;
xlim([m,M]);
ylim([1e-2,1]);
xlabel('Iteration')
ylabel('Relative L^2 Output Error');
set(gca,'Xtick',[m:M]);
set(gca,'XtickLabel',[m-1:M-1]);
legend('Full-Order Inversion','Original Algorithm','Data-Misfit Enhanced','Monte-Carlo Enhanced','Monte-Carlo & Data-Misfit Enhanced','Location','NorthEast');

print -depsc fmri2.eps;
T = [ 5*(m:M)' , (m:M)' , squeeze(Z(6,2,1))./squeeze(Z(10,2,m:M)) , squeeze(Z(7,1,m:M))./squeeze(Z(10,1,m:M)) , squeeze(Z(6,2,1))./squeeze(Z(10,1,m:M)+Z(10,2,m:M)) ]

