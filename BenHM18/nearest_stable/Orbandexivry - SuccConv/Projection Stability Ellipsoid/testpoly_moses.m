%testpoly_moses runs the tests found in the article "Determining the
%Closest Stable polynomial to an Unstable one" by Randolph Moses et al.
%Some results were saved in a .mat files titled "ixx" where "xx" is a
%number giving the values of tests successfully executed.
%Created : 30/03/2011
%History : 17/04/2012 : enabling the possibility of reusing unstable
%                       instances that were previously generated.
%Last modification : 17/04/2012

close all

%% parameters

%parameter for this file
reuse_old_values = 1;

%for calling StablePolyMain
roc = 0;
elon = 0;
verbose = 0;
compute = 0;
savf = 'none';
path = '';

%original stable instance (from Moses)
a  = [1 -2.7607 3.8106 -2.6535 0.9238];

%% initialisation
nexp = 50;
nc = length(a);
unstable = zeros(nexp,nc);
stable = zeros(nexp,nc);
devbf = zeros(nexp,nc);
devaf = zeros(nexp,nc);
dist = zeros(nexp,1);
distbf = zeros(nexp,1);
distaf = zeros(nexp,1);

if reuse_old_values
    %load values of a,unstable and runstable.
    nexp = 50;
    try
        load 50_original_unstable_polynomials_and_their_roots
    catch
        warning('Tried to open the file but didn''t work')
        rethrow(lasterror);
    end
end

i = 1;

%processing montecarlo simulation
while i <= nexp
    try
        if reuse_old_values
            p = unstable(i,:);
        else
            y = random('norm',0, 1e-2,[1 4]);
            p = a+[0 y];
        end

        b = StablePolyMain(fliplr(p),[],1, compute, verbose,elon,roc,savf,path);
        if any(abs(b-p)>1e2*eps) == 1 || reuse_old_values
            unstable(i,:) = p;
            stable(i,:) = b;
            devbf(i,:) = p-a;
            devaf(i,:) = b-a;
            dist(i) = norm(p-b);
            distbf(i) = norm(p-a);
            distaf(i) = norm(b-a);
            %             PlotStableSet(p,b,a);
            %             if dist(i) > 1
            %                 keyboard %for enabling the analysis of the results
            %             end
            i = i +1;
        end
    catch
        %dont do anything just continue to computing the statistics
        if strcmp(getfield(lasterror,'identifier'),'MATLAB:license:checkouterror')
            rethrow(lasterror);
        end
        warning('Error in stablepolymain');
        i
    end
end

%computing the statistics
mdist = mean(dist(1:nexp))
mdistbf = mean(distbf(1:nexp))
mdistaf = mean(distaf(1:nexp))
stdevbf = std(devbf(1:nexp,[2:nc]))
stdevaf = std(devaf(1:nexp,[2:nc]))
meanp = mean(unstable(1:nexp,[2:nc]))
meanb = mean(  stable(1:nexp,[2:nc]))

%% plot of the roots of the polynomials
runstable = zeros(nexp,nc-1);
rstable = zeros(nexp,nc-1);
figure;
hold on;
for i = 1 : nexp
    runstable(i,:) = roots(unstable(i,:));
    h(i) = plot(real(runstable(i,:)),imag(runstable(i,:)),'o');
    set(h(i),'markersize',5,'markeredgecolor','black','linewidth',1.5);
end
%mise en page
box on;
set(gca,'fontsize',14)
xlabel('Re(z)')
ylabel('Im(z)')
set(gca,'xlim',[0.5 0.85],'ylim',[0.5 0.85]);
set(gca,'xtick',[0.5 0.85],'ytick',[0.5 0.85]);
x = [0:0.1:360];
h1 = plot(cosd(x),sind(x));
set(h1,'Color','black','Linewidth',3)
set_plot(5,5,1)

figure;
hold on;
for i = 1 : nexp
    rstable(i,:) = roots(stable(i,:));
    hh(i) = plot(real(rstable(i,:)),imag(rstable(i,:)),'o');
    set(hh(i),'markersize',5,'markeredgecolor','black','linewidth',1.5);
end
%adding the starting point
r = roots(a);
sp = plot(real(r),imag(r),'o')
set(sp,'markersize',5,'linewidth',3,'markeredgecolor',[0.75 0.75 0.75]);
%mise en page
box on;
set(gca,'fontsize',14)
xlabel('Re(z)')
ylabel('Im(z)')
set(gca,'xlim',[0.5 0.85],'ylim',[0.5 0.85]);
set(gca,'xtick',[0.5 0.85],'ytick',[0.5 0.85]);
x = [0:0.1:360];
h2 = plot(cosd(x),sind(x));
set(h2,'Color','black','Linewidth',3)
set_plot(5,5,1)

