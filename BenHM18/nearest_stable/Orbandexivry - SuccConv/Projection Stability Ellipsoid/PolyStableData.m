%Data for polynomials

%1.2-13.1x+8.85x^2+20.65x^3+8.4x^4+x^5 (roots are -2 -3 -4 0.5 0.1)
% a = [1.2 -13.1 8.85 20.65 8.4 1];
%(roots are -2 -3 -4 -0.5 0.1)

%%


%%
%first case (continuous)
% a = poly([-2 -3 4 -0.5 0.1]); 
% a = a*i

%second case (continuous)
% a = poly([-2.5 -3 4 -0.5 0.1]);
% a = a%+i*ones(1,n)*1e-9;
% path =  '../../results/projection_stability_ellipsoid/Polynomials/Continuous/deg5_realcoeff_pbisrealequaltooneenforced/'
% path =  '../../results/projection_stability_ellipsoid/Polynomials/Continuous/deg5_realcoeff_becomescomplexcoeff/'

%third case (discrete)
% a = poly([0.1 1.1 -0.2 -0.8 0.9]);
% a = a+i*ones(1,n)*1e-9;

% % fourth case(discrete)
% a = poly([0.1+0.3i -0.7-0.5i -0.5-0.95i -0.8 0.9-0.05i]);
% % path = '../../results/projection_stability_ellipsoid/Polynomials/Discrete/deg5_cmplxcoeff/'
% path = 'D:\Documents and Settings\fxorban\My Documents\Dropbox\Doctorat\notes\figures'
% namefig = 'sol_barrier_beta-nu'

% % fifth case (discrete)
% a = poly([-2 -3 4 -0.5 0.1]); %does not work with SPReflection because -2 and 0.5 are mirror roots with respect to the unit circle

% sixth case (continuous)
a = poly([-1 -1 -1 0.1]);%work in the continuous time case because of the bijection that occurs between the continuous and discrete time case
path = '.'

% % seventh case (continuous)
% a = poly([1 -3 -3.5 -6 -7]);%neds to do something for the roots '1', in the bijection discrete-continuous

% % Heighth case (continuous)
n = 4;
delta = -0.1;
a = [1 zeros(1,n-1) -delta];
% % path = 'D:\Documents and Settings\fxorban\My Documents\fxorban\'
% % path = 'D:\Documents and Settings\fxorban\My Documents\Dropbox\Doctorat\notes\figures'
% % namefig = 'sol_companion_barrier_beta-101'

% % Ninth case (continuous)
% r = [-1 -5 2 4 6 8];
% a = poly(r);



%%

isdiscrete = 0;
verbose = 1;
compute = 0;
% S = poly([-1 -1 -1 -1]).';
% S = S(end:-1:1);
S = [];
e_lon =0;%cannot be used with polynomials yet!
roc  = 0;
savf = '';
namefig='';
