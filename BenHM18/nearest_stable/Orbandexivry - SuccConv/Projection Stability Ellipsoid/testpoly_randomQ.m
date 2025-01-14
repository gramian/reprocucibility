%test for the efficiecy of choosing a random Q

%It repeats the experiment "nite" times so
%that Q takes differents random values.
%% params for the tests
nite = 100;

%% definition of the problem
% delta = 0.1;
% a = [1 0 0 0 delta];
% % path = 'D:\Documents and Settings\fxorban\My Documents\fxorban\'
% 


a = poly([-2 -3 4 -0.5 0.1]);

a = a(end:-1:1);
a = a.';

n = length(a)-1;%polynomials are supposed monic so you don't take that coeff
Z = [zeros(1,n-1) 0;eye(n-1) zeros(n-1,1)];
en = [zeros(n-1,1);1];
A = Z - a(1:n)*en';


%% params for the algorithm
isdiscrete = 0;
verbose = 0;
compute = 0;
S = [];
e_lon =0;%cannot be used with polynomials yet!
roc  = 0;
savf = 'None';
path = '';

%% Processing

objval = zeros(nite,1);

for k = 1 : nite
    
    %the line for choosing Q random should be activated
    [bi]= StablePolyMain(a,S,isdiscrete,compute,verbose,e_lon,roc,savf,path);
    
    objval(k) = norm(a-flipud(bi.'),'fro');
end

mobj = mean(objval)
