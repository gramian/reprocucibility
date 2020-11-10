function [time] = plot_alliterates(path,filename)

if nargin == 0
    path =  '../../results/projection_stability_ellipsoid/Polynomials/Discrete/alliterates/'
    filename = 'stock';
end

if ~isdir(path)
        error('The given directory does not exist')
end

tic;

S =load(filename);
Bmat = getfield(S,'Bmat');

[dim,nit] = size(Bmat);
dim =sqrt(dim);

if nit >50
    warning('The number of iterates is quite high')
    nit
    if nit >100
%        error('The number of iterates is really too high!!')
        nit = 100;
    end
end

%taking the original unstable problem
A = reshape(Bmat(:,1),dim,dim);
a = poly(A);

%taking the starting point of the problem
S = reshape(Bmat(:,2),dim,dim);
s = poly(S);

%cycling trough the iterates
for i =  98:nit
    i
    It = reshape(Bmat(:,i),dim,dim);
    it = poly(It);
    
    PlotStableSet(a,it,s);
    
    print(gcf,'-dpdf','-r300',sprintf('%s%s%i%s',path,'it_',i-2,'.pdf'))
    close;
    s = it;
end
time =toc;
    
    