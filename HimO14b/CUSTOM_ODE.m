function y = CUSTOM_ODE(f,h,T,z,u,p)
% CUSTOM_ODE (sample custom ode solver for emgr)
% by Christian Himpe, 2014 ( http://gramian.de )
% released under BSD 2-Clause License ( opensource.org/licenses/BSD-2-Clause )
%*

 if(exist('OCTAVE_VERSION'))
  y = lsode(@(y,t) f(y,ones(size(u,1),1)*(t<=h),p),z,linspace(0,h*T,T))';
 else,
  y = deval(ode45(@(t,y) f(y,ones(size(u,1),1)*(t<=h),p),[0 h*T],z),linspace(0,h*T,T));
 end;
