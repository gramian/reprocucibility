%%
% This code is published under the BSD3-Clause License.
% Copyright (c) 2016, Jens.Saak
%%
%   All rights reserved.
%  Redistribution and use in source and binary forms, with or without modification,
%  are permitted provided that the following conditions are met:
%  1. Redistributions of source code must retain the above copyright
%   notice, this list of conditions and the following disclaimer.
%  2. Redistributions in binary form must reproduce the above copyright
%   notice, this list of conditions and the following disclaimer in the
%    documentation and/or other materials provided with the distribution.
%  3. Neither the name of the copyright holder nor the names of its
%   contributors may be used to endorse or promote products derived from
%    this software without specific prior written permission.
%  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
%  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
%  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
%  A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
%  HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
%  iNCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
%  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
%  OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
%  ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR
%  TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
%  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
%  DAMAGE.
% All rights reserved. 
% Jens Saak, MPI Magdeburg (c) 2016


% tested with 
%    Matlab 8.6.0.267246 (R2015b) Linux/distru 23-May-2016

% The example used to test the IRKA algorithm is the
% FOM example from the SLICOT Library
% http://slicot.org/20-site/126-benchmark-examples-for-model-reduction
% downloaded at 21.12.2015 by Joerg Fehr, ITM Uni Stuttgart
FOM_Ex = load('fom.mat');
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('%')
disp('% The FOM example from the Slicot collection')
disp('%')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

%% Test IRKA with odd shift number
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('% Multishift example -- reduced order 15')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
close all
nred = 15;
b = randn(size(FOM_Ex.B,2),nred);
c = randn(nred,size(FOM_Ex.C,1));
shifts = 1:nred-1;
shifts = [shifts 243 ]';

[V5,W5,A_red5,B_red5,C_red5,E_red5,s5]=IRKA(FOM_Ex.A,... % system Mat. FOM ex.
    FOM_Ex.B,...                                  % input  Mat. FOM ex.
    FOM_Ex.C,...                                  % output Mat. FOM ex.
    speye(size(FOM_Ex.A,1)),...                   % descriptor Matrix
    shifts',...                           % initial shift values
    b,...                                         % tangential direction for B
    c,...                                         % tangential direction for C
    25,...                                        % maximum iteration number
    1.0e-2,...                                    % convergence tolerance
    1);                                           % debug output on
%%
for i=1:5, 
    fig=figure(i);
    if exist('matlab2tikz','file')
        matlab2tikz(sprintf('figure_%d.tikz',i),...
            'width','.4\linewidth',...
            'height','.3\linewidth'); 
    else
        saveas(fig,sprintf('figure_%d',i),'epsc');
    end
end
