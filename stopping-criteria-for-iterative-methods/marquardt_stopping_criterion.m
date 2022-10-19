 function stop=marquardt_stopping_criterion(x,step,eps)
%
% marquard_stopping_criterion - Marquardt's stopping criterion for detecting convergence.
% 
% The input parameters indicate the (old) iteration point 'x' and the 
% proposed 'step' (which is equal to -alpha*s), both as column vectors,
% and an desired accuracy eps.
%
% The output parameter 'stop' indicates whether convergence occurred.
%   stop=0 : continue the iteration process.
%   stop=1 : convergence detected.
%
% By Pieter Collins, Department of Knowledge Engineering, 
%   Maastricht University, 2012.
% 
% Modified from stpcrit.m by Ralf Peeters, Mathematics Dept., 
%   Universiteit Maastricht. Jan. 1999.
% 
%

stop=1;
n=length(x);
for i=1:n,
   stop=stop*(abs(step(i))<=eps*(abs(x(i))+eps*8));
end;

%
% End of function marquard_stopping_criterion.
