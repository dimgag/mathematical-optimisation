function [points,values]=rank_one_quasi_newton_method(f,x0,eps,maxiter,debug)
%
% rank_one_quasi_newton_method - The rank-one quasi-Newton method for unconstrained minimization.
%
% Use is made of a routine fletcher_line_search that carries out an efficient
% approximate line minimization as described in Fletcher, pp.33-40.
%
% The stopping criterion employed is that of Marquardt.
%
% Function:                      f.
% Starting point:                x0 (column vector).
% Desired accuracy in x:         eps (number).
% Maximum number of iterations:  maxiter.
%
% Output printing for debugging purposes is controlled by the
% string variable 'debug'. When debug='off', no output is printed.
%
% The output parameters 'points' and 'values' contain the points that 
% constitute the optimization path and the corresponding function 
% values in their columns. 
%
% By Pieter Collins, Department of Knowledge Engineering,
% Maastrict University, 2012
%
% Modified from sd.m By Ralf Peeters, Mathematics Dept., Universiteit Maastricht.
% Jan. 1999.
%

% Set default arguments if necessary
if nargin < 5, debug='off'; end
if nargin < 4, maxiter=128; end

% initialization.
x=x0;
dfx=f(gradientinit(x)); % compute value and gradient
fx=dfx.x;
gx=dfx.dx';

n=size(x,1);
H=eye(n);

points=[x];
values=[fx];

rho=0.25; sigma=0.5; % Accuracy parameter for line search

i=0;
stop=0;
while ((i<maxiter) && (stop==0)), % steepest descent iteration.
   i=i+1;
   s= - H * gx; % Search direction.

   if strcmp(debug,'on'),
      i,x,fx,gx,norm(gx),H,s
   end;
   assert(s'*gx<0);
   
   alpha=line_search(f,x,s,rho,sigma,debug); % line search.
   assert(alpha>0);
   new_x=x+alpha*s; % iteration step.
   new_dfx=f(gradientinit(new_x));
   new_fx=new_dfx.x;
   new_gx=new_dfx.dx';
   assert(new_fx<fx)
   
   points=[points,new_x];
   values=[values,new_fx];
   
   stop=marquardt_stopping_criterion(x,new_x-x,eps); % convergence detection.
   
   if norm(new_gx)<1e-15, stop=1; end; % zero gradient detection
   
   % Update H
   delta = alpha*s;
   gamma = new_gx - gx;
   u = delta - H * gamma;
   
   % Update H, checking for division by zero
   if u'*gamma ~= 0,
       H = H + u*u'/(u'*gamma);
   end
   
   x=new_x;
   fx=new_fx;
   gx=new_gx;

end;

%
% End of function steepest_descent.
