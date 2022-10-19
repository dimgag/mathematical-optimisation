function [points,values]=steepest_descent_method(f,x0,eps,maxiter,debug)
%
% SD - The steepest descent algorithm for unconstrained minimization.
%
% Use is made of the routine LNSRCH that carries out an efficient
% approximate line minimization as described in Fletcher, pp.33-40.
%
% The stopping criterion employed is that of Marquardt, see STPCRIT.
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
fx=f(x);
points=[x];
values=[fx];

rho=1/64; sigma=2*rho; % Accuracy parameter for line search

i=0;
stop=0;
while ((i<maxiter) && (stop==0)), % steepest descent iteration.
   i=i+1;
   if strcmp(debug,'on'),
      disp(['i=',num2str(i),', x=[',num2str(x'),'], f(x)=',num2str(fx)]);
   end;
   gfx=f(gradientinit(x)); % compute value and gradient
   fx=gfx.x; % compute derivative and gradient
   s=-gfx.dx'; % steepest descent search direction.
   alpha=line_search(f,x,s,rho,sigma,debug); % line search.
   new_x=x+alpha*s; % iteration step.
   new_fx=f(new_x);
   points=[points,new_x];
   values=[values,new_fx];
   stop=marquardt_stopping_criterion(x,new_x-x,eps); % convergence detection.
   x=new_x;
   fx=new_fx;
end;

%
% End of function steepest_descent.
