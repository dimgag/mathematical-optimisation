function alpha=fletcher_line_search(f,x,s,rho,sigma,debug)
%
% fletcher_line_search - Line search for the function specified by f, 
%          starting from the iteration point x (column vector)
%          in the descent direction specified by s (column vector).
%          Acceptability of alpha is controlled by the line search
%          parameters rho and sigma, which must satisfy 0<rho<1/2 and
%          rho<sigma<1.
%
% Use is made of the sophisticated approach described in 
%     R. Fletcher, (1987). Practical Methods of Optimization.
%     John Wiley and Sons, Chichester (etc.). 2nd edition, pp. 33-40. 
% This approach requires function values as well as (directional)
% derivatives of the function to be minimized, as supplied by the 
% routines indicated by 'fname' and 'dfname'. 
%
% The procedure consists of two phases: a bracketing phase and a 
% sectioning phase. Both involve the minimization of a cubic 
% interpolation polynomial over a bounded interval, which is carried 
% out by the subroutine CUBICMIN.
%
% The string variable 'debug' controls the output generated for the
% purpose of easy debugging. If debug='off', no printed debugging 
% output is generated.
%
% Upon leaving this routine, the output parameter 'alpha' yields an 
% acceptable stepsize parameter for the given search direction s.
%
% The parameters tau1, tau2 and tau3 control 
% the initialization (tau1) and the sectioning (tau2 and tau3) of the 
% bracket that contains acceptable values for alpha.
%
% By Ralf Peeters, Mathematics Dept., Universiteit Maastricht.
% Jan. 1999.
%

% computation of initial values.
gf0=f(x+gradientinit(0)*s);
f0=gf0.x; f0prime=gf0.dx;

% Display arguments
if strcmp(debug,'on'),
    disp(['  fletcher_line_search:']);
    disp(['    x=[',num2str(x'),'], s=[',num2str(s'),'], rho=',num2str(rho),', sigma=',num2str(sigma)]);
    disp(['    f(0)=',num2str(f0),', df(0)=',num2str(f0prime)]);
end

if nargin < 4, rho=1/16; end % it should hold that 0.0 < rho < 0.5.
if nargin < 5, sigma=2*rho; end % it should hold that rho < sigma < 1.0.
if nargin < 6, debug='off'; end
   
% specification of fixed parameters of the line search routine.
tau1=9.0;    % it should hold that tau1 > 1.0.
tau2=0.1;    % it should hold that 0.0 < tau2 < tau3 <= 0.5,
tau3=0.5;    % and it is advisable that tau2 <= sigma.


fbar = f0 + 1024 * f0prime;
mu=(fbar-f0)/(rho*f0prime); % upper bound on the bracket size.

alpha0=0;
alpha1=1; % initial estimate for alpha.

terminate=0; % (0: continue, 1: end b-phase, 2: end s-phase).
fold=f0;
foldprime=f0prime;
j=1;
while terminate==0, % bracketing phase of the line search routine.
   if strcmp(debug,'on'),
      disp(['      bracketing step ',num2str(j)]);
   end;
   f1=f(x+alpha1*s);
   if f1<=fbar,
      terminate=2;
   else
      if ((f1>f0+alpha1*rho*f0prime) || (f1>=fold)),
         a=alpha0;
         fa=fold;
         faprime=foldprime;
         b=alpha1;
         fb=f1;
         gfb=f(x+gradientinit(alpha1)*s); 
         fbprime=gfb.dx;
         terminate=1;
      else
         gf1=f(x+gradientinit(alpha1)*s); 
         f1prime=gf1.dx;
         if abs(f1prime)<=-sigma*f0prime,
            terminate=2;
         else
            if f1prime>=0,
               a=alpha1;
               fa=f1;
               faprime=f1prime;
               b=alpha0;
               fb=fold;
               fbprime=foldprime;
               terminate=1;
            else
               if mu<=2*alpha1-alpha0,
                  alpha0=alpha1;
                  alpha1=mu;
                  fold=f1;
               else
                  low=2*alpha1-alpha0;
                  upp=min(mu,alpha1+tau1*(alpha1-alpha0));
                  alphanew=cubicmin(alpha0,alpha1,fold,f1,foldprime,f1prime,low,upp);
                  alpha0=alpha1;
                  fold=f1;
                  foldprime=f1prime;
                  alpha1=alphanew;
               end;
            end;
         end;
      end;
   end;
   j=j+1;
end;

j=1;
while terminate==1, % sectioning phase of the line search routine.
   if strcmp(debug,'on'),
      disp(['      sectioning step ',num2str(j),': a=',num2str(a),' b=',num2str(b)]);
   end;
   low=a+tau2*(b-a);
   upp=b-tau3*(b-a);
   alphanew=cubicmin(a,b,fa,fb,faprime,fbprime,low,upp);
   alpha1=alphanew;
   f1=f(x+alpha1*s);
   if ((f1>f0+rho*alpha1*f0prime) || (f1>=fa)),
      b=alpha1;
      fb=f1;
   else
      gf1=f(x+gradientinit(alpha1)*s); 
      f1prime=gf1.dx;
      if abs(f1prime)<=-sigma*f0prime,
         terminate=2;
      else
         if (b-a)*f1prime>=0,
            b=a;
            fb=fa;
         end;
         a=alpha1;
         fa=f1;
      end;
   end;
   j=j+1;
end;

if strcmp(debug,'on'),
    gf1=f(x+gradientinit(alpha1)*s); 
    disp(['    result: alpha=',num2str(alpha1),', f(alpha)=',num2str(gf1.x),', df(alpha)=',num2str(gf1.dx)]);
end;

alpha=alpha1; % output of this line search routine.

%
% End of function LNSRCH.

function alpha=cubicmin(a,b,fa,fb,faprime,fbprime,low,upp)
%
% CUBICMIN - Minimization of a cubic interpolation polynomial over a 
%            bounded interval [low,upp].
%
% The cubic interpolation polynomial passes through the points (a,fa)
% and (b,fb), where it has the derivatives faprime and fbprime, 
% respectively. The routine handles both the situations a<b and b<a.
%
% By Ralf Peeters, Mathematics Dept., Universiteit Maastricht.
% Jan. 1999.
%

% conversion of input arguments to the interval [0,1] instead of [a,b].
f0=fa;
f0p=faprime*(b-a);
f1=fb;
f1p=fbprime*(b-1);
p(1)=(low-a)/(b-a);
p(2)=(upp-a)/(b-a);

% computation of coefficients of the cubic interpolation polynomial.
c0=f0;
c1=f0p;
c2=3*(f1-f0)-2*f0p-f1p;
c3=f0p+f1p-2*(f1-f0);

% computation of the values at the bounds and at the stationary points.
v(1)=c0+c1*p(1)+c2*p(1)^2+c3+p(1)^3;
v(2)=c0+c1*p(2)+c2*p(2)^2+c3+p(2)^3;
discr=c2^2-3*c3*c1;
if ((abs(c3)>1e-8) && (discr>=0)), % safeguard against the case c3=0.
   p(3)=(-c2-sqrt(discr))/(3*c3);
   p(4)=(-c2+sqrt(discr))/(3*c3);
   v(3)=c0+c1*p(3)+c2*p(3)^2+c3+p(3)^3;
   v(4)=c0+c1*p(4)+c2*p(4)^2+c3+p(4)^3;
end;

% determination of the minimum point alpha in the interval [low,upp].
minpoint=p(1);
vmin=v(1);
for i=2:length(v),
   if a<b, % usual ordering of interval points a and b.
      if ((p(i)>=p(1)) && (p(i)<=p(2)) && (v(i)<vmin)),
         minpoint=p(i);
         vmin=v(i);
      end;
   else % order of interval points a and b is interchanged.
      if ((p(i)<=p(1)) && (p(i)>=p(2)) && (v(i)<vmin)),
         minpoint=p(i);
         vmin=v(i);
      end;
   end;
end;
alpha=a+minpoint*(b-a); % conversion of 'minpoint' to original coordinates.

%
% End of function CUBICMIN.
