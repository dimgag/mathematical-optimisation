function [ alpha ] = bisection_line_search( f, x, s, rho, sigma, debug )
% bisection_line_search -  Find alpha>0 such that x+alpha*s approximately 
%    minimises f using a bisection method in [0,1].
%    No checks are made concerning the line search parameters rho and
%    sigma.

    if nargin < 4, rho=1/16; end
    if nargin < 6, debug='off'; end
        
    f0=f(x);
    a=0; b=1; 
    fa=f(x+a*s); fb=f(x+b*s);
    
    c=(a+b)/2; fc=f(x+c*s);
    if strcmp(debug,'on'),
        disp(['f(',num2str(a),')=',num2str(fa)]);
        disp(['f(',num2str(b),')=',num2str(fb)]);
        disp(['f(',num2str(c),')=',num2str(fc)]);
    end
    while abs(fb - fa) > rho * ( f0-fc )
        if fa<fb,
            b=c; fb=fc;
        else,
            a=c; fa=fc;
        end
        c=(a+b)/2; fc=f(x+c*s);
        if strcmp(debug,'on'),
            disp(['f(',num2str(c),')=',num2str(fc)]);
        end
    end
    
    alpha=c;
end
