function [ alpha ] = line_search( f, x, s, rho, sigma, debug )
% line_search: Find alpha>0 such that x+alpha*s approximately minimises f
% Dispatches to one of several possible functions for performing the 
% actual search

    alpha = fletcher_line_search(f,x,s,rho,sigma,debug);
   %alpha = bisection_line_search(f,x,s,rho,sigma,debug);

end
