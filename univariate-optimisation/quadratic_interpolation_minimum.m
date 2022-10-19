function m = quadratic_interpolation_minimum(f,a,c,b)
    fa=f(a);
    fc=f(c);
    fb=f(b);

    m=c - ((c-a)*(c-a)*(f(c)-f(b))-(c-b)*(c-b)*(f(c)-f(a)))/((c-a)*(f(c)-f(b))-(c-b)*(f(c)-f(a)))/2;
