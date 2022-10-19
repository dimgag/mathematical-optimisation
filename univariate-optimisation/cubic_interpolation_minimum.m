function m = cubic_interpolant_critical(f,df,a,b)
    fa=f(a);
    fb=f(b);
    dfa=df(a);
    dfb=df(b);

    al=(fb-fa)/(b-a);
    be=(dfb-dfa)/2;
    ga=(dfa+dfb)/2;

    m=(a+b)/2 - (b-a)/2 * (sqrt(be^2+3*(al-ga)*(3*al-ga))-be)/(3*(al-ga));
