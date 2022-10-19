function m = golden_section(f,a,c,b,e)
    w=(3-sqrt(5))/2;
    fa=f(a);
    fb=f(b);
    fc=f(c);

    assert(a<c)
    assert(c<b)
    assert(fa>fc)
    assert(fb>fc)

    assert(e>0)

    while(b-a>e),

        if (c-a)>(b-c) || ((c-a)==(b-c) && fa<fb),
            x=c-w*(c-a);
            fx=f(x);
            if fx<fc,
                b=c; fb=fc; c=x; fc=fx;
            else
                a=x; fa=fx;
            end
        else
            x=c+w*(b-c);
            fx=f(x);
            if fx<fc,
                a=c; fa=fc; c=x; fc=fx;
            else
                b=x; fb=fx;
            end
        end
        [x;fx]
        [a,c,b;fa,fc,fb]
    end

    m=c;
