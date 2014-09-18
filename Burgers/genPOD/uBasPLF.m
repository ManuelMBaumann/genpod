function [ uBasPLF ] = uBasPLF( n,x0,xe,N )
%UBASPLF gives the n-th function corresponding to a basis of  
%L^2(x0,xe) of hierarchical piecewise linear functions, see example in
%MSchmidts diss
%   n       = logical index of basis function
%   x0,xe   = starting and endpoint of the interval
%   N       = number of sample points

x = linspace(x0,xe,N);

% uBasPLF = polyval(ChebyshevPoly(n-1),x);

if n==1
    uBasPLF = linterp([x0 xe],[ 1  0],x);
    return
end

if n==2
    uBasPLF = linterp([x0 xe],[0 1],x);
    return
end

if n==3
    uBasPLF = linterp([x0 (x0+xe)*0.5 xe],[0 1 0],x);
    return
end

% uBasPLF = sin((n-2)*pi*x/(xe-x0));

l2 = floor(log2(n-2));
absInt = linspace(-x0,xe,2^(l2+1)+1);
ordInt = zeros(2^(l2+1)+1,1);
ordInt(2) = 1;
ordInt = circshift(ordInt,(n-2-2^l2)*2);
uBasPLF = linterp(absInt,ordInt,x);
    
return
