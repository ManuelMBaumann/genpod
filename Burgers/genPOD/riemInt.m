function [ int ] = riemInt( fvec,t0,te )
%[ int ] = riemInt( fvec,t0,te ) computes the integral of the vector fvec 
%over the intervall t0,te by means of interpolated left and right hand Riemann sums
%fvec is vector containing EQUALLY distributed function values

n = length(fvec)-1;
int = 0.5*sum(fvec(1:n)+fvec(2:n+1))*(te-t0)/n;

end
