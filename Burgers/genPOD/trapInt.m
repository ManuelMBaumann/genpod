function err = trapInt( fvec,t0,te )


n = length(fvec)-1;
dt = (te-t0)/n;


dtvec = dt*ones(1,n);
trapvec = 0.5*(fvec(1:end-1) + fvec(2:end));
err = sum(dtvec.*trapvec);

end
