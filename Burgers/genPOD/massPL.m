function [ massPL, mY, NU ] = massPL( N,nb,t0,te,FUN)
%MY computes the mass matrix for the spatial discretization functions of
%the output space y_h, if hierarchical piecewise linear functions are used

%       t0,te - inner product of L^2(t0,te)
%       N     - spatial discretization
%       nb    - no. basis functions

M=size(FUN,2); % =Nt
massPL = zeros(N,nb);
mY = zeros(nb);
NU = zeros(nb,M);

for j=1:nb
    %NU(j,:)= HaarWavelet(j,t0,te,M);
    NU(j,:)= uBasPLF(j,t0,te,M);
end

for i=1:N
    for j=1:nb
        massPL(i,j) = trapInt( FUN(i,:).*NU(j,:),t0,te );
    end
end

for i=1:nb
    for j=1:i
        mY(i,j) = trapInt( NU(j,:).*NU(i,:),t0,te );
        mY(j,i) = mY(i,j);
    end
end

return