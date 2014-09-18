
ns = [17 33 65];
k1 = [2 4 8 16 32];

errM = zeros(length(ns),length(k1));
gerrM_pwl = zeros(length(ns),length(k1));
%gerrM_haar = zeros(1,length(k1));

c1 = 1;
leg = ['^', 'o', 's'];
for n=ns
    c2=1;
    for k=k1
        [errM(c1,c2), gerrM_pwl(c1,c2)] = MOR_driver(n,k);
        c2 = c2 + 1;
    end
    figure(666)
    semilogy(k1,errM(c1,:),['r-' num2str(leg(c1))]);
    hold on
    semilogy(k1,gerrM_pwl(c1,:),['k-' num2str(leg(c1))]);
    c1 = c1 +1 ;
end
    