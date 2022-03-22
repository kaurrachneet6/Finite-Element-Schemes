function [uhh] = compute_GLOBAL_MATRIX(h,N)
K = zeros(N);
M= zeros(N);
m = 0:h:1
f=zeros(N,1);
for i=1:N
     k1 = zeros(N);
    if i==1
    k1(i,i)=1;
    
    k1
    else
    k1(i-1,i-1) =1;
    k1(i-1,i)=-1;
    k1(i,i-1)=k1(i-1,i);
    k1(i,i)=1;
           
    end
 K = K+k1; 
end

for i=1:N
     m1 = zeros(N);
    if i==1
    m1(i,i)=1/3;
    m1
    else
    m1(i-1,i-1) =1/3;
    m1(i-1,i)=1/6;
    m1(i,i-1)=m1(i-1,i);
    m1(i,i)=1/3;
    end
 M = M+m1; 
end

for k = 1:N
 kEl = k:k+1;
 fe = getElementForceVector4(m(kEl));
 if k==1 f(k) = f(k)+fe(2);
 else
 f(kEl-1) = f(kEl-1)+fe;
 end
end
K = K./h
M= M.*h
A=K+M
f
uh=A^(-1)*f
uhh = [ 0; uh ]

end

