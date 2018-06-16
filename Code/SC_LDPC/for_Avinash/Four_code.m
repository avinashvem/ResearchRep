M=2^14;
J=750;
A=zeros(J,M);
%A=zeros(750,2^12);
Index = randperm(M,J);
for i=1:J
    for j=1:M
        A(i,j)=exp(1i*2*pi*(Index(i)-1)*(j-1)/(M));
    end
end