n=63;
k=36;
Gx=bchgenpoly(n,k);
gx=de2bi(double(Gx.x));
Px=zeros(1,n+1);
Px([1 n+1])=1;
[hx,rem]=gfdeconv(Px,gx',2);


H=zeros(n-(length(hx)-1),n);

for i=1:n-length(hx)+1
   H(i,i:i-1+length(hx))=hx;
end
k=n-length(H(:,1));

filename=sprintf('BCH_%d_%d.txt',n,k);
fid=fopen(filename,'w');
fprintf(fid,'[');
for i=1:numel(H(:,1))
  fprintf(fid,'[');
    for j=1:numel(H(1,:))-1
    fprintf(fid,'%d,',H(i,j));    
    end
    if i~=numel(H(:,1))
     fprintf(fid,'%d],',H(i,j+1));    
    else 
     fprintf(fid,'%d]]',H(i,j+1));   
    end
end
fclose(fid);
