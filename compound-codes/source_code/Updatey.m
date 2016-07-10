function [yout]=Updatey(yin,Vcon,Ccon2,nd,ud)
if(ud==0)
    yout=yin;
elseif(ud==1)
    e=numel(Vcon);
    E=zeros(e+1,1);
    E(Vcon(nd,:))=1;
    yout=mod(sum([E(Ccon2) yin],2),2);
end
