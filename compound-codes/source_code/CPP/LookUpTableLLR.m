% Function Definition
clear all
clc
step=.002;
LLR_MAX=25;
x=step/2:step:LLR_MAX;
y=log((exp(x)+1)./(exp(x)-1));
%y(1)=y(2);
%plot(x,y)
%%
% t=0.01;
% y(floor(t/.01)+1)
% log((exp(t)+1)./(exp(t)-1))
z=y';
size_of_LUT=length(z);

fid = fopen('LUT','w');
%fprintf( fid,'%f %d %d\n',step,LLR_MAX,size_of_LUT);
fprintf( fid,'%f\n', z);
fclose(fid);

fid = fopen('LUT_Parameters','w');
fprintf( fid,'%f %d %d\n',step,size_of_LUT,LLR_MAX);
fclose(fid);
