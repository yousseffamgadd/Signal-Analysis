clc;
clear all;
close all;
a=9;
b=1;

Ns=2;
t=0:0.1:30;
x= a*sin(0.5*b*pi*t);
N=1:5;
l=2.^N;

 Xs=x(1:Ns:end);%sampling
 ts=t(1:Ns:end);
 
 for i=1:length(l)
   [q(i),var_theo(i),var_prac(i),Sqnr_theo(i),sqnr_prac(i),finalquantized]=Q(l(i),Xs,a)
 
 end 

plot(t,x);
hold on
plot(ts,finalquantized);
title('Input and Output Signals');

 figure();
 plot(l,q);
 title('mean');
 figure();
 plot(l,var_prac);
 title('practical variance');
 figure();
 plot(l,var_theo);
 title('theortical variance');
 figure();
 plot(l,sqnr_prac);
 title('practical sqnr');
 figure();
 plot(l,Sqnr_theo);
 title('theortical sqnr');
 function  [q,var_theo,var_prac,Sqnr_theo,sqnr_prac,finalquantized]=Q(l,Xs,a)
  vmax=a
  vmin=-a
  delta=(vmax-vmin)/l
  xq=vmin+(delta/2):delta:vmax-(delta/2)%quantization
  for i=1:length(Xs)
    error=abs(Xs(i)-xq)
    y=find(error==min(error))
    selectedlevel=y(1)
    quantizedvalue(i)=xq(selectedlevel)
     encoded(i)=selectedlevel-1
  end %errors and encoding
q=mean(abs(Xs-quantizedvalue));%mean
var_theo=(delta)^2/12;%theortical variance
var_prac=var(Xs-quantizedvalue);%practical variance
Sqnr_theo=a^2/var_theo;% theortical sqnr
sqnr_prac=a^2/var_prac;%practical sqnr
msg=unique(encoded);
for i=1 : length(msg)
prob(i)=length(find(encoded==msg(i)))/length(encoded);
end
[dict,avg]=huffmandict(msg,prob);
sourceenc=huffmanenco(encoded,dict);% Huffman Source encoder.
sourcedec=huffmandeco(sourceenc,dict);%Huffman Source decoder
for i=1:length(Xs)
decoded(i)=sourcedec(i)+1;
finalquantized(i)=xq(decoded(i));
end 
SI=-log2(prob);
H=-sum(prob.*log2(prob));
v=ceil(SI);
L=sum(prob.*v);
E=(H/L)*100 ;%efficiency
comp=log2(l)/L ;%compression rate
end
function [Xs,ts]= Sampler(Ns,t,x) 
a=9;
b=1;
t=0:0.1:30;
x=a*sin(0.5*b*pi*t);
Ns=2;
ts=t(1:Ns:end);
Xs=x(1:Ns:end);

end %sampling function




