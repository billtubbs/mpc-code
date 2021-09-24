function [F,M]=Diophantine_equation(C,D,Hp)
%% Recursive solution of the Diophantine equation
%
%     [F,M]=Diophantine_equation(C,D,Hp):
%
%           C/D = F + z^(-j)*M/D, for j = 1 to Hp
%           C and D in decreasing power of z^-1
%
%     Example for 2-step ahead: 
%        (1 - 0.2^z-1) / (1 -0.8z^-1) = 1 + 0.6z^-1 + 0.48z^-2/(1 - 0.8z^-1)
%        [F,M]=Diophantine_equation([1 -0.2],[1 -0.8],2)
%        We obtain the result for 1-step ahead and 2-step ahead
%
% André Desbiens and Kaddour Najim, Copyright (c), April 1993


Atild=1;
Pn=C;
Pd=D;
N2=Hp;

E=zeros(N2,N2);
E(1,1)=Pn(1)/(Pd(1)*Atild(1));
AtPd=conv(Atild,Pd);
temp=subspoly(Pn,E(1,1)*AtPd);
temp=temp(2:length(temp));
cF=length(temp);
if cF==0
    temp=0;
    cF=1;
end
F(1,1:cF)=temp;
for j=2:N2
    r=F(j-1,1)/(AtPd(1));
    for i=1:cF
        if i+1>length(AtPd)
            atild=0;
        else
            atild=AtPd(i+1);
        end
        if i+1>cF
            f=0;
        else
            f=F(j-1,i+1);
        end
        F(j,i)=f-atild*r;
    end
    E(j,1:j)=[E(j-1,1:j-1) r];
end

M=F;
F=E;
 
function P=subspoly(A,B)
%% Returns the subtraction of two polynomials.
%
%           P=subspoly(A,B) returns the polynomial P = A - B.  The polynomials
%           A and B do not need to be the same length.
% André Desbiens and Kaddour Najim, Copyright (c), April 1993

lA=length(A);
lB=length(B);
if lA>lB
    B=[B zeros(1,lA-lB)];
else
    A=[A zeros(1,lB-lA)];
end
P=A-B;
