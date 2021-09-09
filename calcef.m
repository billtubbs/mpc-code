%%%%%%%%%%%%%%%%%%%%
%%% É. Plamondon %%%
%%% A. Kettani   %%%
%%% A. Desbiens  %%%
%%%%%%%%%%%%%%%%%%%%

function [E,F]=calcEF(A,B,C,D,j)
%	Calcule les matrices de prediction deterministe E et F

[lA,cA]=size(A);
[lB,cB]=size(B);
[lC,cC]=size(C);
[lD,cD]=size(D);

if any(any(D))

	E=zeros(lC*j,cA);
	F=zeros(lC*j,cB*j+cD);
	
	E(1:lC,:)=C*A;
	F(1:lC,1:cB+cD)=[C*B D];
	
	for i=1:(j-1)
		E(lC*i+1:lC*(1+i),:)=E(lC*(i-1)+1:lC*i,:)*A;	
		F(lC*i+1:lC*(i+1),1:j*cB+cD)=[E(lC*(i-1)+1:lC*i,:)*B F(lC*(i-1)+1:lC*i,1:(j-1)*cB+cD)];
	end

else
	
	E=zeros(lC*j,cA);
	F=zeros(lC*j,cB*j);
	
	E(1:lC,:)=C*A;
	F(1:lC,1:cB)=[C*B];
	
	for i=1:(j-1)
		E(lC*i+1:lC*(1+i),:)=E(lC*(i-1)+1:lC*i,:)*A;	
		F(lC*i+1:lC*(i+1),1:j*cB)=[E(lC*(i-1)+1:lC*i,:)*B F(lC*(i-1)+1:lC*i,1:(j-1)*cB)];
	end

end