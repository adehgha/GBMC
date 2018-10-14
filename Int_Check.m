function Cond = Int_Check(Var,Precis)
%%% This function checks the elements of Matrix var to be 
%%% integral with the Precision of 1e-(Precis).
Tval = 1/10^(Precis); T1 = abs(Var); Cond = (abs(T1 - round(T1)) < Tval);
% % Cond=rem(round(Var*Tval)/Tval,1) == 0;
end