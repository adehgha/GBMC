function [Alcm] = LCM_array(Array,Order)
Tmp=0;
Array(Array==0)=1;
if(nargin > 1)
    if any(strcmpi(Order,{'rows', 'columns', 'col'}))==1
        if any(strcmpi(Order,{'columns', 'col'}))==1
            Array = Array'; Alcm = (LCM_array(Array, 'rows'))';
            Tmp=1;
        end
    else
        error('Not a valid input. Please choose either "rows" or "columns" keys for this function.');
    end
end
%%%%% Check for passing a column!
if size(Array,2)==1
    Array=Array';
end
%---------------------------------------------------------------------
Alcm = lcm(Array(:,1), Array(:,2));
for i=1:size(Array,2)-2
    Alcm = lcm(Alcm, Array(:,i+2));
end
%---------------------------------------------------------------------
%%%%% Check for passing a column!
if Tmp==1
    Alcm=Alcm';
end
end