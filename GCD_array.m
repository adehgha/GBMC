function [Agcd] = GCD_array(Array,Order)
Tmp=0;
if(nargin > 1)
    if any(strcmpi(Order,{'rows', 'columns', 'col'}))==1
        if any(strcmpi(Order,{'columns', 'col'}))==1
            Array = Array'; Agcd = (GCD_array(Array, 'rows'))';
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
Agcd = gcd(Array(:,1), Array(:,2));
for i=1:size(Array,2)-2
    Agcd = gcd(Agcd, Array(:,i+2));
end
%---------------------------------------------------------------------
%%%%% Check for passing a column!
if Tmp==1
    Agcd=Agcd';
end
end
