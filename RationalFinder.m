
function [output]=RationalFinder(input,tolerance,Order)

% % -------Handling the case of "close to integer" matrices---------
check=Int_Check(input,17);
if all( check == true(size(input)))
    input=round(input);
    output=input;
else
    
    %%% By default it flattens the array (if nargin < 3)
    DefaultTol=1e-6;
    %--------------------------
    switch nargin
        case 1
            tol=DefaultTol; Sz=size(input);
            input = reshape(input,1,Sz(1)*Sz(2));
        case 2
            tol = tolerance; Sz=size(input);
            input = reshape(input,1,Sz(1)*Sz(2));
        case 3
            Switch=0;
            tol=tolerance;
            ErrMsg=strcat(['Not a valid input. Please choose either',...
                ' "rows" or "columns" keys for this function.']);
            Keys=find(strcmpi(Order,{'rows', 'columns', 'col'}));
            if ~Keys
                error(ErrMsg);
            end
            if (Keys)==2 || (Keys)==3
                input = input';
            end
            %handling the case of asking a column vector with the 'row' key by
            %mistake.
            if Keys==1 && size(input,2)==1
                input = input';Switch=1;
            end
    end
    
    % % input = (round(input./tol).*tol);
    %--------------------------
    %%%%% Check for passing a column!
    
    % % tmp=(input~=0);
    Tol1 = 1e-6; tmp=(abs(input)> Tol1);
    % % Vec=1e20*ones(size(input,1),size(input,2));
    Vec=2*max(abs(input(:)))*ones(size(input,1),size(input,2));
    Vec(tmp)=input(tmp);
    MIN=min(abs(Vec),[],2);
    
    input = input./repmat(MIN,1,size(input,2));
    [N,D]=rat(input,tol);
    LCM_rows = LCM_array(D, 'rows');
    LCM_mat=repmat(LCM_rows,1,size(input,2));
    Rounded=(N.*LCM_mat)./D;
    
    output = Rounded;
    % % %%Eliminating the multiples.
    % % [ArrayGCD] = GCD_array(Rounded,'rows');
    % % GCDMatrix=repmat(ArrayGCD,1,size(input,2));
    % % output = Rounded./GCDMatrix;
    %--------------------------
    if nargin < 3
        output=reshape(output,Sz);
    else
        if (Keys)==2 || (Keys)==3
            output = output';
        end
        if Keys==1 && Switch==1
            output = output';
        end
    end
end
end