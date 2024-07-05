function [expr,renamedFunctions,renamedStrings]=funsToSyms(expr,maxorder)

    if ~exist('maxorder','var') || isempty(maxorder)
        maxorder=5;
    end

    if isa(expr,"symfun")
        expr=formula(expr);
    end

    %Find the symfuns. 
    funlist=findSymType(expr,"diff");

    renamedFunctions=dictionary();
    renamedStrings=dictionary();
    
    expr=mapSymType(expr,"diff",@renamingDiffs);

    expr=mapSymType(expr,"symfun",@renamingSymFuns);



    function out=renamingDiffs(A)
            %Find the name of the function which we are differentiating. 
            f=findSymType(A,'symfun');
            if numel(f) ~= 1
                error('Error in renaming %s, %d symbolic functions found when one expected.\n',string(A),numel(f));
            end
            %Get the name of the function in a string. 
            fname=symFunType(f);
            %Find the order of the derivative.  Interstingly, you can just
            %count the commas in the string version of A.  Maybe.  
            order=sum(char(string(A))==',');
            
           
            %Now make the string representation. 
            out=['D',char(fname),repmat('t',1,order)];
            renamedFunctions(A)=out;
            renamedStrings(out)=A;
    end

    function out=renamingSymFuns(A)
        
        %Ok, here it gets hard.  Undefined functions are of type symtype,
        %those are easy.  However, if we do diff(f(x(t)),t) we get back a
        %thing like D(f)(x(t)) whose type is also symfun.  :(

        %We might be able to tell the two apart using the class of A, but
        %I'm not sure.  For now, first get the string representaiton of A. 

        str=char(symFunType(A));
        
        %If str contains "D(" then we know we have a derivative.

        if str(1)=='D' && str(2) =='('
            %Then we need to figure out the order of the derivative. Count
            %the number of open parens.
            order=sum(str=='(');

            %Get the name of the thing we are differntiating. Find the
            %first close paren and the last open paren. 
            lo=find(str=='(');
            lo=lo(end);

            fc=find(str==')');
            fc=fc(1);

            %The string between lo and fc is the name of the function being
            %differentiated. 

            fname=str(lo+1:fc-1);

            %Append the appropriate number of p's for prime.
            out=[fname,repmat('p',1,order)];


        else
            %This case is easy.  We let symFunType get the name for us and
            %we are done. 
            out=str;
        end

        renamedFunctions(A)=out;
        renamedStrings(out)=A;
    end


end

    


