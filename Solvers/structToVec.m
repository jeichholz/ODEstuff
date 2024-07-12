function [vec,isSet]=structToVec(structs,varToIdx)
    if isstruct(structs)
        structs={structs};
    end
    vec=zeros(numel(varToIdx.keys),1);
    isSet=vec;
    for i=1:numel(structs)
        struc=structs{i};
        fields=fieldnames(struc);
        for j=1:numel(fields)
            f=fields{j};
            vec(varToIdx(f))=getfield(struc,f);
            isSet(varToIdx(f))=1;
        end
    end

