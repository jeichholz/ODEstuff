function expr=structsubs(expr,struct)
    fields=fieldnames(struct);
    for i=1:numel(fields)
        expr=subs(expr,fields{i},getfield(struct,fields{i}));
    end
end