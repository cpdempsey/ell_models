function cellname = makeName(modeltype)
    intypes = unique(modeltype);
    intypes(strcmp(intypes,'none'))=[];
    
    cellname = [intypes{:}];