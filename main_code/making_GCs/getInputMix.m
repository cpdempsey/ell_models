function [MF_input, Ws, cellname, modeltype] = getInputMix(GCparams,MF_indices)


MF_input   = [0 0 0];
Ws         = [0 0 0 0 0 0];


%pick the input types and set the probability of input spiking
clear modeltype;
emptycount=3;
modeltype=cell(3,1);
while(emptycount>2)
    emptycount=0;
    for i=1:3
        inclass = rand(1);
        if(inclass<GCparams.theta.e)
            modeltype{i} = 'early';
            MF_prob(i) = 1;
        elseif(inclass<(GCparams.theta.e+GCparams.theta.m))
            modeltype{i} = 'med';
            MF_prob(i) = 1;
        elseif(inclass<(GCparams.theta.e+GCparams.theta.m+GCparams.theta.t))
            modeltype{i} = 'tonic';
            emptycount=emptycount+1;
        elseif((inclass<(GCparams.theta.e+GCparams.theta.m+GCparams.theta.t+GCparams.theta.l))) 
            modeltype{i} = 'late';
            MF_prob(i) = rand(1);
        elseif((inclass<(GCparams.theta.e+GCparams.theta.m+GCparams.theta.t+GCparams.theta.l+GCparams.theta.p)))
            modeltype{i} = 'pause';
            MF_prob(i) = 1;
        else
            modeltype{i} = 'none';
            emptycount=emptycount+1;
        end
    end
end

cellname = makeName(modeltype);


Winds=sum(GCparams.Wsparse,3)~=0;
%pick the input cells
for i=1:3
        switch(modeltype{i})
            case 'none'  %generated GC receives no input!
                
                MF_input(i) = 0;
                Ws(i) = 0;
                Ws(i+3) = 0;
                
            case 'tonic' %generated GC receives only tonic input!
                
                MF_input(i) = MF_indices.(modeltype{i})(ceil(rand(1)*length(MF_indices.(modeltype{i}))));
                [Ws(i) , Ws(i+3)] = get_tonic_weights_fixed(GCparams);
                
            otherwise
                
                MF_input(i) = MF_indices.(modeltype{i})(ceil(rand(1)*length(MF_indices.(modeltype{i}))));

                Winds_filt=Winds;
                Winds_filt(:,setdiff(1:size(Winds_filt,2),MF_indices.(modeltype{i})))=0;
                [Wi,Wj]=find(Winds_filt);
                usecell = ceil(rand(1)*length(Wi)); %pick an input of that class to use for weights.

                Ws(i)=GCparams.Wsparse(Wi(usecell),Wj(usecell),1);
                Ws(i+3)=GCparams.Wsparse(Wi(usecell),Wj(usecell),2);
                
        end
end

