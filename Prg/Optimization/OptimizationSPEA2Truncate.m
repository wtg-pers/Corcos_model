function [P_arch,P_arch_obj] = OptimizationSPEA2Truncate(P_arch,P_arch_obj,D,Na,Ntot)

D_copy_2 = D;
s = size(P_arch,1);

for i = 1 : s-Na %one iteration per extra element in archive

    D_copy = D_copy_2(:,1:s-i+1);
    k=1;
    b=ones(1,1000);
    
    while length(b) > 1 && k <= Ntot % as long as their is more than one min
        clear b;
        b = [];
        a = min(D_copy(k,:)); %look for min
       
        for j = 1 : size(D_copy,2) % look for all possible equals distance
            if D_copy(k,j) == a
                b = [b j];
            end
        end
        
        clear D_copy;
        D_copy = [];
        parfor j2 = 1:length(b) % save only minimal k-th distance elements
            D_copy =[D_copy D(:,b(j2))];
        end
        k=k+1;
    end
    
    if k<Ntot+1 %if terminated normally, delete the minimal distance element.
        P_arch(b,:) = [];
        P_arch_obj(b,:) = [];
        D_copy_2(b,:) = [];
    else %if not, it means that 2 elements are identical so delete the first one
        P_arch(b(1),:) = [];
        D_copy_2(b(1),:) = [];
        P_arch_obj(b(1),:) = [];
    end

end
