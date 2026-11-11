function [K] = Assembly(Ke,ndof,Mesh_Information,Element_Information)


Type=Mesh_Information.Type;
nel=Mesh_Information.nel_Model;
nnp=Mesh_Information.nnp_Model;
nnp_Elemental=Element_Information.nnp_Elemental;

if ndof==1
    LM=Mesh_Information.IEN;
elseif ndof==3
    LM=Mesh_Information.LM_Model;
else
    error('This assembly file is not built for any DOF except 1 and 3')
end


K 	= zeros(nnp*ndof);

%% Assembly of K and M
%fprintf('\nAssemblying...')

for e = 1:nel
    nn=nnp_Elemental(Type(1,e));
    for loop1 = 1:nn*ndof
        i = LM(loop1,e);
        for loop2 =  1:nn*ndof
            j = LM(loop2,e);
            K(i,j) = K(i,j) + Ke(loop1,loop2,e);
        end
    end
    %fprintf('Element=%d   Type=%d\n',e,Type(e));
end
%fprintf('Done\n')


end
