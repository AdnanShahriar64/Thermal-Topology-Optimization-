function[Ke]=Elemental_Conductivity_Matrix_Generator(Element_Information,Mesh_Information,Material_Information,FE_Information_Region,Dm)
%%

% rho=Material_Information.rho;
B=FE_Information_Region.B;
dJB=FE_Information_Region.dJB;
W=Element_Information.W_Elemental;
% Shape=Element_Information.Shape_Elemental;
Type=Mesh_Information.Type;
nel=Mesh_Information.nel_Model;
% ndof=Mesh_Information.ndof_Model;
% neq=Mesh_Information.neq_Model;
% LM=Mesh_Information.LM_Model;

nnp_Max=Mesh_Information.nnp_max_Model;
ngp_Max=Mesh_Information.ngp_max_Model;
ngp_Elemental=Element_Information.ngp_Elemental;
nnp_Elemental=Element_Information.nnp_Elemental;

Ke=zeros(nnp_Max*1,nnp_Max*1,nel);   % Elemental Stiffness Matrix preallocation


for e=1:nel
    
    gpLimit=ngp_Elemental(Type(1,e));
    nn=nnp_Elemental(Type(1,e));

    for gp=1:gpLimit
        
        B_matrix=B(gp*3-2:gp*3,1:nn,e);
        
        Ke(1:nn,1:nn,e)=Ke(1:nn,1:nn,e)+B_matrix'*Dm*B_matrix*dJB(gp,e)*W(gp,Type(1,e));
      
    end

    
end


end
