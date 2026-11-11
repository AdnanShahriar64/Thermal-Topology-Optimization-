function[Me,H_Source_per_H]=Elemental_Capacitance_Matrix_Generator(Element_Information,Mesh_Information,Material_Information,FE_Information_Region)

    SpcH=Material_Information.SpcH;
    rho=Material_Information.rho;
    B=FE_Information_Region.B;
    dJB=FE_Information_Region.dJB;
    W=Element_Information.W_Elemental;
    Shape=Element_Information.Shape_Elemental;
    Type=Mesh_Information.Type;
    nel=Mesh_Information.nel_Model;
    ndof=Mesh_Information.ndof_Model;
    neq=Mesh_Information.neq_Model;
    LM=Mesh_Information.LM_Model;
    
    nnp_Max=Mesh_Information.nnp_max_Model;
    ngp_Max=Mesh_Information.ngp_max_Model;
    ngp_Elemental=Element_Information.ngp_Elemental;
    nnp_Elemental=Element_Information.nnp_Elemental;


    Me=zeros(nnp_Max,nnp_Max,nel);   % Elemental Mass Matrix preallocation

    H_Source_per_H=zeros(nnp_Max,nnp_Max,nel);
    
    for e=1:nel
        
        gpLimit=ngp_Elemental(Type(1,e));
        nn=nnp_Elemental(Type(1,e));

        for gp=1:gpLimit
            
            N=Shape(gp,1:nn,Type(1,e));

            Me(1:nn,1:nn,e)=Me(1:nn,1:nn,e)+SpcH*rho*(N'*N)*dJB(gp,e)*W(gp,Type(1,e));
            
            H_Source_per_H(1:nn,gp,e)=N'*dJB(gp,e)*W(gp,Type(1,e));
            
        end
        
    end
    
end
