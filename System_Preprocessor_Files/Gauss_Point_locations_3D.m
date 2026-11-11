function[gp_X,gp_Y,gp_Z]=Gauss_Point_locations_3D(Element_Information,Mesh_Information)

    nnp_Elemental=Element_Information.nnp_Elemental;
    ngp_Elemental=Element_Information.ngp_Elemental;
    Type=Mesh_Information.Type;
    IEN=Mesh_Information.IEN;
    Nodes=Mesh_Information.Nodes;
    Shape=Element_Information.Shape_Elemental;
    
    SIEN=size(IEN);
    nel 	= SIEN(2);

%% Memory Preallocation
    gp_X=zeros(ngp_Elemental(Type(1,1)),nel);
    gp_Y=gp_X;
    gp_Z=gp_X;

%% Each element calculation  

for e=1:nel

    nn=nnp_Elemental(Type(1,e)); % No of Nodes
    xyz=[Nodes(IEN(1:nn,e),2) Nodes(IEN(1:nn,e),3) Nodes(IEN(1:nn,e),4)];
    
    gp_X(1:ngp_Elemental(Type(1,e)),e)=Shape(1:ngp_Elemental(Type(1,e)),1:nnp_Elemental(Type(1,e)),Type(1,e))*xyz(:,1);
    gp_Y(1:ngp_Elemental(Type(1,e)),e)=Shape(1:ngp_Elemental(Type(1,e)),1:nnp_Elemental(Type(1,e)),Type(1,e))*xyz(:,2);
    gp_Z(1:ngp_Elemental(Type(1,e)),e)=Shape(1:ngp_Elemental(Type(1,e)),1:nnp_Elemental(Type(1,e)),Type(1,e))*xyz(:,3);
    
        
end  

end
