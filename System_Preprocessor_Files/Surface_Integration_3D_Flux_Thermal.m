function[FE_Information_Surface]=Surface_Integration_3D_Flux_Thermal(Element_Information,Mesh_Information)
%%
Line='---------------------------------------------------------------------------------------------';
% fprintf('\n%s\n',Line)
% fprintf('%s\n','Calculating Surface Integration parameters');
% fprintf('%s\n',Line)

% fprintf('%15s   %15s    %25s   %25s\n','Element no','Local Surface no','Integration point no','Jacobian determinant') 
% fprintf('%s\n',Line)



S_Type=Mesh_Information.S_Type;
SS=size(S_Type);
Nodes=Mesh_Information.Nodes;
IEN=Mesh_Information.IEN;
Type=Mesh_Information.Type;
Surf_Dshape=Element_Information.Surf_Dshape_Elemental;
Surf_W=Element_Information.Surf_W_Elemental;
Surf_Shape=Element_Information.Surf_Shape_Elemental;
ndof=Mesh_Information.ndof_Model;

SNs=max(S_Type(1,:));
T_surf_gp=zeros(Mesh_Information.nnp_max_Model,Mesh_Information.ngp_surf_max_Model,SS(2));
dTeS_gp=zeros(3,Mesh_Information.ngp_surf_max_Model,SS(2));

T_surf=zeros(Mesh_Information.nnp_max_Model,max(Mesh_Information.nsurf_Model),SNs);
Ke_surf=zeros(Mesh_Information.nnp_max_Model,Mesh_Information.nnp_max_Model,Mesh_Information.nel_Model,SNs);
dTeS=zeros(3,max(Mesh_Information.nsurf_Model),SNs);


%%
SNo=0;
Scount=0;
for i=1:SS(2)

        SN=S_Type(1,i);     % Extract Surface Number of the Model

        e=S_Type(2,i);      % Extract Element no
        nn=Element_Information.nnp_Elemental(Type(1,e));       % No of Nodes in the elements
        xyz=[Nodes(IEN(1:nn,e),2) Nodes(IEN(1:nn,e),3) Nodes(IEN(1:nn,e),4)];   % Extract xyz
        Surf=S_Type(3,i);    % Extract Which surface of the element

        if SN>SNo
            Scount=0;
        end
            Scount=Scount+1;


        Surf_gp_for_Surf=Element_Information.Surf_ngp_Elemental(Surf,Type(1,e));
        

        Par_Coordinate=Element_Information.Par_Coord_Elemental(:,:,Surf,Type(1,e));

            
        dTea=[0 0 0]';
        for gp=1:Surf_gp_for_Surf
           
            J=Surf_Dshape(:,1:nn,gp,Surf,Type(1,e))*xyz;
            
            if det(J)<=0
                fprintf(2,'Warning! Mesh is invalid, negative determinant of jacobian on surface detected, Element no=%d, Surface no=%d, Integration point no=%d\n',e,Surf,gp);
            end
            
            Vectors = J.' * Par_Coordinate;

            dTe = cross(Vectors(:,1),Vectors(:,2));
            
            dA = norm(dTe);
            
            T_surf(1:nn,Scount,SN)=T_surf(1:nn,Scount,SN)+Surf_W(gp,Surf,Type(1,e))*(Surf_Shape(gp,1:nn,Surf,Type(1,e)))'*dA;
            
            T_surf_gp(1:nn,gp,i)=Surf_W(gp,Surf,Type(1,e))*(Surf_Shape(gp,1:nn,Surf,Type(1,e)))'*dA;


            Ke_surf(1:nn,1:nn,e,SN)=Ke_surf(1:nn,1:nn,e,SN)+(Surf_Shape(gp,1:nn,Surf,Type(1,e)))'*(Surf_Shape(gp,1:nn,Surf,Type(1,e)))*Surf_W(gp,Surf,Type(1,e))*dA;
            
            dTea=dTea+dTe;
            
            dTeS_gp(:,gp,i)=dTe;

            % fprintf('%15d   %15d    %25d   %25f\n',e,SN,gp,det(J))

        end
            
        dTeS(:,Scount,SN)=dTea/norm(dTea);
        
            
        %% https://scicomp.stackexchange.com/questions/26886/evaluating-the-surface-integral-in-an-fem-finite-elements-method-procedure

        SNo=SN;

end

FE_Information_Surface.T_surf_gp=T_surf_gp;
FE_Information_Surface.dTeS_gp=dTeS_gp;
FE_Information_Surface.dTeS=dTeS;
FE_Information_Surface.Ke_surf=Ke_surf;
FE_Information_Surface.T_surf=T_surf;


% fprintf('%s\n\n',Line)

end
