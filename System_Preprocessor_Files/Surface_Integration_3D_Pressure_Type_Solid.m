function[FE_Information_Surface]=Surface_Integration_3D_Pressure_Type_Solid(Element_Information,Mesh_Information)
%%
% Line='---------------------------------------------------------------------------------------------';
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
% P_surf_gp=zeros(Mesh_Information.nnp_max_Model*Mesh_Information.ndof_Model,Mesh_Information.ngp_max_Model,max(Mesh_Information.nsurf_Model),SNs);
dTeS_gp=zeros(3,Mesh_Information.ngp_max_Model,max(Mesh_Information.nsurf_Model),SNs);

P_surf=zeros(Mesh_Information.nnp_max_Model*ndof,max(Mesh_Information.nsurf_Model),SNs);
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
            
            Vectors=(transpose(J))*Par_Coordinate;

            dTe=(cross(Vectors(:,1),Vectors(:,2)));
            % dTe=dTe/(norm(dTe))^2;
            
            
            % Valx=Valx+(Surf_Shape(gp,1:nn,Surf,Type(1,e)))'*(Surf_Shape(gp,1:nn,Surf,Type(1,e)))*P;
            % hold on
            % plot(Valx)
            % pause
            % Valx=Valx+Surf_W(gp,Surf,Type(1,e))*(Surf_Shape(gp,1:nn,Surf,Type(1,e)))'*dTe(1);
            % Valy=Surf_W(gp,Surf,Type(1,e))*(Surf_Shape(gp,1:nn,Surf,Type(1,e)))'*dTe(2);
            % Valz=Surf_W(gp,Surf,Type(1,e))*(Surf_Shape(gp,1:nn,Surf,Type(1,e)))'*dTe(3);
            % plot(Valx)
            % hold on
            
            P_surf(1:3:nn*ndof-2,Scount,SN)=P_surf(1:3:nn*ndof-2,Scount,SN)+Surf_W(gp,Surf,Type(1,e))*(Surf_Shape(gp,1:nn,Surf,Type(1,e)))'*dTe(1);
            P_surf(2:3:nn*ndof-1,Scount,SN)=P_surf(2:3:nn*ndof-1,Scount,SN)+Surf_W(gp,Surf,Type(1,e))*(Surf_Shape(gp,1:nn,Surf,Type(1,e)))'*dTe(2);
            P_surf(3:3:nn*ndof-0,Scount,SN)=P_surf(3:3:nn*ndof-0,Scount,SN)+Surf_W(gp,Surf,Type(1,e))*(Surf_Shape(gp,1:nn,Surf,Type(1,e)))'*dTe(3);

            % P_surf_gp(1:3:nn*ndof-2,gp,Scount,SN)=Surf_W(gp,Surf,Type(1,e))*(Surf_Shape(gp,1:nn,Surf,Type(1,e)))'*dTe(1);
            % P_surf_gp(2:3:nn*ndof-1,gp,Scount,SN)=Surf_W(gp,Surf,Type(1,e))*(Surf_Shape(gp,1:nn,Surf,Type(1,e)))'*dTe(2);
            % P_surf_gp(3:3:nn*ndof-0,gp,Scount,SN)=Surf_W(gp,Surf,Type(1,e))*(Surf_Shape(gp,1:nn,Surf,Type(1,e)))'*dTe(3);

            dTeS_gp(:,gp,Scount,SN)=dTe;

            dTea=dTea+dTe;
            
            % fprintf('%15d   %15d    %25d   %25f\n',e,SN,gp,det(J))
            

        end
            
        dTeS(:,Scount,SN)=dTea/norm(dTea);
        
            
        %% https://scicomp.stackexchange.com/questions/26886/evaluating-the-surface-integral-in-an-fem-finite-elements-method-procedure

        SNo=SN;

end

% FE_Information_Surface.P_surf_gp=P_surf_gp;

FE_Information_Surface.dTeS_gp=dTeS_gp;
FE_Information_Surface.dTeS=dTeS;
FE_Information_Surface.P_surf=P_surf;



% fprintf('%s\n\n',Line)

end
