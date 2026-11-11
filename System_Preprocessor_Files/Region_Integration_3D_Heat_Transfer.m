function[FE_Information_Region]=Region_Integration_3D_Heat_Transfer(Element_Information,Mesh_Information)
%%
Line='-------------------------------------------------------------------------------------';
% fprintf('\n%s\n',Line)
% fprintf('%s\n','Calculating Region Integration parameters...');
% fprintf('%s\n',Line)
% fprintf('%15s   %25s    %25s\n','Element no','Integration point no','Jacobian determinant') 
% fprintf('%s\n',Line)

Nodes=Mesh_Information.Nodes;
IEN=Mesh_Information.IEN;
Type=Mesh_Information.Type;
DShape_Elemental=Element_Information.DShape_Elemental;
W_Elemental=Element_Information.W_Elemental;
Shape_Elemental=Element_Information.Shape_Elemental;
nel 	= Mesh_Information.nel_Model;

%% B matrix structure
%   [dN1/dx 0       0       dN2/dx 0      0      0...|
%    0      dN1/dy  0       0      dN1/dy 0      0...| Integration point 1
%    0      0       dN1/dz  0      0      dN2/dz 0...|

%    dN1/dx 0       0       dN2/dx 0      0      0...|
%    0      dN1/dy  0       0      dN1/dy 0      0...| Integration point 2
%    0      0       dN1/dz  0      0      dN2/dz 0...|

%% Memory Preallocation
    B=zeros(length(W_Elemental(:,end))*3,Mesh_Information.nnp_max_Model,nel);
    dJB=zeros(length(W_Elemental(:,end)),nel);    

%% Each element calculation

for e=1:nel

        nn=Element_Information.nnp_Elemental(Type(1,e));                                                         % No of Nodes
        xyz=[Nodes(IEN(1:nn,e),2) Nodes(IEN(1:nn,e),3) Nodes(IEN(1:nn,e),4)];   % Setting up xyz vector

        
        for gp=1:Element_Information.ngp_Elemental(Type(1,e))                            % Goes through all Gauss points
                % Determine B matrix ---------------------------------
                % Jacobian Matrix
                    J=DShape_Elemental(:,1:Type(2,e),gp,Type(1,e))*xyz;
                % Determinant of Jacobian Matrix
                    dJ=det(J);
                    if dJ<=0
                        error('Warning! Mesh is invalid, negative determinant of jacobian detected, element no=%d, integration point no=%d\n',e,gp);
                    end
                % Derivative matrix in global coordinate
                    CDsx3=inv(J)*DShape_Elemental(:,1:Type(2,e),gp,Type(1,e));
                % B Matrix
                    %    Be=zeros(3,nn);
                    %for k=1:Type(2,e)
                    %    Be(:,3*k-2:3*k)= [CDsx3(1,k)       0           0      ;
                    %                          0       CDsx3(2,k)       0      ;
                    %                          0             0       CDsx3(3,k)];
                    %end
                    
                    
                    if sum(abs(Shape_Elemental(gp,:,Type(1,e))))<10^-3 || dJ<=0 || sum(find(isnan(CDsx3)))>0
                        fprintf(2,'Warning! Weight is set to 0, element no=%d, integration point no=%d\n',e,gp);
                    else
                    
                        B(gp*3-2:gp*3,1:nn,e)=CDsx3;
                        dJB(gp,e)=dJ;
                    end
                    
                   % fprintf('%15d   %25d    %25f\n',e,gp,dJ)
                    
        end

end

FE_Information_Region.B=B;
FE_Information_Region.dJB=dJB;

% fprintf('%s\n\n',Line)

end
