function[Type,gp_X,gp_Y,gp_Z]=Gauss_Point_Update_locations_3D_L(Nodes,IEN,Type,Shape,W)




%% Update Type by putting gauss point numbers
    SIEN=size(IEN);
    nel 	= SIEN(2);
    for e=1:nel

        eT=Type(1,e);

        gpn=length(W(:,1));
        for i=1:length(W(:,1))
            if W(i,eT)==0
                gpn=i-1;
                break
            end
        end

        Type(3,e)=gpn;

    end
  
    
%% Get gp locations    
% Memory Preallocation
    gp_X=zeros(max(Type(3,:)),nel);
    gp_Y=gp_X;
    gp_Z=gp_X;

% Each element calculation  

for e=1:nel

    nn=Type(2,e); % No of Nodes
    xyz=[Nodes(IEN(1:nn,e),2) Nodes(IEN(1:nn,e),3) Nodes(IEN(1:nn,e),4)];
    
    gp_X(1:Type(3,e),e)=Shape(1:Type(3,e),1:Type(2,e),Type(1,e))*xyz(:,1);
    gp_Y(1:Type(3,e),e)=Shape(1:Type(3,e),1:Type(2,e),Type(1,e))*xyz(:,2);
    gp_Z(1:Type(3,e),e)=Shape(1:Type(3,e),1:Type(2,e),Type(1,e))*xyz(:,3);
    
end  


end
