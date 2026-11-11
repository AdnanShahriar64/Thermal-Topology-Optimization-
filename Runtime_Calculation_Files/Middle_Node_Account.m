function[Nodes]=Middle_Node_Account(Nodes,Pair_for_Interpolate)


for n=1:size(Pair_for_Interpolate,1)

    Upper_Node=Pair_for_Interpolate(n,1);
    Middle_Node=Pair_for_Interpolate(n,2);
    Lower_Node=Pair_for_Interpolate(n,3);

    Node_Upper_xyz=Nodes(Upper_Node,2:4);

    Node_Lower_xyz=Nodes(Lower_Node,2:4);

    Node_Middle_xyz=(Node_Upper_xyz+Node_Lower_xyz)/2;

    Nodes(Middle_Node,2:4)=Node_Middle_xyz;

end



end
