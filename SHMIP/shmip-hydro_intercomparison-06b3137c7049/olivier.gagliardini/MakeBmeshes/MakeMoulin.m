% Matlab Script to insert the 101 elements needed for the Moulin
% in the mesh.boundary file 
% 
clear;
Err = 1.0e-2 ;

for bb=1:5
    mesh = ['mesh_B' num2str(bb)] ;

    % Read the moulin coordinates
    Moulin=dlmread(['./B' num2str(bb) '_M.xy']);
    % Open the mesh.nodes file to find the nodes number
    nodes=dlmread([mesh,'/mesh.nodes']);
    % Read the BC file
    bc=dlmread([mesh,'/mesh.boundary']);
    % Read the header file
    header=dlmread([mesh,'/mesh.header']);

    ms=size(Moulin,1);
    Nnode=size(nodes,1);
    nBC=size(bc,1);
    nheader=size(header,1);

    NodeMoulin(1:ms) = 0 ;
    for ii=1:ms
        xm = Moulin(ii,1) ;
        ym = Moulin(ii,2) ;
        for jj=1:Nnode
           xn = nodes(jj,3) ;
           yn = nodes(jj,4) ;
           if (abs((xm-xn)^2)+abs((ym-yn)^2)) < Err
               NodeMoulin(ii) = jj 
               break ;
           end ;
        end ;
    end ;

    % Test if each Moulin has been associated a mesh node
    if any(NodeMoulin==0)
        Disp('Error - No nodes correspond to some moulin locations')
        break ;
    end ;

    % Write the 101 BC at the end of the mesh.boundary file
    MaxBC = max(bc(:,2)) 

    % Rewrite the file and add the 101 elements
    fid1=fopen([mesh,'/mesh.boundary'],'w');
    for ii=1:nBC
        fprintf(fid1,'%g %g %g %g %g %g %g \n', bc(ii,:));
    end ;
    jj=1 ;
    for ii=1:ms
        fprintf(fid1,'%g %g %g %g %g %g \n', nBC+ii,MaxBC+jj,1,0,101,NodeMoulin(ii));
        jj = jj + 1 ;
    end ;

    % Change the header file
    fid1=fopen([mesh,'/mesh.header'],'w');
    fprintf(fid1,'%g %g %g \n', header(1,1),header(1,2),header(1,3)+ms);
    fprintf(fid1,'%g \n', header(2,1)+1);
    fprintf(fid1,'%g %g \n',101,ms);
    for ii=3:nheader
        fprintf(fid1,'%g %g \n',header(ii,1),header(ii,2));
    end;
end


