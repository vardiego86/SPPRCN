function visualize_migration_trange_color(in_mat,conc,molecule,t1,t2)
% Plots position (XYZ) of particles in cartesian space for the specified
% time period. Color is set based on concentration of the specified 
% molecule.
%
% Arguments:
% - in_mat:     Matrix with individual cell positions in time
% - conc:       Matrix with molecular concentrations in time
% - molecule:   Molecule number (1-22)
% - t1:         starting time step
% - t2:         ending time step
%
% Diego A. Vargas
% March 12, 2015

global L X Y Z

scene = 1;

hfig = figure;
set(gcf,'Color',[1,1,1])

for t = t1:t2
    
    if molecule <= 22
        scatter3(in_mat(:,X,t),in_mat(:,Y,t),in_mat(:,Z,t),50,conc(:,molecule,t)','filled');
    elseif molecule >= 23   % [AJ] or sigmaAJ
        scatter3(in_mat(:,X,t),in_mat(:,Y,t),in_mat(:,Z,t),50,conc(:,t)','filled');
    end
    
    xlabel('X (\mum)','FontSize',12,'FontWeight','bold');
    ylabel('Y (\mum)','FontSize',12,'FontWeight','bold');
    zlabel('Z (\mum)','FontSize',12,'FontWeight','bold');
    
%     title(strcat('time:',{' '},num2str(ceil(scene/6)),' h'),'FontSize',12,'FontWeight','bold');
    title(strcat('time step: ',num2str(t)),'FontSize',12,'FontWeight','bold');
    
    axis([0 L 0 L 0 L]);
    colorbar;

    F(scene) = getframe(hfig);
    scene = scene+1;
    
end

%movie2avi(F,'a_3_b_3_video.avi','fps',6);