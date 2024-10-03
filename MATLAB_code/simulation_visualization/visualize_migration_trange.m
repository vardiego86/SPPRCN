function visualize_migration_trange(in_mat,t1,t2)
% Plots position (XYZ) of particles in cartesian space for the specified
% time period

global L X Y Z

figure;

for t = t1:t2
    
    scatter3(in_mat(:,X,t),in_mat(:,Y,t),in_mat(:,Z,t),4,'filled');
    xlabel('x position');
    ylabel('y position');
    zlabel('z position')
    title(strcat('time step:',{' '},num2str(t)));
    axis([0 L 0 L 0 L]);
    
    pause(0.05);
end