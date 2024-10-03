function [new_conc,junc] = find_junc(dist,conc,junc,sigma_aj)
% find junctions between proteins based on proximity of cells, allowed
% number of neighbors, membrane E-cadherin/B-catenin complexes per cell,
% and adhesivity factor(sigma) of the complexes.
%
% update (E-cad/B-cat)M in concentration matrix
%
% February 2, 2015
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global N r_o a_neigh

X18 = 18;                   % variable number in concentratin matrix for (E-cad/B-cat)M
SigM = 21;                  % variable number in concentration matrix for adhessiveness of (E-cad/B-cat)M

r_n = r_o;                  % distance above which two cells are not considered as neighbors

new_conc = conc;

[sorted_dist,dist_i] = sort(dist,2,'ascend');      % dist_i holds cell2 number (column), row indicated cell1

% defining which AJs between cells are lost with none gained (i.e. lost cell-cell connections):
delta_junc = zeros(N,N)-junc(:,:,1);

% for j = 2:(a_neigh+1)    % 2:(a_neigh+1) (to exclude cell itself)
%    for cell = 1:N   
%       
%        if sorted_dist(cell,j)<r_n % maybe add check to see if junc is not 0 and do fill in matrtix based on symmetry
%            
%            % New junction formation (or dissociation)
%            new_aj = ((k24(conc(cell,SigM))+k24(conc(dist_i(cell,j),SigM)))/2)*(conc(cell,X18))*(conc(dist_i(cell,j),X18))...
%                - ((k_24(sigma_aj(cell))+k_24(sigma_aj(dist_i(cell,j))))/2)*(junc(cell,dist_i(cell,j),1));
%            
%            junc(cell,dist_i(cell,j),2) = junc(cell,dist_i(cell,j),2) + new_aj;
%            
%            % Updating concentration matrix (X18)
%            new_conc(cell,X18) = new_conc(cell,X18) - new_aj;
%            
%            % Making sure broken AJs are not counted twice in formation of(E-cad/B-cat)M
%            delta_junc(cell,dist_i(cell,j)) = 0;
%            
%        end
%        
%    end
% end

for j = 2:(a_neigh+1)    % 2:(a_neigh+1) (to exclude cell itself) 
      
   ind_n = find(sorted_dist(:,j) < r_n);    % Indices of cells within r_n (to loop through less cells).
   
   for cell = 1:length(ind_n)

       % New junction formation (or dissociation)
       new_aj = ((k24(conc(ind_n(cell),SigM))+k24(conc(dist_i(ind_n(cell),j),SigM)))/2)*(conc(ind_n(cell),X18))*(conc(dist_i(ind_n(cell),j),X18))...
           - ((k_24(sigma_aj(ind_n(cell)))+k_24(sigma_aj(dist_i(ind_n(cell),j))))/2)*(junc(ind_n(cell),dist_i(ind_n(cell),j),1));

       junc(ind_n(cell),dist_i(ind_n(cell),j),2) = junc(ind_n(cell),dist_i(ind_n(cell),j),2) + new_aj;

       % Updating concentration matrix (X18)
       new_conc(ind_n(cell),X18) = new_conc(ind_n(cell),X18) - new_aj;

       % Making sure broken AJs are not counted twice in formation of(E-cad/B-cat)M
       delta_junc(ind_n(cell),dist_i(ind_n(cell),j)) = 0;
       
   end
       
end

% From AJs to membrane cadherins
% for i = 1:N
%     for j = i:N
%    
%     if delta_junc(i,j) ~= 0
%         new_conc(i,X18) = new_conc(i,X18) - delta_junc(i,j);
%         new_conc(j,X18) = new_conc(j,X18) - delta_junc(j,i);
% 
%         junc(i,j,2) = 0;
%         junc(j,i,2) = 0;
%     end
%     
%     end
% end

% for i = 1:N
%    
%     C_n = delta_junc(i,i:N)~=0;
%     new_conc(i,X18) = new_conc(i,X18) - sum(delta_junc(i,C_n));
%     new_conc(C_n,X18) = new_conc(C_n,X18) - delta_junc(i,C_n)';
% 
%     junc(i,C_n,2) = 0;
%     junc(C_n',i,2) = 0;
%     
% end

for j = 1:N
   
    C_n = delta_junc(j:N,j)~=0;
    new_conc(j,X18) = new_conc(j,X18) - sum(delta_junc(C_n,j));
    new_conc(C_n,X18) = new_conc(C_n,X18) - delta_junc(C_n,j);

    junc(C_n,j,2) = 0;
    junc(j,C_n',2) = 0;
    
end


