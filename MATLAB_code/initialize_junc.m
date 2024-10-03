function [junc] = initialize_junc(dist,conc)
% find junctions between proteins based on proximity of cells, allowed
% number of neighbors, membrane E-cadherin/B-catenin complexes per cell,
% and adhesivity factor (sigma) of the complexes.
%
% February 2, 2015
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global N r_o a_neigh

X18 = 18;                   % variable number in concentration matrix for (E-cad/B-cat)M
SigM = 21;                  % variable number in concentration matrix for adhesvity factor of (E-cad/B-cat)M

r_n = r_o;                  % distance above which two cells are not considered as neighbors

junc = zeros(N,N,2);

[sorted_dist,dist_i] = sort(dist,2,'ascend');      % dist_i holds cell2 number (column), row indicated cell1

% for j = 2:(a_neigh+1)    % 2:(a_neigh+1) (to exclude cell itself)
%    for cell = 1:N   
%       
%        if sorted_dist(cell,j)<r_n
%            
%            junc(cell,dist_i(cell,j),2) = ((k35(conc(cell,SigM))+k35(conc(dist_i(cell,j),SigM)))/2)*(conc(cell,X27))*(conc(dist_i(cell,j),X27));
%            
%        end
%        
%    end
% end

for j = 2:(a_neigh+1)    % 2:(a_neigh+1) (to exclude cell itself) 
   
   ind_n = find(sorted_dist(:,j) < r_n);    % Indices of cells within r_n (to loop through less cells).
    
   for cell = 1:length(ind_n)
        junc(ind_n(cell),dist_i(ind_n(cell),j),2) = ( (k24(conc(ind_n(cell),SigM)) + k24(conc(dist_i(ind_n(cell),j),SigM)))/2 ).*(conc(ind_n(cell),X18)).*(conc(dist_i(ind_n(cell),j),X18));
   end
       
end
