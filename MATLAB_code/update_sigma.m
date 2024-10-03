function [sigma_aj_out,conc] = update_sigma(sigma_aj_in,conc,junc)
% updates sigma (factor describing degree of N-glycosylation) in membrane 
% and AJ E-cadherin pools after
%
% February 2, 2015

X18 = 18;
SigM = 21;

junc(junc==0) = nan;
avg_junc = nanmean(junc,2);
avg_junc(isnan(avg_junc)) = 0;

C_nz = find(avg_junc~=0);

sigma_aj_out = sigma_aj_in;
delta_sigma_aj = ((k24(conc((C_nz),SigM)).*conc(C_nz,X18))./(avg_junc(C_nz))).*(conc(C_nz,SigM)-sigma_aj_in(C_nz));
sigma_aj_out(C_nz) = sigma_aj_in(C_nz) + delta_sigma_aj;

delta_sigma_M = ((k_24(sigma_aj_in(C_nz)).*avg_junc(C_nz))./(conc(C_nz,X18))).*(sigma_aj_in(C_nz)-conc(C_nz,SigM));
conc(C_nz,SigM) = conc(C_nz,SigM) + delta_sigma_M;
