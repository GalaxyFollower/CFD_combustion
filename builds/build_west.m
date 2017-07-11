

% Nomenclature:
%
%     N(i-1,j) -  Ne     NE(i-1,j+1)
%
%       |         |         |
%
%      n ------  ne - - - nE
%                 |                 |
%       |         |        |        |       |
%                 |                 |
%   P (i,j) - -   e - -  E (i,j+1)
%                 |                 |
%       |         |        |        |       |
%                 |                 |
%       s ------ se - - - sE
%
%       |         |        |        |
%
%     S(i+1,j)  - Se      SE(i+1,j+1)
%
% Indexing of stecil: 

%   D_1 - D2 
%     |     | 
%   D_0 - D3
%     |     | 
%    D1 - D4

% Stecil 

% East 
D3=((dx_ne_n*(dx_P_e/4 + dx_e_Ne/4)*(C_E/8 + (3*C_N)/8 + C_NE/8 + (3*C_P)/8)^2)/S_nepsilon + (dx_s_se*(dx_e_P/4 + dx_Se_e/4)*(C_E/8 + (3*C_P)/8 + (3*C_S)/8 + C_SE/8)^2)/S_sepsilon + (dy_ne_n*(dy_P_e/4 + dy_e_Ne/4)*(C_E/8 + (3*C_N)/8 + C_NE/8 + (3*C_P)/8)^2)/S_nepsilon + (dy_s_se*(dy_e_P/4 + dy_Se_e/4)*(C_E/8 + (3*C_P)/8 + (3*C_S)/8 + C_SE/8)^2)/S_sepsilon - (dx_se_ne*(C_E/2 + C_P/2)^2*(dx_nE_n/4 + dx_s_sE/4 + dx_sE_nE))/S_e + (dy_se_ne*(C_E/2 + C_P/2)^2*(dy_nE_n/4 + dy_s_sE/4 + dy_sE_nE))/S_e)/S_epsilon; 

% West 
D_3=0; 

% South 
D1=((dx_s_se*(dx_P_S/2 + (3*dx_S_Se)/4 + dx_Se_e/4)*(C_E/8 + (3*C_P)/8 + (3*C_S)/8 + C_SE/8)^2)/S_sepsilon + (dy_s_se*(dy_P_S/2 + (3*dy_S_Se)/4 + dy_Se_e/4)*(C_E/8 + (3*C_P)/8 + (3*C_S)/8 + C_SE/8)^2)/S_sepsilon - (dx_s_sE*dx_se_ne*(C_E/2 + C_P/2)^2)/(4*S_e) + (dy_s_sE*dy_se_ne*(C_E/2 + C_P/2)^2)/(4*S_e))/S_epsilon; 

% North 
D_1=((dx_ne_n*(dx_N_P/2 + (3*dx_Ne_N)/4 + dx_e_Ne/4)*(C_E/8 + (3*C_N)/8 + C_NE/8 + (3*C_P)/8)^2)/S_nepsilon + (dy_ne_n*(dy_N_P/2 + (3*dy_Ne_N)/4 + dy_e_Ne/4)*(C_E/8 + (3*C_N)/8 + C_NE/8 + (3*C_P)/8)^2)/S_nepsilon - (dx_nE_n*dx_se_ne*(C_E/2 + C_P/2)^2)/(4*S_e) + (dy_nE_n*dy_se_ne*(C_E/2 + C_P/2)^2)/(4*S_e))/S_epsilon; 

% NW 
D_4=0; 

% NE 
D2=((dx_ne_n*(dx_Ne_N/4 + dx_e_Ne/4)*(C_E/8 + (3*C_N)/8 + C_NE/8 + (3*C_P)/8)^2)/S_nepsilon + (dy_ne_n*(dy_Ne_N/4 + dy_e_Ne/4)*(C_E/8 + (3*C_N)/8 + C_NE/8 + (3*C_P)/8)^2)/S_nepsilon - (dx_nE_n*dx_se_ne*(C_E/2 + C_P/2)^2)/(4*S_e) + (dy_nE_n*dy_se_ne*(C_E/2 + C_P/2)^2)/(4*S_e))/S_epsilon; 

% SW 
D_2=0; 

% SE 
D4=((dx_s_se*(dx_S_Se/4 + dx_Se_e/4)*(C_E/8 + (3*C_P)/8 + (3*C_S)/8 + C_SE/8)^2)/S_sepsilon + (dy_s_se*(dy_S_Se/4 + dy_Se_e/4)*(C_E/8 + (3*C_P)/8 + (3*C_S)/8 + C_SE/8)^2)/S_sepsilon - (dx_s_sE*dx_se_ne*(C_E/2 + C_P/2)^2)/(4*S_e) + (dy_s_sE*dy_se_ne*(C_E/2 + C_P/2)^2)/(4*S_e))/S_epsilon; 

% P 
D0=((dx_ne_n*(dx_N_P/2 + (3*dx_P_e)/4 + dx_e_Ne/4)*(C_E/8 + (3*C_N)/8 + C_NE/8 + (3*C_P)/8)^2)/S_nepsilon + (dx_s_se*(dx_P_S/2 + (3*dx_e_P)/4 + dx_Se_e/4)*(C_E/8 + (3*C_P)/8 + (3*C_S)/8 + C_SE/8)^2)/S_sepsilon + (dy_ne_n*(dy_N_P/2 + (3*dy_P_e)/4 + dy_e_Ne/4)*(C_E/8 + (3*C_N)/8 + C_NE/8 + (3*C_P)/8)^2)/S_nepsilon + (dy_s_se*(dy_P_S/2 + (3*dy_e_P)/4 + dy_Se_e/4)*(C_E/8 + (3*C_P)/8 + (3*C_S)/8 + C_SE/8)^2)/S_sepsilon - (dx_se_ne*(C_E/2 + C_P/2)^2*(dx_n_s + dx_nE_n/4 + dx_s_sE/4))/S_e + (dy_se_ne*(C_E/2 + C_P/2)^2*(dy_n_s + dy_nE_n/4 + dy_s_sE/4))/S_e)/S_epsilon; 

