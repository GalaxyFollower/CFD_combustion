function [stencil] = stamp_south(s, X_vec, Y_vec,dimX,dimY,lambda,alpha)

            % Nomenclature:
            %
            %    NW(s-dimY-1)      Nw     N(s-1)      Ne       NE(s+dimY-1)
            %                      |                   |
            %                      |                   |
            %                      |                   |
            %            nW ------ nw ------ n ------ ne - - - nE
            %            |         |         |         |        |
            %            |         |         |         |        |
            %            |         |         |         |        |
            %        W(s-dimY) - - w ----- P (s) ----- e - -  E (s+dimY)
            %

            %% Calculation of Coordinates
            
            x_P = X_vec(s);
            x_N = X_vec(s-1);
            x_NW = X_vec(s-dimY-1);
            x_W = X_vec(s-dimY);
            x_E = X_vec(s+dimY);
            x_NE = X_vec(s+dimY-1);

            y_P = Y_vec(s);
            y_N = Y_vec(s-1);
            y_NW = Y_vec(s-dimY-1);
            y_W = Y_vec(s-dimY);
            y_E = Y_vec(s+dimY);
            y_NE = Y_vec(s+dimY-1);

            x_Nw = (x_NW+x_N)/2;
            x_Ne = (x_N+x_NE)/2;
            x_nW = (x_W+x_NW)/2;
            x_nE = (x_E+x_NE)/2;

            y_Nw = (y_NW+y_N)/2;
            y_Ne = (y_N+y_NE)/2;
            y_nW = (y_W+y_NW)/2;
            y_nE = (y_E+y_NE)/2;

            x_n = (x_P+x_N)/2;
            x_w = (x_W+x_P)/2;
            x_e = (x_P+x_E)/2;
            x_nw = (x_w+x_Nw)/2;
            x_ne = (x_e+x_Ne)/2;

            y_n = (y_P+y_N)/2;
            y_w = (y_W+y_P)/2;
            y_e = (y_P+y_E)/2;
            y_nw = (y_w+y_Nw)/2;
            y_ne = (y_e+y_Ne)/2;
           
            %% Calculation of Coordinate Differences
            
            %  Around eta
            dx_e_ne = x_ne-x_e;
            dx_ne_nw = x_nw-x_ne;
            dx_nw_w = x_w-x_nw;

            dy_e_ne = y_ne-y_e;
            dy_ne_nw = y_nw-y_ne;
            dy_nw_w = y_w-y_nw;
            
            %  Around etae
            dx_P_E = x_E-x_P;
            dx_E_nE = x_nE-x_E;
            dx_nE_n = x_n-x_nE;
            dx_n_P = x_P-x_n;

            dy_P_E = y_E-y_P;
            dy_E_nE = y_nE-y_E;
            dy_nE_n = y_n-y_nE;
            dy_n_P = y_P-y_n;
            
            %  Around n
            dx_w_e = x_e-x_w;
            dx_e_Ne = x_Ne-x_e;
            dx_Ne_Nw = x_Nw-x_Ne;
            dx_Nw_w = x_w-x_Nw;

            dy_w_e = y_e-y_w;
            dy_e_Ne = y_Ne-y_e;
            dy_Ne_Nw = y_Nw-y_Ne;
            dy_Nw_w = y_w-y_Nw;
            
            % Around etaw
            dx_W_P = x_P-x_W;
            dx_P_n = x_n-x_P;
            dx_n_nW = x_nW-x_n;
            dx_nW_W = x_W-x_nW;
            
            dy_W_P = y_P-y_W;
            dy_P_n = y_n-y_P;
            dy_n_nW = y_nW-y_n;
            dy_nW_W = y_W-y_nW;


            %% Computation of Control Surface Areas
            S_eta  = abs((x_ne*y_e-x_e*y_ne) + (x_e*y_w-x_w*y_e) + (x_w*y_nw-x_nw*y_w) + (x_nw*y_ne-x_ne*y_nw)) / 2;
            S_etae = abs((x_nE*y_E-x_E*y_nE) + (x_E*y_P-x_P*y_E) + (x_P*y_n-x_n*y_P)   + (x_n*y_nE-x_nE*y_n))   / 2;        
            S_n    = abs((x_Ne*y_e-x_e*y_Ne) + (x_e*y_w-x_w*y_e) + (x_w*y_Nw-x_Nw*y_w) + (x_Nw*y_Ne-x_Ne*y_Nw)) / 2;
            S_etaw = abs((x_n*y_P-x_P*y_n)   + (x_P*y_W-x_W*y_P) + (x_W*y_nW-x_nW*y_W) + (x_nW*y_n-x_n*y_nW))   / 2;


            %% Initialize Stencil
            stencil = zeros(1,dimX*dimY);
                       
            %%  Build Stencil
            a_1 = (lambda*dy_e_ne)  / (S_eta*S_etae);
            a_2 = (lambda*dx_e_ne)  / (S_eta*S_etae);
            a_3 = (lambda*dy_ne_nw) / (S_eta*S_n);
            a_4 = (lambda*dx_ne_nw) / (S_eta*S_n);
            a_5 = (lambda*dy_nw_w)  / (S_eta*S_etaw);
            a_6 = (lambda*dx_nw_w)  / (S_eta*S_etaw);

            %  Computation of Stencil Values
            %  P
            stencil(s) = a_1*(dy_P_E/2+dy_nE_n/4+dy_n_P*0.75)+a_2*(dx_P_E/2+dx_nE_n/4+dx_n_P*0.75)+a_3*(dy_w_e+dy_e_Ne/4+dy_Nw_w/4)+a_4*(dx_w_e+dx_e_Ne/4+dx_Nw_w/4)+a_5*(dy_W_P/2+dy_P_n*0.75+dy_n_nW/4)+a_6*(dx_W_P/2+dx_P_n*0.75+dx_n_nW/4);            
            %  North
            stencil (s-1) = a_1*(dy_nE_n/4+dy_n_P/4)+a_2*(dx_nE_n/4+dx_n_P/4)+a_3*(dy_e_Ne/4+dy_Ne_Nw+dy_Nw_w/4)+a_4*(dx_e_Ne/4+dx_Ne_Nw+dx_Nw_w/4)+a_5*(dy_P_n/4+dy_n_nW/4)+a_6*(dx_P_n/4+dx_n_nW/4);            
            %  West
            stencil (s-dimY) = a_3*(dy_Nw_w/4)+a_4*(dx_Nw_w/4)+a_5*(dy_W_P/2+dy_n_nW/4+dy_nW_W*0.75)+a_6*(dx_W_P/2+dx_n_nW/4+dx_nW_W*0.75);            
            %  East
            stencil (s+dimY) = a_1*(dy_P_E/2+dy_E_nE*0.75+dy_nE_n/4)+a_2*(dx_P_E/2+dx_E_nE*0.75+dx_nE_n/4)+a_3*(dy_e_Ne/4)+a_4*(dx_e_Ne/4);
            %  North-West
            stencil(s-dimY-1) = a_3*(dy_Nw_w/4)+a_4*(dx_Nw_w/4)+a_5*(dy_n_nW/4+dy_nW_W/4)+a_6*(dx_n_nW/4+dx_nW_W/4);      
            %  North-East
            stencil(s+dimY-1) = a_1*(dy_E_nE/4+dy_nE_n/4)+a_2*(dx_E_nE/4+dx_nE_n/4)+a_3*(dy_e_Ne/4)+a_4*(dx_e_Ne/4);
            
end



