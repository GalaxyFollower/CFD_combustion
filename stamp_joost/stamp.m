function [stencil] = stamp(s, X_vec, Y_vec,dimX,dimY,lambda)

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
            %            |         |         |         |        |
            %            |         |         |         |        |
            %            |         |         |         |        |
            %           sW ------ sw ------  s ------ se - - - sE
            %                      |                   |
            %                      |                   |
            %                      |                   |
            %   SW(s-dimY+1)       Sw     S(s+1)      Se       SE(s+dimY+1)
            %

            %% Calculation of Coordinates

            x_P  = X_vec(s);
            x_N  = X_vec(s-1);
            x_NW = X_vec(s-dimY-1);
            x_W  = X_vec(s-dimY);
            x_SW = X_vec(s-dimY+1);
            x_S  = X_vec(s+1);
            x_SE = X_vec(s+dimY+1);
            x_E  = X_vec(s+dimY);
            x_NE = X_vec(s+dimY-1);

            y_P  = Y_vec(s);
            y_N  = Y_vec(s-1);
            y_NW = Y_vec(s-dimY-1);
            y_W  = Y_vec(s-dimY);
            y_SW = Y_vec(s-dimY+1);
            y_S  = Y_vec(s+1);
            y_SE = Y_vec(s+dimY+1);
            y_E  = Y_vec(s+dimY);
            y_NE = Y_vec(s+dimY-1);

            x_Nw = (x_NW+x_N)/2;
            x_Ne = (x_N+x_NE)/2;
            x_Sw = (x_SW+x_S)/2;
            x_Se = (x_S+x_SE)/2;
            x_nW = (x_W+x_NW)/2;
            x_sW = (x_SW+x_W)/2;
            x_nE = (x_E+x_NE)/2;
            x_sE = (x_SE+x_E)/2;

            y_Nw = (y_NW+y_N)/2;
            y_Ne = (y_N+y_NE)/2;
            y_Sw = (y_SW+y_S)/2;
            y_Se = (y_S+y_SE)/2;
            y_nW = (y_W+y_NW)/2;
            y_sW = (y_SW+y_W)/2;
            y_nE = (y_E+y_NE)/2;
            y_sE = (y_SE+y_E)/2;

            x_n  = (x_P+x_N)/2;
            x_w  = (x_W+x_P)/2;
            x_s  = (x_S+x_P)/2;
            x_e  = (x_P+x_E)/2;
            x_nw = (x_w+x_Nw)/2;
            x_sw = (x_Sw+x_w)/2;
            x_se = (x_Se+x_e)/2;
            x_ne = (x_e+x_Ne)/2;

            y_n  = (y_P+y_N)/2;
            y_w  = (y_W+y_P)/2;
            y_s  = (y_S+y_P)/2;
            y_e  = (y_P+y_E)/2;
            y_nw = (y_w+y_Nw)/2;
            y_sw = (y_Sw+y_w)/2;
            y_se = (y_Se+y_e)/2;
            y_ne = (y_e+y_Ne)/2;


            %% Calculation of Coordinate Differences
            
            %  Around P
            dx_sw_se = x_se - x_sw;
            dx_se_ne = x_ne - x_se;
            dx_ne_nw = x_nw - x_ne;
            dx_nw_sw = x_sw - x_nw;

            dy_sw_se = y_se - y_sw;
            dy_se_ne = y_ne - y_se;
            dy_ne_nw = y_nw - y_ne;
            dy_nw_sw = y_sw - y_nw;

            %  Around s
            dx_Sw_Se = x_Se - x_Sw;
            dx_Se_e  = x_e  - x_Se;
            dx_e_w   = x_w  - x_e;
            dx_w_Sw  = x_Sw - x_w;

            dy_Sw_Se = y_Se - y_Sw;
            dy_Se_e  = y_e  - y_Se;
            dy_e_w   = y_w  - y_e;
            dy_w_Sw  = y_Sw - y_w;

            %  Around e
            dx_s_sE  = x_sE - x_s;
            dx_sE_nE = x_nE - x_sE;
            dx_nE_n  = x_n  - x_nE;
            dx_n_s   = x_s  - x_n;

            dy_s_sE  = y_sE - y_s;
            dy_sE_nE = y_nE - y_sE;
            dy_nE_n  = y_n  - y_nE;
            dy_n_s   = y_s  - y_n;

            %  Around n
            dx_w_e   = x_e  - x_w;
            dx_e_Ne  = x_Ne - x_e;
            dx_Ne_Nw = x_Nw - x_Ne;
            dx_Nw_w  = x_w  - x_Nw;

            dy_w_e   = y_e  - y_w;
            dy_e_Ne  = y_Ne - y_e;
            dy_Ne_Nw = y_Nw - y_Ne;
            dy_Nw_w  = y_w  - y_Nw;

            %  Around w
            dx_sW_s  = x_s  - x_sW;
            dx_s_n   = x_n  - x_s;
            dx_n_nW  = x_nW - x_n;
            dx_nW_sW = x_sW - x_nW;

            dy_sW_s  = y_s  - y_sW;
            dy_s_n   = y_n  - y_s;
            dy_n_nW  = y_nW - y_n;
            dy_nW_sW = y_sW - y_nW;



            %% Computation of Control Surface Areas
            S_P = abs((x_ne*y_se-x_se*y_ne) + (x_se*y_sw-x_sw*y_se) + (x_sw*y_nw-x_nw*y_sw) + (x_nw*y_ne-x_ne*y_nw)) / 2;
            S_s = abs((x_e*y_Se-x_Se*y_e)   + (x_Se*y_Sw-x_Sw*y_Se) + (x_Sw*y_w-x_w*y_Sw)   + (x_w*y_e-x_e*y_w))     / 2;
            S_e = abs((x_nE*y_sE-x_sE*y_nE) + (x_sE*y_s-x_s*y_sE)   + (x_s*y_n-x_n*y_s)     + (x_n*y_nE-x_nE*y_n))   / 2;
            S_n = abs((x_Ne*y_e-x_e*y_Ne)   + (x_e*y_w-x_w*y_e)     + (x_w*y_Nw-x_Nw*y_w)   + (x_Nw*y_Ne-x_Ne*y_Nw)) / 2;
            S_w = abs((x_n*y_s-x_s*y_n)     + (x_s*y_sW-x_sW*y_s)   + (x_sW*y_nW-x_nW*y_sW) + (x_nW*y_n-x_n*y_nW))   / 2;


            %% Initialize Stencil
            stencil = zeros(1,dimX*dimY);
            
            %%  Build Stencil
            a_1 = (lambda*dy_sw_se)/(S_P*S_s);
            a_2 = (lambda*dx_sw_se)/(S_P*S_s);
            a_3 = (lambda*dy_se_ne)/(S_P*S_e);
            a_4 = (lambda*dx_se_ne)/(S_P*S_e);
            a_5 = (lambda*dy_ne_nw)/(S_P*S_n);
            a_6 = (lambda*dx_ne_nw)/(S_P*S_n);
            a_7 = (lambda*dy_nw_sw)/(S_P*S_w);
            a_8 = (lambda*dx_nw_sw)/(S_P*S_w);

            %  Computation of Stencil Values
            %  P
            stencil(s) = a_1*(dy_e_w+dy_Se_e/4+dy_w_Sw/4)+a_2*(dx_e_w+dx_Se_e/4+dx_w_Sw/4)+a_3*(dy_n_s+dy_s_sE/4+dy_nE_n/4)+a_4*(dx_n_s+dx_s_sE/4+dx_nE_n/4)+a_5*(dy_w_e+dy_e_Ne/4+dy_Nw_w/4)+a_6*(dx_w_e+dx_e_Ne/4+dx_Nw_w/4)+a_7*(dy_s_n+dy_sW_s/4+dy_n_nW/4)+a_8*(dx_s_n+dx_sW_s/4+dx_n_nW/4);
            %  North
            stencil(s-1) = a_3*(dy_nE_n/4)+a_4*(dx_nE_n/4)+a_5*(dy_Ne_Nw+dy_e_Ne/4+dy_Nw_w/4)+a_6*(dx_Ne_Nw+dx_e_Ne/4+dx_Nw_w/4)+a_7*(dy_n_nW/4)+a_8*(dx_n_nW/4);
            %  West
            stencil(s-dimY) = a_1*(dy_w_Sw/4)+a_2*(dx_w_Sw/4)+a_5*(dy_Nw_w/4)+a_6*(dx_Nw_w/4)+a_7*(dy_nW_sW+dy_n_nW/4+dy_sW_s/4)+a_8*(dx_nW_sW+dx_sW_s/4+dx_n_nW/4);
            %  South
            stencil(s+1) = a_1*(dy_Sw_Se+dy_Se_e/4+dy_w_Sw/4)+a_2*(dx_Sw_Se+dx_Se_e/4+dx_w_Sw/4)+a_3*(dy_s_sE/4)+a_4*(dx_s_sE/4)+a_7*(dy_sW_s/4)+a_8*(dx_sW_s/4);
            %  East
            stencil(s+dimY) = a_1*(dy_Se_e/4)+a_2*(dx_Se_e/4)+a_3*(dy_s_sE/4+dy_sE_nE+dy_nE_n/4)+a_4*(dx_s_sE/4+dx_sE_nE+dx_nE_n/4)+a_5*(dy_e_Ne/4)+a_6*(dx_e_Ne/4);
            %  North-West
            stencil(s-dimY-1) = a_5*(dy_Nw_w/4)+a_6*(dx_Nw_w/4)+a_7*(dy_n_nW/4)+a_8*(dx_n_nW/4);
            %  South-West
            stencil(s-dimY+1) = a_1*(dy_w_Sw/4)+a_2*(dx_w_Sw/4)+a_7*(dy_sW_s/4)+a_8*(dx_sW_s/4);
            %  South-East
            stencil(s+dimY+1) = a_1*(dy_Se_e/4)+a_2*(dx_Se_e/4)+a_3*(dy_s_sE/4)+a_4*(dx_s_sE/4);
            %  North-East
            stencil(s+dimY-1) = a_3*(dy_nE_n/4)+a_4*(dx_nE_n/4)+a_5*(dy_e_Ne/4)+a_6*(dx_e_Ne/4);

end



