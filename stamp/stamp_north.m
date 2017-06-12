function [stencil,l_e_w,S_sigma] = stamp_north(s, X_vec, Y_vec,dimX,dimY,lambda,alpha)

            % Nomenclature:
            %
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
            x_W  = X_vec(s-dimY);
            x_SW = X_vec(s-dimY+1);
            x_S  = X_vec(s+1);
            x_SE = X_vec(s+dimY+1);
            x_E  = X_vec(s+dimY);

            y_P  = Y_vec(s);
            y_W  = Y_vec(s-dimY);
            y_SW = Y_vec(s-dimY+1);
            y_S  = Y_vec(s+1);
            y_SE = Y_vec(s+dimY+1);
            y_E  = Y_vec(s+dimY);

            x_Sw = (x_SW+x_S)/2;
            x_Se = (x_S+x_SE)/2;
            x_sW = (x_SW+x_W)/2;
            x_sE = (x_SE+x_E)/2;

            y_Sw = (y_SW+y_S)/2;
            y_Se = (y_S+y_SE)/2;
            y_sW = (y_SW+y_W)/2;
            y_sE = (y_SE+y_E)/2;

            x_w  = (x_W+x_P)/2;
            x_s  = (x_S+x_P)/2;
            x_e  = (x_P+x_E)/2;
            x_sw = (x_Sw+x_w)/2;
            x_se = (x_Se+x_e)/2;

            y_w  = (y_W+y_P)/2;
            y_s  = (y_S+y_P)/2;
            y_e  = (y_P+y_E)/2;
            y_sw = (y_Sw+y_w)/2;
            y_se = (y_Se+y_e)/2;
                    
            %% Calculation of Coordinate Differences
            
            %  Around sigma
            dx_sw_se = x_se - x_sw;
            dx_se_e  = x_e  - x_se;
            dx_w_sw  = x_sw - x_w;
            
            dy_sw_se = y_se - y_sw;
            dy_se_e  = y_e  - y_se;
            dy_w_sw  = y_sw - y_w;
            
            %  Around sigmaw
            dx_sW_s  = x_s  - x_sW;
            dx_s_P   = x_P  - x_s;
            dx_P_W   = x_W  - x_P;
            dx_W_sW  = x_sW - x_W;
            
            dy_sW_s  = y_s  - y_sW;
            dy_s_P   = y_P  - y_s;
            dy_P_W   = y_W  - y_P;
            dy_W_sW  = y_sW - y_W;
            
            %  Around s
            dx_Sw_Se = x_Se - x_Sw;
            dx_Se_e  = x_e  - x_Se;
            dx_e_w   = x_w  - x_e;
            dx_w_Sw  = x_Sw - x_w;
            
            dy_Sw_Se = y_Se - y_Sw;
            dy_Se_e  = y_e  - y_Se;
            dy_e_w   = y_w  - y_e;
            dy_w_Sw  = y_Sw - y_w;
            
            %  Around sigmae
            dx_s_sE = x_sE - x_s;
            dx_sE_E = x_E  - x_sE;
            dx_E_P  = x_P  - x_E;
            dx_P_s  = x_s  - x_P;
            
            dy_s_sE = y_sE - y_s;
            dy_sE_E = y_E  - y_sE;
            dy_E_P  = y_P  - y_E;
            dy_P_s  = y_s  - y_P;
            

            %% Computation of Control Surface Areas
            S_sigma  = abs((x_e*y_se - x_se*y_e) + (x_se*y_sw - x_sw*y_se) + (x_sw*y_w - x_w*y_sw) + (x_w*y_e - x_e*y_w)) / 2;
            S_sigmaw = abs((x_P*y_s  - x_s*y_P)  + (x_s*y_sW  - x_sW*y_s)  + (x_sW*y_W - x_W*y_sW) + (x_W*y_P - x_P*y_W)) / 2;
            S_s      = abs((x_e*y_Se - x_Se*y_e) + (x_Se*y_Sw - x_Sw*y_Se) + (x_Sw*y_w - x_w*y_Sw) + (x_w*y_e - x_e*y_w)) / 2;
            S_sigmae = abs((x_E*y_sE - x_sE*y_E) + (x_sE*y_s  - x_s*y_sE)  + (x_s*y_P  - x_P*y_s)  + (x_P*y_E - x_E*y_P)) / 2;
            

            %% Initialize Stencil
            stencil = zeros(1,dimX*dimY);
                       
            %%  Build Stencil
            a_1 = (lambda*dy_w_sw)  / (S_sigma*S_sigmaw);
            a_2 = (lambda*dx_w_sw)  / (S_sigma*S_sigmaw);            
            a_3 = (lambda*dy_sw_se) / (S_sigma*S_s);
            a_4 = (lambda*dx_sw_se) / (S_sigma*S_s);            
            a_5 = (lambda*dy_se_e)  / (S_sigma*S_sigmae);
            a_6 = (lambda*dx_se_e)  / (S_sigma*S_sigmae);
            
            %  Computation of Node (Face) Length
            l_e_w = ((dx_e_w)^2+(dy_e_w)^2)^(0.5);

            %  Computation of Stencil Values
            %  P
            stencil(s)        = a_1*(dy_sW_s/4+dy_s_P*0.75+dy_P_W*0.5)+a_2*(dx_sW_s/4+dx_s_P*0.75+dx_P_W*0.5)+a_3*(dy_Se_e/4+dy_e_w+dy_w_Sw/4)+a_4*(dx_Se_e/4+dx_e_w+dx_w_Sw/4)+a_5*(dy_s_sE/4+dy_E_P*0.5+dy_P_s*0.75)+a_6*(dx_s_sE/4+dx_E_P*0.5+dx_P_s*0.75)-((alpha*l_e_w)/(S_sigma));
            %  South
            stencil (s+1)     = a_1*(dy_sW_s/4+dy_s_P*0.25)+a_2*(dx_sW_s/4+dx_s_P*0.25)+a_3*(dy_Sw_Se+dy_Se_e/4+dy_w_Sw/4)+a_4*(dx_Sw_Se+dx_Se_e/4+dx_w_Sw/4)+a_5*(dy_s_sE/4+dy_P_s*0.25)+a_6*(dx_s_sE/4+dx_P_s*0.25);
            %  West
            stencil (s-dimY)  = a_1*(dy_sW_s/4+dy_P_W*0.5+dy_W_sW*0.75)+a_2*(dx_sW_s/4+dx_P_W*0.5+dx_W_sW*0.75)+a_3*(dy_w_Sw/4)+a_4*(dx_w_Sw/4);           
            %  East
            stencil (s+dimY)  = a_3*(dy_Se_e/4)+a_4*(dx_Se_e/4)+a_5*(dy_s_sE/4+dy_sE_E*0.75+dy_E_P*0.5)+a_6*(dx_s_sE/4+dx_sE_E*0.75+dx_E_P*0.5);            
            %  South-West
            stencil(s-dimY+1) = a_1*(dy_sW_s/4+dy_W_sW*0.25)+a_2*(dx_sW_s/4+dx_W_sW*0.25)+a_3*(dy_w_Sw/4)+a_4*(dx_w_Sw/4);                      
            %  South-East
            stencil(s+dimY+1) = a_3*(dy_Se_e*0.25)+a_4*(dx_Se_e*0.25)+a_5*(dy_s_sE/4+dy_sE_E*0.25)+a_6*(dx_s_sE/4+dx_sE_E*0.25);
                      
end



