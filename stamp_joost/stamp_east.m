function [stencil,l_s_n,S_gamma] = stamp_east(s, X_vec, Y_vec,dimX,dimY,lambda,alpha)

            % Nomenclature:
            %
            %    NW(s-dimY-1)      Nw     N(s-1)
            %                      |                   
            %                      |                   
            %                      |                   
            %            nW ------ nw ------ n 
            %            |         |         |                
            %            |         |         |                
            %            |         |         |
            %        W(s-dimY) - - w ----- P (s) 
            %            |         |         |
            %            |         |         |
            %            |         |         |
            %           sW ------ sw ------  s 
            %                      |                   
            %                      |                   
            %                      |                   
            %   SW(s-dimY+1)       Sw     S(s+1)
            %

            %% Calculation of Coordinates
            
            x_P  = X_vec(s);
            x_N  = X_vec(s-1);
            x_NW = X_vec(s-dimY-1);
            x_W  = X_vec(s-dimY);
            x_SW = X_vec(s-dimY+1);
            x_S  = X_vec(s+1);

            y_P  = Y_vec(s);
            y_N  = Y_vec(s-1);
            y_NW = Y_vec(s-dimY-1);
            y_W  = Y_vec(s-dimY);
            y_SW = Y_vec(s-dimY+1);
            y_S  = Y_vec(s+1);

            x_Nw = (x_NW+x_N)/2;
            x_Sw = (x_SW+x_S)/2;
            x_nW = (x_W+x_NW)/2;
            x_sW = (x_SW+x_W)/2;

            y_Nw = (y_NW+y_N)/2;
            y_Sw = (y_SW+y_S)/2;
            y_nW = (y_W+y_NW)/2;
            y_sW = (y_SW+y_W)/2;

            x_n  = (x_P+x_N)/2;
            x_w  = (x_W+x_P)/2;
            x_s  = (x_S+x_P)/2;
            x_nw = (x_w+x_Nw)/2;
            x_sw = (x_Sw+x_w)/2;

            y_n  = (y_P+y_N)/2;
            y_w  = (y_W+y_P)/2;
            y_s  = (y_S+y_P)/2;
            y_nw = (y_w+y_Nw)/2;
            y_sw = (y_Sw+y_w)/2;

            %% Calculation of Coordinate Differences
            
            %  Around gamma
            dx_sw_s  = x_s  - x_sw;
            dx_n_nw  = x_nw - x_n;
            dx_nw_sw = x_sw - x_nw;

            dy_sw_s  = y_s  - y_sw;
            dy_n_nw  = y_nw - y_n;
            dy_nw_sw = y_sw - y_nw;
            
            % Around ngamma
            dx_w_P   = x_P  - x_w;
            dx_P_N   = x_N  - x_P;
            dx_N_Nw  = x_Nw - x_N;
            dx_Nw_w  = x_w  - x_Nw;

            dy_w_P   = y_P  - y_w;
            dy_P_N   = y_N  - y_P;
            dy_N_Nw  = y_Nw - y_N;
            dy_Nw_w  = y_w  - y_Nw;
            
            %  Around w
            dx_s_n   = x_n  - x_s;   
            dx_n_nW  = x_nW - x_n;
            dx_nW_sW = x_sW - x_nW;
            dx_sW_s  = x_s  - x_sW;

            dy_s_n   = y_n  - y_s;   
            dy_n_nW  = y_nW - y_n;
            dy_nW_sW = y_sW - y_nW;
            dy_sW_s  = y_s  - y_sW;

            %  Around sgamma
            dx_S_P   = x_P  - x_S;
            dx_P_w   = x_w  - x_P;
            dx_w_Sw  = x_Sw - x_w;
            dx_Sw_S  = x_S  - x_Sw;

            dy_S_P   = y_P  - y_S;
            dy_P_w   = y_w  - y_P;
            dy_w_Sw  = y_Sw - y_w;
            dy_Sw_S  = y_S  - y_Sw;
            

            %% Computation of Control Surface Areas
            S_gamma  = abs((x_n*y_s - x_s*y_n) + (x_s*y_sw - x_sw*y_s)  + (x_sw*y_nw - x_nw*y_sw)  + (x_nw*y_n - x_n*y_nw)) / 2;
            S_ngamma = abs((x_N*y_P - x_P*y_N) + (x_P*y_w  - x_w*y_P)   + (x_w*y_Nw  - x_Nw*y_w)   + (x_Nw*y_N - x_N*y_Nw)) / 2;
            S_w      = abs((x_n*y_s - x_s*y_n) + (x_s*y_sW - x_sW*y_s)  + (x_sW*y_nW - x_nW*y_sW)  + (x_nW*y_n - x_n*y_nW)) / 2;
            S_sgamma = abs((x_P*y_S - x_S*y_P) + (x_S*y_Sw - x_Sw*y_S)  + (x_Sw*y_w  - x_w*y_Sw)   + (x_w*y_P  - x_P*y_w))  / 2;
           
            %% Initialize Stencil
            stencil = zeros(1,dimX*dimY);
                       
            %%  Build Stencil
            a_1 = (lambda*dy_n_nw)  / (S_gamma*S_ngamma);
            a_2 = (lambda*dx_n_nw)  / (S_gamma*S_ngamma);            
            a_3 = (lambda*dy_nw_sw) / (S_gamma*S_w);
            a_4 = (lambda*dx_nw_sw) / (S_gamma*S_w);            
            a_5 = (lambda*dy_sw_s)  / (S_gamma*S_sgamma);
            a_6 = (lambda*dx_sw_s)  / (S_gamma*S_sgamma);
            
            %  Computation of Node (Face) Length
            l_s_n = ((dx_s_n)^2+(dy_s_n)^2)^(0.5);

            %  Computation of Stencil Values
            %  P
            stencil(s)        = a_1*(dy_P_N/2+dy_Nw_w/4+dy_w_P*0.75)+a_2*(dx_P_N/2+dx_Nw_w/4+dx_w_P*0.75)+a_3*(dy_s_n+dy_n_nW/4+dy_sW_s/4)+a_4*(dx_s_n+dx_n_nW/4+dx_sW_s/4)+a_5*(dy_S_P/2+dy_P_w*0.75+dy_w_Sw/4)+a_6*(dx_S_P/2+dx_P_w*0.75+dx_w_Sw/4)-((alpha*l_s_n)/(S_gamma));
            %  South
            stencil (s+1)     = a_3*(dy_sW_s/4)+a_4*(dx_sW_s/4)+a_5*(dy_S_P/2+dy_w_Sw/4+dy_Sw_S*0.75)+a_6*(dx_S_P/2+dx_w_Sw/4+dx_Sw_S*0.75);
            %  North
            stencil (s-1)     = a_1*(dy_P_N/2+dy_N_Nw*0.75+dy_Nw_w/4)+a_2*(dx_P_N/2+dx_N_Nw*0.75+dx_Nw_w/4)+a_3*(dy_n_nW/4)+a_4*(dx_n_nW/4);  
            %  West
            stencil (s-dimY)  = a_1*(dy_Nw_w/4+dy_w_P*0.25)+a_2*(dx_Nw_w/4+dx_w_P*0.25)+a_3*(dy_n_nW/4+dy_nW_sW+dy_sW_s/4)+a_4*(dx_n_nW/4+dx_nW_sW+dx_sW_s/4)+a_5*(dy_P_w*0.25+dy_w_Sw/4)+a_6*(dx_P_w*0.25+dx_w_Sw/4);
            %  South-West
            stencil(s-dimY+1) = a_3*(dy_sW_s/4)+a_4*(dx_sW_s/4)+a_5*(dy_w_Sw/4+dy_Sw_S*0.25)+a_6*(dx_w_Sw/4+dx_Sw_S*0.25); 
            %  North-West
            stencil(s-dimY-1) = a_1*(dy_N_Nw*0.25+dy_Nw_w/4)+a_2*(dx_N_Nw*0.25+dx_Nw_w/4)+a_3*(dy_n_nW/4)+a_4*(dx_n_nW/4);
                                   
end



