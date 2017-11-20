function [ dist_norm ] = DistCalc( calcPos,pSelectCalc,calcNparticles,config )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
             %% DISTANCES CALCULATION
                    del_r = calcPos - repmat(calcPos(pSelectCalc,:),[calcNparticles,1]);

                    del_r = abs (del_r);
                    del_rc = abs(repmat(config.SysLen,[calcNparticles,1]) - del_r) ;
                    del_r = min(del_r,del_rc) ;
            %         disp_norm_temp = del_r .* del_r ;
            %         disp_norm = (disp_norm_temp(:,1) + disp_norm_temp(:,2) + disp_norm_temp(:,3)).^.5 ;
                    del_r(pSelectCalc,:) = [];
                    del_r = del_r' ;
                    dist_norm = dot(del_r, del_r).^0.5;

end

