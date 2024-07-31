classdef lowcost_functions
    methods     ( Static = true )
       function  cost = R(X, sources, sound_speed, delays)
           % Cost function for unknown receivers locations.
           % 
           % Parameters
           % --------- 
           % X : 1-D array size (3*Nr) 
           %   Receivers coords: [x1 y1 z1 ... xNr yNr zNr] 
           % sources : 2-D array size Ns x 3 
           %   Sources coords: [x1 y2 y3; ... ; xNs yNs zNs]
           % sound_speed: float
           %
           % delays : 3-D array size (Ns x Nr x Nr)
           %
           % Returns
           % -------
           %   sqrt( sum((delays_X - delays)^2) )

           receivers = reshape(X, length(X)/3, 3);
           taus_x = pdist2(receivers, sources) / sound_speed;
           delays_x = tools.compute_receivers_delays(taus_x);

           cost = sqrt(sum((delays_x(:) - delays(:)).^2));
       end

%         function  cost = RC(X, sources, delays)
%             % Unknown Receivers and sound speed
%             rX = X(1:end-1);
%             sound_speed = X(end);
% 
%             receivers = reshape(rX, length(rX)/3, 3);
%             taus_x = pdist2(receivers, sources) / sound_speed;
%             delays_x = tools.compute_receivers_delays(taus_x);
% 
%             cost = sqrt(sum((delays_x(:) - delays(:)).^2));
%         end
    end
end
