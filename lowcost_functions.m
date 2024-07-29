classdef lowcost_functions
    methods     ( Static = true )
%         function  C = R(receivers, sources, sound_speed, delays)
%             %S, sources, sound_speed and delays are known.
%             Nr = size(receivers, 2);
%             Ns = size(sources, 2);
%             
%             % get distance along each axis
%             dx = (receivers(1,:).' - sources(1,:)); 
%             dy = (receivers(2,:).' - sources(2,:));
%             dz = (receivers(3,:).' - sources(3,:));
%         
%             % get travel times and delays
%             tau = sqrt(dx.^2 + dy.^2 + dz.^2) / sound_speed;
%             
%             delays_computed = NaN(Nr,Nr,Ns);
%             for i_s = 1:Ns
%                 delays_computed(:,:,i_s) = tau(:,i_s) - tau(:,i_s).';
%             end
%         
%             %% get cost function
%             C = sqrt(sum((delays_computed(:) - delays(:)).^2));
%         end

        function  cost = R(receivers, sources, sound_speed, delays)
            taus_x = tools.compute_taus_receivers_to_sources(receivers,sources,sound_speed);
            delays_x = tools.compute_receivers_delays(taus_x);

            cost = sqrt(sum((delays_x(:) - delays(:)).^2));
        end
    end
end
