classdef lowcost_functions
    methods     ( Static = true )
       function  cost = R(X, sources, sound_speed, delays)
           % Unknown receivers
            receivers = reshape(X, length(X)/3, 3);
            taus_x = pdist2(receivers, sources) / sound_speed;
            delays_x = tools.compute_receivers_delays(taus_x);

            cost = sqrt(sum((delays_x(:) - delays(:)).^2));
       end

        function  cost = RC(X, sources, delays)
            % Unknown Receivers and sound speed
            rX = X(1:end-1);
            sound_speed = X(end);

            receivers = reshape(rX, length(rX)/3, 3);
            taus_x = pdist2(receivers, sources) / sound_speed;
            delays_x = tools.compute_receivers_delays(taus_x);

            cost = sqrt(sum((delays_x(:) - delays(:)).^2));
        end
    end
end
