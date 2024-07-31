classdef tools
    methods ( Static = true )
        
        function delays = compute_receivers_delays(taus)
            % Compute the delays between sources-receivers travel times
            %
            % Parameters
            % ----------
            %   taus: Nr x Ns
            %
            % Returns
            % -------
            %   delays: (Ns x Nr x Nr)
            Ns = size(taus, 2);
            Nr = size(taus, 1);

            delays = NaN(Ns,Nr,Nr);
            for i_e = 1:Ns
                delays(i_e,:,:) = taus(:, i_e) - taus(:, i_e).';
            end
        end

        function r = rand0(dim, norm)
            r = 2 * norm * (rand(dim) - .5);
        end
    end
end

