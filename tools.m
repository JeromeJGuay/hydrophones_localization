classdef tools
    methods ( Static = true )
        function taus = compute_taus_receivers_to_sources(receivers, sources, sound_speed)
            % Compute the travel times (taus) between receivers and sources.
            %
            % Parameters
            % ----------
            %   receivers: Nr x 3, [rx1 ry1 rz1; ... ; rxN ryN rzN]
            %   sources: Ns x 3, [sx1 sy1 sz1; ... ; sxN syN szN]
            %   sound_speed: Scalar.
            %
            % Returns
            % -------
            %   distances (Nr x Ns)
            taus = compute_distances_sources_to_receivers(receivers, sources) / sound_speed;
        end

        function distances = compute_distances_receivers_to_sources(receivers, sources)
            % Compute the distance between receivers and sources.
            %
            % Parameters
            % ----------
            %   receivers: Nr x 3, [rx1 ry1 rz1; ... ; rxN ryN rzN]
            %   sources: Ns x 3, [sx1 sy1 sz1; ... ; sxN syN szN]
            %
            % Returns
            % -------
            %   distances (Nr x Ns)
            Nr = length(receivers); % Ns x 3
            Ns = length(sources); %  Nr x 3
            
            
            s_expanded = reshape(sources, [1, 3, Ns]); % 3 x Ns x 1
            r_expanded = reshape(receivers', [Nr, 3, 1]); % 3 x 1 x Nr
            
            % Compute pairwise differences (sources to receivers wise)
            diffs = r_expanded - s_expanded; % 3 x Ns x Nr
            
            % Compute squared distances
            squared_distances = sum(diffs.^2, 2); % 1 x Ns x Nr
            
            % Compute distances by taking the square root of squared distances
            distances = reshape(sqrt(squared_distances), [4,6]); % Ns x Nr
        end
    end
end

