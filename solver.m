classdef solver
    methods     ( Static = true )
        
        function xsol = R(r0, sources, sound_speed, dtaus)
            X0 = r0(:);

            fmincon_options = optimoptions('fmincon','display','off');
            
            cfun = @(X) lowcost_functions.R(X, sources, sound_speed, dtaus);
            
            X1 = fmincon(cfun, X0,  [],[],[],[],[],[],[], fmincon_options); %

            xsol = reshape(X1, length(X0)/3, 3);
        end
        
    end
end

