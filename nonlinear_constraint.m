classdef nonlinear_constraint
    methods     ( Static = true )
        function [c,ceq] = max_distance(x, distance)
            c(1) = sum((x(:,1) - x(:,2)).^2).^(1/2) - distance;
            c(2) = sum((x(:,1) - x(:,3)).^2).^(1/2) - distance;
            c(3) = sum((x(:,1) - x(:,4)).^2).^(1/2) - distance;
            c(4) = sum((x(:,2) - x(:,3)).^2).^(1/2) - distance;
            c(5) = sum((x(:,2) - x(:,4)).^2).^(1/2) - distance;
            c(6) = sum((x(:,3) - x(:,4)).^2).^(1/2) - distance;
            ceq = [];
        end
    end
end
