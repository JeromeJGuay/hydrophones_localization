classdef figures
    methods ( Static = true )
        function pyramid_3d_solution(real, guess, solution)
            % Plots 3D pyramids for the real, guess and solution.
            labels = 1:size(real, 2);
            figure
            hold on
           
            figures.pyramid_3d_plot(real, 'b', labels);
            figures.pyramid_3d_plot(guess,'k', labels);
            figures.pyramid_3d_plot(solution, 'r', labels);

            hold off
            set(gca,'ZDir','reverse')
            grid on
            legend('Guess','', 'Real', '', 'Solution')
            view(3)
            daspect([1 1 1])
            shg
        end
    end

    methods (Static= true, Access= private)
        
        function pyramid_3d_plot(receivers, color, labels)
            
            plot3( ...
                receivers(:, 1), ...
                receivers(:, 2), ...
                receivers(:, 3), ...
                strcat(color,'.'),'MarkerSize',10)
            
            trisurf( ...
                boundary(receivers, 1), ...
                receivers(:, 1), ...
                receivers(:, 2), ...
                receivers(:, 3), ...
                'FaceColor',color,'FaceAlpha',0.1)
            
            for i=1:size(receivers,2)
               text(receivers(1,i),receivers(2,i),receivers(3,i),['   ' ...
                num2str(labels(i))],'HorizontalAlignment','left','FontSize',8);
            end    
        end
    
    end
end

