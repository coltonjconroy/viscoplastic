if  plot_n == n
    %             time;
    figure('Renderer','zbuffer');
    set(gca,'NextPlot','replaceChildren'); hold on;  axis([x0, xN, -umax-.2, umax+.2])
    set(gcf,'position',get(0,'screensize'))
    set(gca,'FontSize',13,'FontWeight','demi')
    box on
    subplot(2,1,1)
    for j = 1:nelems
        x = [PSI.xa]*X(ELEM(j).nodes)';
        xplot = [  X(j); x ; X(j+1)];
        %         SEp = zeros(length(xplot),1);
        %         for i = 1:length(xplot)
        %             [gradse,~,se,vavg] = problem_data(xplot(i),time,prob,type);
        %             SEp(i) = se;
        %         end
        %         plot(xplot,SEp,'r:','LineWidth',5)
        hold on
        uplot = [ [PHI.b1]*h(:,j,1); [PHI.L2]*h(:,j,1); [PHI.b2]*h(:,j,1) ];
        
        plot([xplot(1),xplot(end)],[uplot(1),uplot(end)],'bo',...
            'LineWidth',2,...
            'MarkerEdgeColor','b',...
            'MarkerFaceColor',[.49 1 .63],...
            'MarkerSize',10)
        plot(xplot,uplot,'b-')
        
    end
    xlabel(['x-coordinate (m)'],'FontSize',13,'FontWeight','demi');
    ylabel('Thickness h_h (m)','FontSize',13,'FontWeight','demi')%,'FontName','SansSerif'
    axis([x0 xN -.25 1.6])
    % legend('  Analytic Solution',['  DG, p = ',num2str(p)],'Location','NorthWest')
    % text(-0.90*umax,-4.5,['t = ',num2str(time/3600),' hrs'],'FontSize',13,'FontWeight','demi')
    subplot(2,1,2)
    for j = 1:nelems
        x = [PSI.xa]*X(ELEM(j).nodes)';
        xplot = [  X(j); x ; X(j+1)];
        %         vp = zeros(length(xplot),1);
        %         for i = 1:length(xplot)
        %             [~,~,~,vavg] = problem_data(xplot(i),time,prob,type);
        %             vp(i) = vavg;
        %         end
        %         plot(xplot,vp,'r:','LineWidth',5)
        hold on
        
        if strcmp(u_fig,'on')
            uplot = [ [PHI.b1]*hu(:,j,1); [PHI.L2]*hu(:,j,1); [PHI.b2]*hu(:,j,1) ];
        elseif strcmp(p_fig,'on')
            uplot = [ [PHI.b1]*pb(:,j); [PHI.L2]*pb(:,j); [PHI.b2]*pb(:,j) ];
        end
        
        
        plot([xplot(1),xplot(end)],[uplot(1),uplot(end)],'bo',...
            'LineWidth',2,...
            'MarkerEdgeColor','b',...
            'MarkerFaceColor',[.49 1 .63],...
            'MarkerSize',10)
        plot(xplot,uplot,'b-')
        
    end
    xlabel(['x-coordinate (m)'],'FontSize',13,'FontWeight','demi');
    ylabel('Depth Avg. momentum (m/s)','FontSize',13,'FontWeight','demi')
    axis([x0 xN -1.5 10])
    %legend('  Analytic Solution',['  DG, p = ',num2str(p)],'Location','SouthWest')
    % text(-0.90*umax,-4.5,['t = ',num2str(time/3600),' hrs'],'FontSize',13,'FontWeight','demi')
    frame = getframe(gcf);
    writeVideo(writerObj,frame);
    plot_n = plot_n + floor(.25/dt);
    close
end