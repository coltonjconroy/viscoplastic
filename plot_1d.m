subplot(2,2,1:2)
for j = 1:nelems
    x = [PSI.xa]*X(ELEM(j).nodes)';
    xplot = [  X(j); x ; X(j+1)];
%     SEp = zeros(length(xplot),1);
%     for i = 1:length(xplot)
%         [gradse,~,se,vavg] = problem_data(xplot(i),time,prob,type);
%         SEp(i) = se;
%     end
%     plot(xplot,SEp,'r:','LineWidth',5)
    hold on
    uplot = [ [PHI.b1]*h(:,j,irk); [PHI.L2]*h(:,j,irk); [PHI.b2]*h(:,j,irk) ];
    
    plot([xplot(1),xplot(end)],[uplot(1),uplot(end)],'bo',...
        'LineWidth',2,...
        'MarkerEdgeColor','b',...
        'MarkerFaceColor',[.49 1 .63],...
        'MarkerSize',10)
    plot(xplot,uplot,'b-')

end
xlabel(['x-coordinate (m)'],'FontSize',13,'FontWeight','demi');
ylabel('Debris thickness h_h (m)','FontSize',13,'FontWeight','demi')%,'FontName','SansSerif'
%legend('  Analytic Solution',['  DG, p = ',num2str(p)],'Location','Best')
%text(-0.90*umax,-4.5,['t = ',num2str(time/3600),' hrs'],'FontSize',13,'FontWeight','demi')
subplot(2,2,3:4)
for j = 1:nelems
    x = [PSI.xa]*X(ELEM(j).nodes)';
    xplot = [  X(j); x ; X(j+1)];
%     SEp = zeros(length(xplot),1);
%     for i = 1:length(xplot)
%         [gradse,~,se,vavg] = problem_data(xplot(i),time,prob,type);
%         SEp(i) = se;
%     end
%     plot(xplot,SEp,'r:','LineWidth',5)
    hold on
    uplot = [ [PHI.b1]*dhdx(:,j); [PHI.L2]*dhdx(:,j); [PHI.b2]*dhdx(:,j) ];
    
    plot([xplot(1),xplot(end)],[uplot(1),uplot(end)],'bo',...
        'LineWidth',2,...
        'MarkerEdgeColor','b',...
        'MarkerFaceColor',[.49 1 .63],...
        'MarkerSize',10)
    plot(xplot,uplot,'b-')

end
xlabel(['x-coordinate (m)'],'FontSize',13,'FontWeight','demi');
ylabel('dh/dx','FontSize',13,'FontWeight','demi')
