function plot_HJ(Mat, sort, Titles, spy_flag,text_flag)

    % sort = [{INT1},{MB1},{MB2},{INT2},{PT1},{PT2}];
    
    INT1 = sort{1};
    MB1 = sort{2};
    MB2 = sort{3};
    INT2 = sort{4};
    PT1 = sort{5};
    PT2 = sort{6};
    
    n1 = numel(PT1);
    n2 = numel(PT2);
    n = n1+n2;
    
    for i =1:numel(Mat)
    
        M = Mat{i};
        titl = Titles{i};
    
        figure;
        imagesc(M);colormap(bluewhitered(32));colorbar;
    
        if exist('spy_flag','var') && spy_flag == 1
            spy (M);
        end
    
        title(titl)
        % figure;imagesc(M-diag(diag(M)));colormap(bluewhitered(32));colorbar;title('J sorted')
        
        hold on; plot([0,n1+numel(MB2)+0.5],[numel(INT1)+0.5,numel(INT1)+0.5],'--k')
        hold on; plot([numel(INT1)+0.5,numel(INT1)+0.5],[0,n1+numel(MB2)+0.5],'--k')
        
        hold on; plot([0,n+0.5],[n1+0.5,n1+0.5],'-k','linewidth',2)
        hold on; plot([n1+0.5,n1+0.5],[0,n+0.5],'-k','linewidth',2)
        
        hold on; plot([numel(INT1)+0.5,n+0.5],[n1+numel(MB2)+0.5,n1+numel(MB2)+0.5],'--k')
        hold on; plot([n1+numel(MB2)+0.5,n1+numel(MB2)+0.5],[numel(INT1)+0.5,n+0.5],'--k')
        
        xlabel('Region/Node')
        ylabel('Region/Node')
        ax = gca;ax.FontWeight = 'bold';
    
        if exist('text_flag','var') && text_flag ==1
            text(numel(INT1)/2 - 2,numel(INT1)/2,'INT1','fontweight','bold','FontSize', 14,'Color','b')
            text(numel(INT1)+numel(MB1)/2 - 2,numel(INT1)+numel(MB1)/2,'B1','fontweight','bold','FontSize', 14,'Color',[0 0.7 0])
            text(n1 + numel(MB2)/2 - 2,n1 + numel(MB2)/2,'B2','fontweight','bold','FontSize', 14,'Color',[0 0.7 0])
            text(n1+n2-numel(INT2)/2 - 3,n1+n2-numel(INT2)/2,'INT2','fontweight','bold','FontSize', 14,'Color','b')
            text(numel(INT1)+numel(MB1)/2 - 6,n1 + numel(MB2)/2,'Cross-PT','fontweight','bold','FontSize', 12)
            text(n1 + numel(MB2)/2 - 6,numel(INT1)+numel(MB1)/2,'Cross-PT','fontweight','bold','FontSize', 12)
            text(numel(INT1)/2 ,n1+numel(MB2)+6,'MB-inducing','fontweight','bold','FontSize', 12,'Color',[0.8 0 0])
            text(n1+numel(MB2)-10,numel(INT1)/2,'MB-inducing','fontweight','bold','FontSize', 12,'Color',[0.8 0 0])
        end
        titl_str = split(convertCharsToStrings(titl),' (');
        saveas(gca,strcat(titl_str(1),'.png'))
    end
