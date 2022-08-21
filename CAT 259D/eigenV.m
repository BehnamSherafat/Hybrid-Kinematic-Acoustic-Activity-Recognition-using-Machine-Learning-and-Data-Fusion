format long g

latent = load('latent');
latent = latent.latent;
explained = load('explained');
explained = explained.explained;

figure(5);
subplot(2,1,1)
component_num = 1:15;
EigenValue = latent(1:15,:);
plot(component_num, EigenValue,'-o')
title('Components Eigenvalues','FontSize',20,'FontWeight','bold');
y = ylabel('Eigenvalue','FontSize',20,'FontWeight','bold');
ylim([0,10]);
x1 = xlabel('Principal Component','FontSize',26,'FontWeight','bold');
xlim([1,15]);
xticks(1:15)
hold on
yline(1,'-.b');
set(gca,'FontSize',20)


subplot(2,1,2)
paretoNew(explained)
title('Variance Explained by Components','FontSize',20,'FontWeight','bold');
y1 = ylabel('Variance Explained','FontSize',20,'FontWeight','bold');
xlabel('Principal Component','FontSize',26,'FontWeight','bold');
yticks(0:10:100)
set(gca,'FontSize',20)
x = 1:13;
y = explained(1:13,:);  
 for i1=1:numel(y)
    text(x(i1),y(i1),num2str(y(i1)),'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',18)
 end

