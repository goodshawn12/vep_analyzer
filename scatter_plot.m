resp = load('average_resp.mat');
acc = load('average_acc.mat');

resp = resp.average_resp;
acc = acc.average_accuracy;

figure()
for stimuli = 1:4
    current_resp = zeros(1,30);
    current_acc = zeros(1,30);
    for i = 1:3
        for j = 1:10
            current_resp((i-1)*10+j) = resp(stimuli,i,j);
            current_acc((i-1)*10+j) = acc(stimuli,i,j);
        end
    end
    
    scatter(current_resp,current_acc,'filled')
    hold on
end

xlabel('Perceived intensity rating','FontSize',16)
ylabel('Cross validation accuracy (%)','FontSize',16)
title('Accuracy vs. Intensity Rating (grouped by stimuli)','FontSize',14)
legend('FMSEQ','MSEQ','LF SSVEP','HF SSVEP','location','northwest')
ylim([0,101]);
grid on

figure()
for contrast = 1:3
    current_resp = zeros(1,40);
    current_acc = zeros(1,40);
    for m = 1:4
        for n = 1:10
            current_resp((m-1)*10+n) = resp(m,contrast,n);
            current_acc((m-1)*10+n) = acc(m,contrast,n);
        end
    end
    
    scatter(current_resp,current_acc,'filled')
    hold on
end

xlabel('Perceived intensity rating','FontSize',16)
ylabel('Cross validation accuracy (%)','FontSize',16)
title('Accuracy vs. Intensity Rating (grouped by contrast)','FontSize',14)
legend('Contrast=2','Contrast=4','Contrast=8','location','northwest')
ylim([0,101]);
grid on


shape_list = ['o','o','o','o','o','o','o','^','^','^'];

figure()
for subject = 1:10
    current_resp = zeros(1,12);
    current_acc = zeros(1,12);
    for x = 1:4
        for y = 1:3
            current_resp((x-1)*3+y) = resp(x, y, subject);
            current_acc((x-1)*3+y) = acc(x, y, subject);
            if resp(x, y, subject) == 0
                fprintf(x)
                fprintf(y)
                fprintf(subject)
            end
        end
    end
    scatter(current_resp, current_acc, 'filled', shape_list(subject))
    hold on
end

xlabel('Perceived intensity rating','FontSize',16)
ylabel('Cross validation accuracy (%)','FontSize',16)
title('Accuracy vs. Intensity Rating (grouped by subject)','FontSize',14)
ylim([0,101])
legend('subject 1','subject 2','subject 3','subject 4','subject 5','subject 6','subject 7','subject 8', 'subject 9', 'subject 10', 'location', 'eastoutside')
grid on
