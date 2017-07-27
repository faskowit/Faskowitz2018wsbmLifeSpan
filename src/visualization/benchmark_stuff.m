subplot(1,3,1)
scatter(numTrials_to_test,median(numTrialsTest_logE,2))
ylabel('log evidence, dots: median, line: mean')
xlabel('num trials')
hold on
plot(numTrials_to_test,mean(numTrialsTest_logE,2))
hold off

subplot(1,3,2)
scatter(numTrials_to_test,std(numTrialsTest_logE,0,2))
ylabel('log evidence std')
xlabel('num trials')
hold on
plot(numTrials_to_test,std(numTrialsTest_logE,0,2))
hold off

subplot(1,3,3)
scatter(numTrials_to_test,median(numTrials_time_results,2)./3600)
ylabel('hours runtime, dots: median, line: mean')
xlabel('num trials')
hold on
plot(numTrials_to_test,mean(numTrials_time_results,2)./3600)
hold off

%%

for idx=1:size(numTrialsTest_logE,2)

    scatter(numTrials_to_test,numTrialsTest_logE(:,idx))
    hold on
    
end

scatter(numTrials_to_test,median(numTrialsTest_logE,2),135,'filled')
ylabel('log evidence, dots: median, line: mean')
xlabel('num trials')
hold on
plot(numTrials_to_test,mean(numTrialsTest_logE,2))
hold off

%%

boxplot(numTrialsTest_logE',numTrials_to_test)


