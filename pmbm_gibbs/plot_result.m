function gen_figures = plot_result( meas, truth, est, ospa_vals )

% OSPA illustration
averOspa = mean(ospa_vals,3);
sumOspa = mean(squeeze(sum(ospa_vals,1)),2)

figure
ospa = gcf;
hold on
subplot(3,1,1); plot(1:meas.K,averOspa(:,1)); grid on; ylabel('OSPA Dist');
subplot(3,1,2); plot(1:meas.K,averOspa(:,2)); grid on; ylabel('OSPA Loc');
subplot(3,1,3); plot(1:meas.K,averOspa(:,3)); grid on; ylabel('OSPA Card');
xlabel('Time');

% Cardinality illustration
figure
card = gcf;
hold on
grid on
plot(truth.N)
plot(mean(est.N,2),'*')
title('Number of target')
xlabel('Time step')
ylabel('N')
legend('real number','estimation','Location','northwest')

gen_figures = [ospa card];

end

