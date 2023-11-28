clear all; close all;

[dmean, dmin, dmax, Rrup, period] = empirical_kevin;
k = 5;

figure(1)

plot(Rrup,dmean(:,k)); hold on;
plot(Rrup,dmax(:,k),':'); hold on;
plot(Rrup,dmin(:,k),':'); hold on;
set(gca, 'YScale', 'log', 'XScale', 'log', 'color', 'white');


