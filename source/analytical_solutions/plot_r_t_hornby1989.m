% plot figure 2 at Hornby 1989.
f    = linspace(0, 1000, 1000);
a   = 0.1;
D   = [0.1:0.1:1]/100;
c    = 1500;
figure
for i =1 : length(D)
    d = D(i);
    [r, t] = reflection_transmission(f, a, d, c);
    plot(f/1e3, abs(r),'k-'); hold on;
    xlim([0, 1.0]);xlabel('frequency (Hz)');
    ylim([0, 1.0]); ylabel('|r|');
end
hold off;
figure
for i =1 : length(D)
    d = D(i);
    [r, t] = reflection_transmission(f, a, d, c);
    plot(f/1e3, abs(t),'k-'); hold on;
    xlim([0, 1.0]);xlabel('frequency (Hz)');
    ylim([0, 1.0]); ylabel('|t|');
end
hold off;

