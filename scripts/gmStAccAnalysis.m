function gmStAccAnalysis(dt, acc, figID, stName)


    % 
    % acc = load("st.x0.y10.accx.txt");
    n = length(acc);
    t = dt:dt:n*dt;
    
    y  = fft(acc); 
    fs = 1/dt; % sample frequency
    f = (0:n-1)*(fs/n); % frequency range
    power = abs(y).^2/n; % power of the DFT
    
    figure(figID)
    subplot(2,1,1)
    plot(t, acc, 'k','LineWidth',1.2);
    xlabel('Time (s)');
    ylabel('Acceleration (m/s/s)');
    title(stName);
    set(gca, 'FontSize', 10, 'FontWeight', 'bold');

    subplot(2,1,2)
    plot(f,power, 'k','LineWidth',1.2);
    xlabel('Frequency');
    ylabel('Power');
    set(gca, 'FontSize', 10, 'FontWeight', 'bold');
    set(gcf,"Color","white");
end