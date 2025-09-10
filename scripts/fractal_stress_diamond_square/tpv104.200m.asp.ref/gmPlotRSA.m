function gmPlotRSA(dx, np, T, range, xLimit, RBinSize, Rmin, Rmax)
    % gmPlotRSA.m plots RSA distribution on the surface.
    % Dunyu Liu <dliu@ig.utexas.edu>, 20210715.
    % It loads rsa*.mat, where * is MPI process ID. rsa*.mat is processed by
    % gmGetRSA.m, which calls ResSpecRotD50.m and ResSpecTimeDomVectorized_3.m by Dr. Steven Day.
    % It calls sa_Rrupt.m to get SA vs. rupture distance.

    %steps = 10;
    %dt = dt*steps;
    NintMax = 5;
    g = 9.8;
    xs0 = range(1); % km
    xs1 = range(2); % km
    ys0 = range(3); % km
    ys1 = range(4); % km
    dx  = dx/1000; % km
    nT  = length(T);

    itag = 0;
    for me = 0: np-1
        fname=strcat('rsa',num2str(me),'.mat');
        if exist(fname, 'file')
            printTag = strcat('Processing rsa.mat from MPI process ', num2str(me))
            itag = itag + 1;
            a = load(fname);
            tmp1 = a.coor;
            tmp2 = a.res_sa;
            if itag == 1
                res_sa = tmp2;
                coor = tmp1;
            else
                res_sa = [res_sa; tmp2;];
                coor = [coor; tmp1;];
            end
            clear a;
        end
    end
    
    % remove redudant data points (rows)
    [coor1, id] =  unique(coor, 'rows');
    for i = 1: nT
        tmp = res_sa(:,i);
        tmp1 = tmp(id);
        SA(:,i) = tmp1;
    end
    
    [xx,yy] = meshgrid(xs0:dx:xs1,ys0:dx:ys1);
    
    h1 = figure(1);
    set(h1, 'Position', [10 10 500 800]);
    for i = 1: nT
        F = scatteredInterpolant(coor1(:,1)/1e3,coor1(:,2)/1e3, SA(:,i),'linear');
        saMap = F(xx,yy);
        subplot(nT,1,i)
        contourf(xx, yy, saMap, 'LineStyle','none'); colorbar;
        axis equal;
        title(strcat('SA (g) at ', num2str(T(i)),' s'));
        if i == floor(nT/2)
            ylabel('Fault-normal (km)');
        elseif i == 3
            xlabel('Strike (km)');
        end
        set(gca, 'Fontsize', 12, 'Fontweight', 'bold');
    end
    set(gcf, 'color', 'white');
    
    %%
    [dmean, dmin, dmax, Rrup, period_list] = empirical_kevin;
    for i = 1: nT
        for j = 1: 21
            if T(i) == period_list(j)
                k(i) = j; % find the index in Kevin empirical period list
            end
        end
    end 

    h2=figure(2);
    set(h2, 'Position', [10 10 600 200*nT]);
    for i = 1: nT
        [stat, R, numSA, r, saRjb] = gmGetSARrupt(coor1, SA(:,i), xLimit, RBinSize, Rmin, Rmax);
    
        subplot(nT,1,i)
        plot(Rrup,dmean(:,k(i)),'k-', 'LineWidth', 1); hold on;
        plot(Rrup,dmax(:,k(i)),'k:', 'LineWidth', 1); hold on;
        plot(Rrup,dmin(:,k(i)),'k:', 'LineWidth', 1); hold on;
        plot(R,stat(:,1),'o-'); hold on;
        plot(R,stat(:,3),'*'); hold on;
        plot(R,stat(:,4),'^'); hold on;
        plot(R,stat(:,5),'^'); hold on;
        plot(R,stat(:,6),'*'); hold on;
        % for iR = 1: length(R)
        %     for iS = 1: numSA(iR)
        %         plot(R(iR), saRjb(iR,iS),'.'); hold on;
        %     end
        % end
        ylim([0.01 1]);
        title(strcat('SA GMRotD50 (g) at ', num2str(T(i)),' s'));
        set(gca, 'YScale', 'log', 'XScale', 'log', 'color', 'white');
        set(gcf, 'color', 'white');
        set(gca, 'Fontsize', 12, 'Fontweight', 'bold');
    end

    h3=figure(3);
    set(h3, 'Position', [10 10 600 200*nT]);
    for i = 1: nT
        subplot(nT,1,i)
        plot(Rrup,0.6,'k-', 'LineWidth', 1); hold on;
        plot(R,stat(:,2),'o-','Linewidth',2);
        ylim([0 1]);
        title(strcat('Intra-event variability of log(SA) at ', num2str(T(i)),' s'));
        set(gcf, 'color', 'white');
        set(gca, 'Fontsize', 12, 'Fontweight', 'bold');
    end
end