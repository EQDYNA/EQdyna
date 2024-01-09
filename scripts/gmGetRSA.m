function gmGetRSA(dx, dt, gamma, np, T, range)
    % gmGetRSA.m converts binary ground motion output gm* from EQdyna to MATLAB
    % .mat format. Dunyu Liu <dliu@ig.utexas.edu>, 20210715.
    % It calls ResSpecRotD50.m and ResSpecTimeDomVectorized_3.m written by Dr. Steven Day.

    NintMax = 5;
    g = 9.8;
    x0 = range(1); %-40; % km
    x1 = range(2); %40; % km
    y0 = range(3); %-40; % km
    y1 = range(4); %41; % km
    dx = dx/1000; % km
    
    for me = 0: np-1
        ntag  = 0;
        fname = strcat('surface_coor.txt',num2str(me));
        if exist(fname, 'file')
            printTag = strcat('Generating rsa.mat from MPI process ', num2str(me))
            loaddata = 1;
            a        = load(fname);
            [n, m]   = size(a);
            xtmp     = a(1:n,1);
            ytmp     = a(1:n,2);
            for i = 1:n
                if abs(xtmp(i)/1e3) < x1 && abs(ytmp(i)/1e3) < y1
                    loaddata = 1;
                end
            end
            if loaddata == 1
                xcoor(1:n) = a(1:n,1);
                ycoor(1:n) = a(1:n,2);
                fname1     = strcat('gm',num2str(me));
                fileID     = fopen(fname1);
                C          = fread(fileID, 'double');
                n1         = size(C,1);
                nt         = n1/n/3;
                
                for i = 1:n % loop over stations
                    %if abs(xcoor(initial_tag+i-1)/1e3) < 1 && (abs(ycoor(initial_tag+i-1)/1e3) < 1 ||abs(ycoor(initial_tag+i-1)/1e3-40) < 1)
                     if abs(xcoor(i)/1e3) < x1 && abs(ycoor(i)/1e3) < y1
                        ntag = ntag + 1;
                        coor(ntag,1) = xcoor(i);
                        coor(ntag,2) = ycoor(i);
                        for j = 1: nt % loop over time steps
                            timeseries_x(i,j) = C((j-1)*n*3 + (i-1)*3 + 1);
                            timeseries_y(i,j) = C((j-1)*n*3 + (i-1)*3 + 2); 
                        end
                        acc_x = vel_to_acc(timeseries_x(i,:)',dt);
                        acc_y = vel_to_acc(timeseries_y(i,:)',dt);
                        
                        sa_rotd50 = ResSpecRotD50(acc_x', acc_y', dt ,T ,gamma, NintMax);
                        res_sa(ntag,:) = sa_rotd50; 
                    end
                end
                res_sa = res_sa/g; % in g
                save(strcat('rsa',num2str(me),'.mat'),'res_sa','coor');
                clear C xcoor ycoor;
            end
            clear a;
        end
    end
end