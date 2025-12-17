%% Creating movie of substrates

timetotal =120; % set the total time steps of the data

%automating the loading of the outputs
A1 = 'output0000000';
A2 = 'output000000';
A3 = 'output00000'; 
A4 = 'output0000'; 
B = '.xml';

index = 1; % set the substrate index of interest
            % 1 = IL-2
            % 2 = chemokine 

% opening a video
v = VideoWriter('video1.avi')
v.FrameRate = 10;
open(v)
figure('Position', [50 50 900 700])
hold on 

for tcount = 2:timetotal

    clf % clear figure

    % set up output name
    if tcount<11
        K = [A1 num2str(tcount-1,'%d') B];
    elseif tcount<101
        K = [A2 num2str(tcount-1,'%d') B];
    elseif tcount<1001
        K = [A3 num2str(tcount-1,'%d') B];
    else
        K = [A4 num2str(tcount-1,'%d') B];
    end

    % reading output
    MCDS = read_MultiCellDS_xml(K,'C:/Users/b14ma/Documents/University/PhysiCell/PhysiCell/output'); 
    k = find( MCDS.mesh.Z_coordinates == 0 ); 
    
    %plotting contour plot of substrate
    contourf( MCDS.mesh.X(:,:,k), MCDS.mesh.Y(:,:,k), ...
    MCDS.continuum_variables(index).data(:,:,k) , 20 ) ;
    axis image;
    colorbar; 
    xlabel( sprintf( 'x (%s)' , MCDS.metadata.spatial_units) ); 
    ylabel( sprintf( 'y (%s)' , MCDS.metadata.spatial_units) ); 
     title( sprintf('%s (%s) at t = %3.2f %s', MCDS.continuum_variables(index).name , ...
     MCDS.continuum_variables(index).units , ...
     MCDS.metadata.current_time , ...
     MCDS.metadata.time_units ) ); 
 
 	 VUC(tcount) = sum(sum(MCDS.continuum_variables(index).data(:,:,k)))*(MCDS.mesh.X_coordinates(2)-MCDS.mesh.X_coordinates(1))*(MCDS.mesh.Y_coordinates(2)-MCDS.mesh.Y_coordinates(1));


    cell_type_0(tcount-1) = length(find(MCDS.discrete_cells.metadata.type==0));
    cell_type_1(tcount-1) = length(find(MCDS.discrete_cells.metadata.type==1));
 	 
    frame = getframe(gcf);
    writeVideo(v,frame);
    

end
close(v)

% IL-2
figure
plot(VUC,'LineWidth',2)
xlabel('Time')
if index == 1
    ylabel('IL-2')
    title('IL-2 Concentration over Time')
elseif index == 2
    title('Chemokine Concentration')
    ylabel('chemokine')
end 
set(gca,'FontSize',18)

% T Cell
figure
plot(cell_type_0,'LineWidth',2)
title('Number of T Cells')
xlabel('Time')
ylabel('Cell type 0')
set(gca,'FontSize',18)

% Rod Cell 
figure
plot(cell_type_1,'LineWidth',2)
title('Number of Rod Cells')
xlabel('Time')
ylabel('Cell type 1')
set(gca,'FontSize',18)

%% Fitting to exponential

P0 = 40;
r = linspace(0, 0.0006, 1000);
t = linspace(0, 1880, length(cell_type_0));
residual = zeros(size(r));

for i = 1:length(r)
    for j = 1:length(t)
        Pbar = P0*exp(r(i)*t(j));
        P = cell_type_0(j);
        diff_sqr = (Pbar-P).^2;
        residual(i) = residual(i) + diff_sqr;
    end
end

r_good = r(find(residual == min(residual)));

figure,
plot(r, residual)
hold on
plot(r_good, min(residual), 'r*')
hold on
plot(0.000478, residual(797), 'b*')
title('Fitting growth parameter to t cell growth')
xlabel('Growth Rate')
ylabel('Sum or Residuals')
legend('Residual', 'Growth Rate (Minimum Residual)', 'r_p^{max}')

figure,
P_plot = P0*exp(r_good*t);
plot(t, P_plot)
hold on
plot(t, cell_type_0)
title('Simulated T Cell Population Against Expected Cell Growth')
xlabel('Time (minutes)')
ylabel('Cell Population')
legend('Expected Cell Growth', 'Simulated Cell Growth', 'Location', 'northwest')
