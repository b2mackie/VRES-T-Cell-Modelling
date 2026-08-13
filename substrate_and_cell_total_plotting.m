%% Getting Data from Output

timetotal =336; % set the total time steps of the data

%automating the loading of the outputs
A1 = 'output0000000';
A2 = 'output000000';
A3 = 'output00000'; 
A4 = 'output0000'; 
B = '.xml';

for tcount = 2:timetotal

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
    MCDS = read_MultiCellDS_xml(K,'C:\Users\b14ma\OneDrive\Documents\University\2026\MY RES\outputs\output_14days');
    k = find( MCDS.mesh.Z_coordinates == 0 ); 

    cell_type_0(tcount-1) = length(find(MCDS.discrete_cells.metadata.type==0));
    cell_type_1(tcount-1) = length(find(MCDS.discrete_cells.metadata.type==1));
 	 

    %active v inactive v exhausted cells 

    types = MCDS.discrete_cells.metadata.type;
    T_c = (types == 0);
    state = MCDS.discrete_cells.custom.state;
    t_states = state(T_c);
    
    ExhaustedT_all(tcount)   = sum(t_states > 1.5);
    InactiveT_all(tcount) = sum(t_states < 0.5);
    ActiveT_all(tcount) = sum(t_states == 1); 

    % Contact time 

    contact_time = MCDS.discrete_cells.custom.attached_time;
    t_contact_time = contact_time(T_c);

end

%% Plotting

% total cells 
figure
plot(cell_type_0,'LineWidth',2)
xlabel('Time')
ylabel('Cell type 0')
set(gca,'FontSize',18)

% total rod 
figure
plot(cell_type_1,'LineWidth',2)
xlabel('Time')
ylabel('Cell type 1')
set(gca,'FontSize',18)

% Active v inactive v exhuasted cells 
figure 
plot(ActiveT_all, 'LineWidth',2)
hold on
plot(InactiveT_all, '--', 'LineWidth',2)
hold on 
plot(ExhaustedT_all, '.', 'Linewidth',2)
xticks([1 25 49 73 97 121 145 169 193 217 241 265 289 313 336])
xticklabels({'0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14'})
title('Active v Inactive v Exhausted T Cells', 'a_{max} = 1, s = 1')
xlabel('Time (days)')
ylabel('T Cell Count')
legend('Active', 'Inactive', 'Exhausted', 'Location','best')
set(gca,'FontSize',18)

% Contact time histogram 
figure, 
histogram(t_contact_time, 20)
hold on 
xline(720)
xlabel('contact time')
ylabel('number of t cells')
title('total contact time of t cells with microrods')

