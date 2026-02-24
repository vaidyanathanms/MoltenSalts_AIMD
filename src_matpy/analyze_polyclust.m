% To analyze the polyclusters (AlClx.Base) formed as a function of time and also the averaged value for the last time block.

clc;
clear all;
close all;
format long;

%% Flags
plot_time  = 1; % Plot all clusters as a function of time
plot_equil = 1; % Plot only equilibrium data

%% Inputs
basename = 'Urea'; baseshort = '.Ur';
molratio = [13, 15, 17];
thresh_cnt  = 100; % Maximum possible clusters
thresh_frac = 0.1; % Plot fractions only above this value
caseid      = 3; % case number

maindirname = sprintf('../../sim_results/allresults_%s',basename);
outmaindir  = sprintf('../../analyzed_results/allresults_%s',basename);
outfigdir   = sprintf('../../figures/%s',basename);

for mcnt = 1:length(molratio)

  printf('Analyzing mol-ratio in Case-: %d \t %d\n', molratio(mcnt), caseid);

  % Check for files
  mainpath = sprintf('%s/results_%d/Case_%d', maindirname,molratio(mcnt),caseid);
  polyclust_fnames = dir(sprintf('%s/polyclust_*',mainpath));

  if isempty(polyclust_fnames)
    printf('No matching files found for polytype in %s \n', mainpath);
    continue
  end

  % Write output to file
  if ~exist(outmaindir,'dir')
    mkdir(outmaindir)
  endif
  fout = fopen(sprintf('%s/polyclust_time_%d_case_%d.dat',outmaindir,molratio(mcnt),caseid),'w');
  fprintf(fout,'%s\t %s\t %s\t %s\t %s\n','ID_Tot', 'ID_1','ID_2','ID_3','Fraction')

  % Generate polyclus-time arrays
  pclust_time  = zeros(thresh_cnt,4+length(polyclust_fnames)); %ID, IDx, IDy, IDz, f1, f2 ...
  pclust_cnt   = 0; % To track how many clusters are formed

  % Analyze files one-by-one
  for fcnt = 1:length(polyclust_fnames)

    polyc_fname   = sprintf('%s',polyclust_fnames(fcnt).name);
    polyclust_arr = dlmread(sprintf('%s/%s',mainpath,polyc_fname),'',0,0); %skip header lines

    for idcnt = 1:length(polyclust_arr(:,1));

      idtot = find(polyclust_arr(idcnt,1) == pclust_time(:,1)); % Check whether id is present in polyclust_arr

      if numel(idtot) > 1
        error('Non-unique value of ID in the array for %d: ', polyclust_arr(idcnt,1));
      end

      if isempty(idtot)
        pclust_cnt++;
        pclust_time(pclust_cnt,1) = polyclust_arr(idcnt,1); % Store ID to pclust_time
        pclust_time(pclust_cnt,2:4) = polyclust_arr(idcnt,2:4); % Store X,Y,Z Ids to 2:4 columns
        pclust_time(pclust_cnt,4+fcnt) = polyclust_arr(idcnt,5); % Store fraction to the corresponding time column
      else
        pclust_time(idtot,4+fcnt) = polyclust_arr(idcnt,5); % Store fraction to the corresponding time column
      end

    end

    if pclust_cnt > thresh_cnt
      disp('Warning: There are clusters more than %d. Consider changing thresh_cnt', thresh_cnt);
    end

  end

  % Remove rows with 0 as first column and then sort w.r.t IDs

  pclust_time_filt = pclust_time(pclust_time(:,1) > 0, :);
  pclust_sorted = sortrows(pclust_time_filt,1);

  % Write to file
  for wcnt = 1:pclust_cnt

    % Write IDs
    for cid = 1:4
      fprintf(fout, '%d\t %d\t %d\t %d\t ', pclust_sorted(wcnt,cid));
    end

    % Write values
    for fcnt = 1:length(polyclust_fnames)
      fprintf(fout,'%g\t ', pclust_time(wcnt,4+fcnt));
    end

    fprintf(fout,'\n');

  end

  fclose(fout); % Close file
end


%% Consolidate averaged data as a function of molratio
pclust_avg = zeros(length(molratio)*thresh_cnt,4+length(molratio));
id_plcnt = 0;

% Analyze files
for mcnt = 1:length(molratio)
  fname = sprintf('%s/polyclust_time_%d_case_%d.dat',outmaindir,molratio(mcnt),caseid);
  polyplot_arr = dlmread(fname,'',1,0); %skip header lines
  len_line = length(polyplot_arr(1,:));

  for plcnt = 1:length(polyplot_arr(:,1))

    if max(polyplot_arr(plcnt,5:len_line) > thresh_frac)
      % Number of aluminums are different. Only way to check is to check all 3 IDs.
      ref_col1 = polyplot_arr(plcnt,2:4); % column in the file
      avg_col2 = pclust_avg(1:id_plcnt+1,2:4); % averaged array
      match_idx = find(ismember(avg_col2,ref_col1,"rows"));
      if numel(match_idx) > 1
        error('Non-unique value of ID in the match_idx');
      end

      if isempty(match_idx) % new id - not found before
        id_plcnt++;
        pclust_avg(id_plcnt,1)      = id_plcnt; % Assign new ID
        pclust_avg(id_plcnt,2:4)    = polyplot_arr(plcnt,2:4); % Store IDs
        pclust_avg(id_plcnt,4+mcnt) = polyplot_arr(plcnt,len_line); % Store last value
      else
        pclust_avg(match_idx,4+mcnt)    = polyplot_arr(plcnt,len_line); % Store value to the corresponding row
      end
    endif

    if id_plcnt > thresh_cnt
      disp('Warning: There are clusters more than %d for all mol-ratios combined. Consider changing thresh_cnt', thresh_cnt);
    end

  end

end

% Remove rows with 0 as first column and then sort w.r.t Aluminum IDs
pclust_avg = pclust_avg(pclust_avg(:,1) > 0, :);
pclust_avg_sorted = sortrows(pclust_avg,2); % Sort w.r.t Al

% Write to o/p file
fconsout = fopen(sprintf('%s/polyclust_avg_all_case_%d.dat',outmaindir,caseid),'w');
fprintf(fout,'%s\t %s\t %s\t %s\t ','ID_Tot', 'ID_1','ID_2','ID_3')
for mcnt = 1:length(molratio) % Write headers
  fprintf(fout,'%g\t ', molratio(mcnt)/10);
end
fprintf(fout,'\n');
fmt = [repmat("%g\t ", 1, columns(pclust_avg_sorted)-1), "%g\n"];
for catcnt = 1:length(pclust_avg_sorted)
    fprintf(fconsout,fmt, pclust_avg_sorted(catcnt,:));
end
fclose(fconsout);


%% Plot data - as a function of time

if ~exist(outfigdir,'dir')
  mkdir(outfigdir)
endif

clr_arr = {'g','b','r','m','#800080','k','#FFA500','#A52A2A','#006400','c','#F987C5','#40E0D0','#BDB76B','#40E0D0', ...
'g','b','r','m','#800080','k','#FFA500','#A52A2A','#006400','c','#F987C5'};

if plot_time
  for mcnt = 1:length(molratio)

    h1 = figure; ax1 = gca; hold on; box on;
    h2 = figure; ax2 = gca; hold on; box on;

    set(ax1, 'FontSize', 16);
    set(ax2, 'FontSize', 16);
    xlabel(ax1,'Time Block','FontSize',20,'Interpreter','Tex')
    ylabel(ax1,'N(s)','FontSize',20,'Interpreter','Tex')
    xlabel(ax2,'Time Block','FontSize',20,'Interpreter','Tex')
    ylabel(ax2,'N(s)','FontSize',20,'Interpreter','Tex')

    leg1_arr = {}; leg2_arr = {};
    leg1cnt = 1; leg2cnt = 1;

    fname = sprintf('%s/polyclust_time_%d_case_%d.dat',outmaindir,molratio(mcnt),caseid);
    polyplot_arr = dlmread(fname,'',1,0); %skip header lines
    len_line = length(polyplot_arr(1,:));
    tarr = 1:1:(len_line-4);

    for plcnt = 1:length(polyplot_arr(:,1))

      if max(polyplot_arr(plcnt,5:len_line) > thresh_frac)
        plot(ax1,tarr,polyplot_arr(plcnt,5:len_line),'color',clr_arr{leg1cnt},'Marker','o','MarkerSize',10,'LineStyle','--')
        leg1_arr{leg1cnt} = ['Al_' num2str(polyplot_arr(plcnt,2)) 'Cl_' num2str(polyplot_arr(plcnt,3)) baseshort '_' num2str(polyplot_arr(plcnt,4))];
        leg1cnt++;
      else
        plot(ax2,tarr,polyplot_arr(plcnt,5:len_line),'color',clr_arr{leg2cnt},'Marker','d','MarkerSize',10,'LineStyle','--')
        leg2_arr{leg2cnt} = ['Al' num2str(polyplot_arr(plcnt,2)) 'Cl' num2str(polyplot_arr(plcnt,3)) baseshort num2str(polyplot_arr(plcnt,4))];
        leg2cnt++;
      end

    end

    if leg1cnt > 1
      legend(ax1,leg1_arr,'location','bestoutside','FontSize',16,'Interpreter','Tex')
      saveas(h1,sprintf('%s/polyclust_%d_gt_thresh_case_%d', outfigdir,molratio(mcnt),caseid),'png')
    end

    if leg2cnt > 1
      legend(ax2,leg2_arr,'location','bestoutside','FontSize',16)
      saveas(h2,sprintf('%s/polyclust_%d_leq_thresh_case_%d', outfigdir,molratio(mcnt),caseid),'png')
    end

    close(h1); close(h2);

  end
end

%% Plot bar chart - for the last case (Plot only values > threshold)
% Group different clusters for all molratios

if plot_equil
  % Plot pclust_avg as a grouped bar chart
  h1 = figure; hold on; box on;
  set(gcf, 'Position', [100, 100, 1600, 1200]);  % [x, y, width, height]
  set(gca, 'FontSize', 12);
  categories = {};
  for catcnt = 1:length(pclust_avg(:,1))
    categories{catcnt} = ['Al_{' num2str(pclust_avg(catcnt,2)) '}' 'Cl_{' num2str(pclust_avg(catcnt,3)) '}' baseshort '_{' num2str(pclust_avg(catcnt,4)) '}'];
  end
  leg_arr = arrayfun(@num2str, molratio/10, 'UniformOutput', false);
  bar(pclust_avg(:,5:4+mcnt)); % Plot all data
  set(gca, "XTick", 1:length(pclust_avg));
  set(gca, "XTickLabel", categories);
  xt  = get(gca,"XTick");
  xt1 = get(gca,"XTickLabel");
  set(gca, "XTickLabel", []);  % Remove default labels
  % Add vertical labels
  for i = 1:length(xt)
    text(xt(i), -0.1, categories{i}, 'Rotation', 45, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle','FontSize',16);
  end
  %xlabel('Cluster Name','FontSize',20,'Interpreter','Tex');
  ylabel('n(s)','FontSize',20,'Interpreter','Tex');
  legend(leg_arr);
  xtickangle(90);  % Rotate x-tick labels by 90 degrees
  saveas(h1,sprintf('%s/polyclust_all_equil_gt_thresh_case_%d', outfigdir,caseid),'png')
end


