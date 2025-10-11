% To analyze the binary clusters formed as a function of time and also the averaged value for the last time block.

clc;
clear all;
close all;
format long;

%% Flags
plot_time  = 0; % Plot all clusters as a function of time
plot_equil = 1; % Plot only equilibrium data

%% Inputs
basename = 'AcONH2';
molratio = [13, 14, 15, 17];
thresh_cnt  = 100; % Maximum possible clusters
thresh_frac = 0.1; % Plot fractions only above this value

maindirname = sprintf('../../sim_results/allresults_%s',basename);
outmaindir  = sprintf('../../analyzed_results/allresults_%s',basename);
outfigdir   = sprintf('../../figures/%s',basename);

for mcnt = 1:length(molratio)

  printf('Analyzing mol-ratio: %d\n', molratio(mcnt));

  % Check for files
  mainpath = sprintf('%s/results_%d', maindirname,molratio(mcnt));
  binclust_fnames  = dir(sprintf('%s/clust_12*.dat',mainpath));

  if isempty(polyclust_fnames)
    printf('No matching files found for polytype in %s \n', mainpath);
    continue
  end

  % Write output to file
  if ~exist(outmaindir,'dir')
    mkdir(outmaindir)
  endif
  fout = fopen(sprintf('%s/polyclust_time_%d.dat',outmaindir,molratio(mcnt)),'w');
  fprintf(fout,'%s\t %s\t %s\t %s\t %s\n','ID_Tot', 'ID_1','ID_2','ID_3','Fraction')

  % Generate polyclus-time arrays
  binclust_time  = zeros(thresh_cnt,4+length(polyclust_fnames)); %ID, IDx, IDy, IDz, f1, f2 ...
  binclust_cnt   = 0; % To track how many clusters are formed

  % Analyze files one-by-one
  for fcnt = 1:length(polyclust_fnames)

    binc_fname   = sprintf('%s',polyclust_fnames(fcnt).name);
    binclust_arr = dlmread(sprintf('%s/%s',mainpath,binc_fname),'',0,0); %skip header lines

    for idcnt = 1:length(binclust_arr(:,1));

      idtot = find(binclust_arr(idcnt,1) == binclust_time(:,1)); % Check whether id is present in binclust_arr

      if numel(idtot) > 1
        error('Non-unique value of ID in the array for %d: ', binclust_arr(idcnt,1));
      end

      if isempty(idtot)
        binclust_cnt++;
        binclust_time(binclust_cnt,1) = binclust_arr(idcnt,1); % Store ID to binclust_time
        binclust_time(binclust_cnt,2:4) = binclust_arr(idcnt,2:4); % Store X,Y,Z Ids to 2:4 columns
        binclust_time(binclust_cnt,4+fcnt) = binclust_arr(idcnt,5); % Store fraction to the corresponding time column
      else
        binclust_time(idtot,4+fcnt) = binclust_arr(idcnt,5); % Store fraction to the corresponding time column
      end

    end

    if binclust_cnt > thresh_cnt
      disp('Warning: There are clusters more than %d. Consider changing thresh_cnt', thresh_cnt);
    end

  end

  % Remove rows with 0 as first column and then sort w.r.t IDs

  binclust_time_filt = binclust_time(binclust_time(:,1) > 0, :);
  binclust_sorted = sortrows(binclust_time_filt,1);

  % Write to file
  for wcnt = 1:binclust_cnt

    % Write IDs
    for cid = 1:4
      fprintf(fout, '%d\t %d\t %d\t %d\t ', binclust_sorted(wcnt,cid));
    end

    % Write values
    for fcnt = 1:length(polyclust_fnames)
      fprintf(fout,'%g\t ', binclust_time(wcnt,4+fcnt));
    end

    fprintf(fout,'\n');

  end

  fclose(fout); % Close file
end


%% Consolidate averaged data as a function of molratio
binclust_avg = zeros(length(molratio)*thresh_cnt,4+length(molratio));

% Analyze files
for mcnt = 1:length(molratio)
  fname = sprintf('%s/polyclust_time_%d.dat',outmaindir,molratio(mcnt));
  binplot_arr = dlmread(fname,'',1,0); %skip header lines
  len_line = length(binplot_arr(1,:));
  id_bycnt = 0;

  for bycnt = 1:length(binplot_arr(:,1))

    if max(binplot_arr(bycnt,5:len_line) > thresh_frac)
      idtot = find(binplot_arr(bycnt,1) == binclust_avg(:,1)); % Check whether id is present in binclust_avg
      if isempty(idtot) % new id - not found before
        id_bycnt++;
        binclust_avg(id_bycnt,1:4)    = binplot_arr(bycnt,1:4); % Store IDs
        binclust_avg(id_bycnt,4+mcnt) = binplot_arr(bycnt,len_line); % Store value
      else
        binclust_avg(idtot,4+mcnt)    = binplot_arr(bycnt,5); % Store value to the corresponding row
      end
    endif

    if id_bycnt > thresh_cnt
      disp('Warning: There are clusters more than %d for all mol-ratios combined. Consider changing thresh_cnt', thresh_cnt);
    end

  end

end

% Write to o/p file
fconsout = fopen(sprintf('%s/polyclust_avg_all.dat',outmaindir),'w');
fprintf(fout,'%s\t %s\t %s\t %s\t ','ID_Tot', 'ID_1','ID_2','ID_3')
for mcnt = 1:length(molratio) % Write headers
  fprintf(fout,'%g\t ', molratio(mcnt)/10);
end
fprintf(fout,'\n');
len_cat = find(binclust_avg(:,1) == 0)(1); % Find appearance of the first zeros
fmt = [repmat("%g\t ", 1, columns(binclust_avg)-1), "%g\n"];
for catcnt = 1:len_cat-1
    fprintf(fconsout,fmt, binclust_avg(catcnt,:));
end
fclose(fconsout);


%% Plot data - as a function of time

if ~exist(outfigdir,'dir')
  mkdir(outfigdir)
endif

clr_arr = {'g','b','r','m','#800080','k','#FFA500','#A52A2A','#006400','c','#F987C5','#40E0D0','#BDB76B','#40E0D0'};

if plot_time
  for mcnt = 1:length(molratio)

    h1 = figure; ax1 = gca; hold on; box on;
    h2 = figure; ax2 = gca; hold on; box on;

    set(ax1, 'FontSize', 16);
    set(ax2, 'FontSize', 16);
    xlabel(ax1,'Time Block','FontSize',20,'Interpreter','Tex')
    ylabel(ax1,'f_{c}','FontSize',20,'Interpreter','Tex')
    xlabel(ax2,'Time Block','FontSize',20,'Interpreter','Tex')
    ylabel(ax2,'f_{c}','FontSize',20,'Interpreter','Tex')

    leg1_arr = {}; leg2_arr = {};
    leg1cnt = 1; leg2cnt = 1;

    fname = sprintf('%s/polyclust_time_%d.dat',outmaindir,molratio(mcnt));
    binplot_arr = dlmread(fname,'',1,0); %skip header lines
    len_line = length(binplot_arr(1,:));
    tarr = 1:1:(len_line-4);

    for bycnt = 1:length(binplot_arr(:,1))

      if max(binplot_arr(bycnt,5:len_line) > thresh_frac)
        plot(ax1,tarr,binplot_arr(bycnt,5:len_line),'color',clr_arr{leg1cnt},'Marker','o','MarkerSize',10,'LineStyle','--')
        leg1_arr{leg1cnt} = ['Al' num2str(binplot_arr(bycnt,2)) 'Cl' num2str(binplot_arr(bycnt,3)) baseshort num2str(binplot_arr(bycnt,4))];
        leg1cnt++;
      else
        plot(ax2,tarr,binplot_arr(bycnt,5:len_line),'color',clr_arr{leg2cnt},'Marker','d','MarkerSize',10,'LineStyle','--')
        leg2_arr{leg2cnt} = ['Al' num2str(binplot_arr(bycnt,2)) 'Cl' num2str(binplot_arr(bycnt,3)) baseshort num2str(binplot_arr(bycnt,4))];
        leg2cnt++;
      end

    end

    if leg1cnt > 1
      legend(ax1,leg1_arr,'location','bestoutside','FontSize',16)
      saveas(h1,sprintf('%s/polyclust_%d_gt_thresh', outfigdir,molratio(mcnt)),'png')
    end

    if leg2cnt > 1
      legend(ax2,leg2_arr,'location','bestoutside','FontSize',16)
      saveas(h2,sprintf('%s/polyclust_%d_leq_thresh', outfigdir,molratio(mcnt)),'png')
    end

    close(h1); close(h2);

  end
end

%% Plot bar chart - for the last case (Plot only values > threshold)
% Group different clusters for all molratios

if plot_equil
  % Plot binclust_avg as a grouped bar chart
  h1 = figure; hold on; box on;
  set(gca, 'FontSize', 16);
  len_cat = find(binclust_avg(:,1) == 0)(1); % Find appearance of the first zeros
  categories = {};
  for catcnt = 1:len_cat-1
    categories{catcnt} = ['Al' num2str(binclust_avg(catcnt,2)) 'Cl' num2str(binclust_avg(catcnt,3)) baseshort num2str(binclust_avg(catcnt,4))];
  end
  leg_arr = arrayfun(@num2str, molratio/10, 'UniformOutput', false);
  bar(binclust_avg(1:len_cat-1,5:4+mcnt)); % Plot all data
  set(gca, "XTick", 1:len_cat-1);
  set(gca, "XTickLabel", categories);
  xt  = get(gca,"XTick");
  xt1 = get(gca,"XTickLabel");
  set(gca, "XTickLabel", []);  % Remove default labels
  % Add vertical labels
%  for i = 1:length(xt)
%    text(xt(i), -0.1, categories{i}, 'Rotation', 90, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle','FontSize',16);
%  end
  xlabel('Cluster Name','FontSize',20,'Interpreter','Tex');
  ylabel('n(s)','FontSize',20,'Interpreter','Tex');
  legend(leg_arr);
  xtickangle(90);  % Rotate x-tick labels by 90 degrees
  saveas(h1,sprintf('%s/polyclust_all_equil_gt_thresh', outfigdir),'png')
end

