classdef misc_emissions_agu_plots
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Property1
    end
    
    methods(Static)
        
        function varargout = add_weekend_series(weekday_fig, weekend_fig, varargin)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            p = advInputParser;
            p.addParameter('plot_args',{});
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            pout = p.Results;
            
            add_args = update_params('remove', varargin, pout);
            plot_args = pout.plot_args;
            
            if iscell(weekday_fig)
                loc_names = weekday_fig;
                wkday_figs = misc_emissions_analysis.plot_oh_conc_by_year('loc_inds', loc_names, 'days_of_week', 'TWRF', plot_args{:});
                wkend_figs = misc_emissions_analysis.plot_oh_conc_by_year('loc_inds', loc_names, 'days_of_week', 'US', plot_args{:});
                
                for i_fig = 1:numel(wkday_figs)
                    misc_emissions_agu_plots.add_weekday_series_main(wkday_figs(i_fig), wkend_figs(i_fig), add_args{:});
                end
                close(wkend_figs);
                varargout{1} = wkday_figs;
            else
                misc_emissions_agu_plots.add_weekday_series_main(weekday_fig, weekend_fig, add_args{:});
            end
            
        end
        
        function save_oh_fig_steps(fig, basename)
            
            % turn off legend, if present
            leg = findobj(fig,'type','legend');
            if ~isempty(leg)
                leg.Visible = 'off';
            end
            
            steps = struct('i', {{'Lifetime + SS'}},...
                'i2', {{'Lifetime + SS', 'Lifetime + SS (wkend)'}},...
                'iw', {{'Lifetime + SS', 'WRF-Chem'}},...
                'iwc', {{'Lifetime + SS', 'WRF-Chem', 'WRF (city center)'}});
            
            ch = fig.Children;
            ax = ch(end);
            names = cell(size(ax.Children));
            n_named = 0;
            for i=1:numel(names)
                names{i} = ax.Children(i).DisplayName;
                n_named = n_named + ~isempty(names{i});
            end
            
            % find the weekend VCD handle
            vcd_ax = ch(1);
            for i=1:numel(vcd_ax.Children)
                if strcmp(vcd_ax.Children(i).Marker, 'diamond')
                    weekend_vcd = vcd_ax.Children(i);
                    break;
                end
            end
            
            suffixes = fieldnames(steps);
            for i=1:numel(suffixes)
                suffix = suffixes{i};
                series = steps.(suffix);
                
                no_weekend = ~any(regcmp(series, 'wkend'));
                if no_weekend
                    weekend_vcd.Visible = 'off';
                else
                    weekend_vcd.Visible = 'on';
                end
                for j=1:numel(ax.Children)
                    if ismember(names{j}, series) || (isempty(names{j}) && ismember(names{j+n_named}, series))
                        ax.Children(j).Visible = 'on';
                    else
                        ax.Children(j).Visible = 'off';
                    end
                    
                end
                
                savename = sprintf('%s-%s.png', basename, suffix);
                saveas(fig, savename);
            end
        end
    end
    
    methods(Static, Access=protected)
        function add_weekday_series_main(weekday_fig, weekend_fig, varargin)
            p = advInputParser;
            p.addParameter('weekend_series', 'Lifetime + SS');
            p.addParameter('keep_series', {'Lifetime + SS', 'WRF-Chem', 'WRF (city center)'});
            p.parse(varargin{:});
            pout = p.Results;
            
            weekend_series = pout.weekend_series;
            keep_series = pout.keep_series;
            
            % first get rid of the second y axes on the bottom fig
            delete(weekday_fig.Children(1));
            
            % then find out how many named children of the top axes there
            % are. 
            n_named = 0;
            ax_oh = weekday_fig.Children(end);
            ax_vcd = weekday_fig.Children(1);
            names = cell(size(ax_oh.Children));
            for i_ch = 1:numel(ax_oh.Children)
                if n_named > 0 && isempty(ax_oh.Children(i_ch).DisplayName)
                    error('Unnamed child follows a named child - not expected!')
                end
                
                n_named = n_named + ~isempty(ax_oh.Children(i_ch).DisplayName);
                names{i_ch} = ax_oh.Children(i_ch).DisplayName;
            end
            
            % find the children to delete.
            xx_del = find(~ismember(names, keep_series) & ~cellfun(@isempty, names));
            % also delete the error bars which are assumed to come first
            % and be in the same order
            n_del = numel(xx_del);
            xx_del = veccat(xx_del, xx_del - n_named);
            delete(ax_oh.Children(xx_del));
            n_named = n_named - n_del;
            
            % get the weekend OH and NO2 VCDs
            ax_end = weekend_fig.Children(end);
            n_named_end = numel(find_named_children(ax_end));
            for i_ch = 1:numel(ax_end.Children)
                ch = ax_end.Children(i_ch);
                if strcmp(ch.DisplayName, weekend_series)
                    x_end = ch.XData;
                    y_end = ch.YData;
                    x_err_end = ax_end.Children(i_ch - n_named_end).XData;
                    y_err_end = ax_end.Children(i_ch - n_named_end).YData;
                    break
                end
            end
            
            no2_end = weekend_fig.Children(2).Children.YData;
            no2_color = weekend_fig.Children(2).Children.MarkerFaceColor;
            % Add the weekend OH and weekend NO2 to the weekday axes
            common_opts = {'marker', 'd', 'markersize', 10, 'markeredgecolor','k', 'linewidth', 2};
            l=line(ax_oh, x_end + 0.2, y_end, 'color', 'g', 'markerfacecolor', 'g', common_opts{:});
            l.DisplayName = sprintf('%s (wkend)', weekend_series);
            lerr=line(ax_oh, x_err_end + 0.2, y_err_end, 'color', 'g', 'linewidth', 2);
            lerr.DisplayName = '';
            line(ax_vcd, x_end + 0.2, no2_end, 'color', no2_color, 'markerfacecolor', no2_color, common_opts{:});
            
            % move the weekend series down so that the relationship between
            % series and error bars in the children is constant
            while unnamed_below(ax_oh.Children, l)
                uistack(l ,'down')
            end
            
            % set the font sizes to 16, remove "(weekdays)" from the title,
            % and get rid of the tick labels on the top plot.
            ax_oh.FontSize = 16;
            ax_oh.XTickLabel = repmat({''}, 1, length(ax_oh.XTickLabel));
            tstr = ax_oh.Title.String;
            ax_oh.Title.String = regexprep(tstr,'\s*\(\w*\)','');
            
            ax_vcd.FontSize = 16;
            
            reset_h(ax_oh)
            reset_h(ax_vcd);
            
            function reset_h(h)
                props = {'YLimMode','YTickMode','YTickLabelMode'};
                for prp=props
                    h.(prp{1}) = 'auto';
                    drawnow
                    h.(prp{1}) = 'manual';
                end
            end
            
            function yn = unnamed_below(ch, l)
                xx = find(ch == l);
                yn = isempty(ch(xx+1).DisplayName);
            end
        end
    end
    
end

