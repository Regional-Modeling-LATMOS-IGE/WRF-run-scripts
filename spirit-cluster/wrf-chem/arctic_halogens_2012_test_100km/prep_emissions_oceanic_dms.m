% function [] = prep_emissions_oceanic_dms()

%--------  Adds oceanic DMS concentrations to WRF-Chem emission files (wrfchemi_d0*) -------%
  %
  % Louis Marelle, 2018/10/04
  %
  % To run from the command line:
  %   $ matlab -nodisplay -singleCompThread -nodesktop -r "prep_wrfchemi_oceanic_dms();exit"
  % I use recent matlab functions (matlab 2016), so make sure a recent version of matlab is loaded, e.g. before use:
  %   $ module unload matlab, module load matlab/2016b
  %
  % Inputs: 
  %  WRFCHEMI_PATH: path to a directory containing wrfchemi_d0* files, to which
  %    DMS oceanic emissions should be added
  %
  % Outputs:
  %   Modifies the original wrfchemi files in WRFCHEMI_PATH/
  %
  % Purpose:
  %   This MATLAB routine (function) creates daily wrfchemi_dXX input emission
  %   files for WRFChem containing E_DMS_OC, which is the oceanic DMS content from
  %   Lana et al. 2011 interpolated on the WRF grid. E_DMS_OC is used to compute
  %   oceanic dms emissions in module_nightingale_dmsemiss.F, a custom online DMS
  %   emission routine that I developped for the CBMZ-MOSAIC and SAPRC-MOSAIC
  %   mechanisms including some DMS chemistry. In the future, E_DMS_OC should
  %   probably be read in a different input file.
  %
  % This routine can be modified easily to create hourly files, by replacing the 
  % "for ihour = 0:0" loop by "for ihour = 0:23"
  %---------------------------------------------------------------------------------------------%

  WRFCHEMI_PATH = './'

  % --------   Parameters --------
  % Emission datasets directories
  LANA_DMS_DIRECTORY = '/data/marelle/EMISSIONS/DMS_LANA/';

  % Species
  % Mechanism species, written in wrfchemi file
  SPECNAMES_WRFCHEMI = {'DMS_OC'};

  % Number of vertical levels in the wrfchemi file
  zdim = 1;


  % -------- Initialization -------- 
  disp(' ')
  disp('-------- prep_wrfchemi_oceanic_dms - produce wrfchemi files containing E_DMS_OC, the DMS oceanic content --------')

  disp(' ')
  disp('Initialization')

  % Include my MATLAB toolbox
  addpath('/home/marelle/scripts/matlab/TOOLBOX');

  % Get domain info from the wrfchemi filenames in WRFCHEMI_PATH/
  disp('  Opening wrfchemi files to get WRF run dates')
  path_files = dir([WRFCHEMI_PATH, '/wrfchemi_d0*00']);
  wrfchemi_filenames = { path_files.name };
  wrfchemi_filenames = cell2mat(wrfchemi_filenames');
  % Get the max number of domains from wrfchemi filenames
  wrfchemi_domains = wrfchemi_filenames(:, 11:12);
  ndomains = max(str2num(wrfchemi_domains));

  % -------- WRF domain Loop -------- 
  for idomain = 1:ndomains
    % Domain name    
    domain = ['d0' num2str(idomain)];
    disp(' ')
    disp(['Initializing for WRF domain ', domain])

    % Get wrfchemi dates from the wrfchemi_<domain>* filenames in WRFCHEMI_PATH/
    path_files = dir([WRFCHEMI_PATH, '/wrfchemi_', domain, '*00']);
    wrfchemi_filenames = { path_files.name };
    wrfchemi_filenames = cell2mat(wrfchemi_filenames');
    wrfchemi_dates = wrfchemi_filenames(:, 14:32);
    wrfchemi_dates = datenum(wrfchemi_dates, 'yyyy-mm-dd_HH:MM:SS');
    years_list = year(wrfchemi_dates)';
    months_list = month(wrfchemi_dates)';
    days_list = day(wrfchemi_dates)';
    % List of individual simulation years
    years_list_y = years_list([true, diff(years_list)~=0]);
    % List of individual simulation months (e.g., for a run from 2012/12/1 to
    % 2013/3/1, months_list_m = [12 1 2 3])
    months_list_m = months_list([true, diff(months_list)~=0]);
    % List of years corresponding to individual simulation months (e.g., same
    % example than above, years_list_m = [2012 2013 2013 2013]
    years_list_m = years_list([true, diff(months_list)~=0]);
    % List of individual simulation days
    days_list_d = days_list([true, diff(days_list)~=0]);
    % List of months corresponding to individual simulation days
    months_list_d = months_list([true, diff(days_list)~=0]);
    % List of years corresponding to individual simulation days
    years_list_d = years_list([true, diff(days_list)~=0]);


    % -------- Initialization for domain idomain -------- 
    % Read projection & grid info from wrfinput file
    disp('  Reading WRF domain info')
    [moad_cen_lat, truelat1, truelat2, stdlon, imax, jmax, kmax, dx, dy, ref_lat, ref_lon, map_proj, hemi, ref_x, ref_y] = get_WRF_grid(['./wrfinput_', domain]);
    wrf_proj = get_WRF_proj(['./wrfinput_', domain]);

    % Get DMS climatology grid dimensions
    disp('  Reading lat & lon for the DMS climatology')
    emissions_lat_DMS = -89.5:1:89.5;
    emissions_lon_DMS = -179.5:1:179.5;

    % Read WRF domain lat, lon and map projection factors from a wrfinput file
    % Find min & max lat & lon from the WRF domain
    disp('  Opening wrfinput file to get WRF lat & lon boundaries')
    ncid = netcdf.open(['./wrfinput_', domain], 'NC_NOWRITE');
    wrf_xlat = ncread(ncid, 'XLAT');
    % Adding -1 & + 1 degrees to avoid missing emissions at the domain
    % boundaries
    lat_min = min(min(wrf_xlat)) - 1;
    lat_max = max(max(wrf_xlat)) + 1;
    wrf_xlong = ncread(ncid, 'XLONG');
    lon_min = min(min(wrf_xlong)) - 1;
    lon_max = max(max(wrf_xlong)) + 1;
    clear data
    % Read the map factors from the wrfinput file
    data_mapfacx_wrf = ncread(ncid, 'MAPFAC_MX');
    data_mapfacy_wrf = ncread(ncid, 'MAPFAC_MY');
    data_landmask_wrf = ncread(ncid, 'LANDMASK');
    netcdf.close(ncid);


    % -------- Loop over months and interpolate monthly DMS data -------- 
    disp(' ')
    disp(['Main program, for domain ', domain])
    disp(['Looping over all times in original wrfchemi, interpolate DMS data for each month, create and write wrfchemi_DMS file'])

    for month_index = 1:length(months_list_m)
      % months_list_m contains all the individual months (e.g., for a run from 2012/12/1 to 2013/3/1, months_list_m = [12 1 2 3])
      month_i = months_list_m(month_index);
      year_i = years_list_m(month_index); % Corresponding year for this month (month_i)
      disp([ '  Month ' num2str(month_i)]);

      % Regrid the monthly DMS oceanic values from Lana et al
      % This is only done once per month because the data is monthly
      disp(['  Regridding DMS'])
      DMS_FILE = [LANA_DMS_DIRECTORY '/DMSclim_' upper(datestr(datenum(1, month_i, 1), 'mmm')) '.csv'];
      disp(['    Opening ' DMS_FILE])
      fid = fopen(DMS_FILE, 'r');
      DMS_string_format = repmat('%f', 1, 360);
      C = textscan(fid, DMS_string_format, 'delimiter', ',');
      data_DMS = flipud((cell2mat(C)))' * 1E-6; % Original data in µmol m-3, converted to mol m-3
      fclose(fid);
      clear C
      data_wrf_DMS = zeros(imax, jmax, zdim);
      for ii = 1:imax
          for jj = 1:jmax
              % If the grid cell is in the ocean, regrid oceanic DMS to the WRF grid by bilinear interpolation
              if(~data_landmask_wrf(ii, jj))
                  data_wrf_DMS(ii, jj, 1) = interp2(emissions_lat_DMS, emissions_lon_DMS, data_DMS, ...
                                           wrf_xlat(ii, jj), wrf_xlong(ii, jj), 'linear');
              end
          end
      end

      % After interpolation, some values of oceanic DMS near coasts are still
      % NaN because of the coarse resolution of the climatology.  If coastal
      % values are NaN, replace them with nearby values until all ocean
      % gridpoints are filled. This is not a very elegant way to do this, but
      % you only need to regrid the DMS oceanic content once per month so it
      % does not matter if this is a bit slow
      is_filling = 1;
      while(is_filling)
        data_wrf_DMS_2 = data_wrf_DMS;
        data_wrf_DMS_2(data_wrf_DMS == 0) = NaN;
        is_filling = 0;

        for ii = 1:imax
          for jj = 1:jmax
            if(isnan(data_wrf_DMS(ii, jj, 1)))
              i1_average = ii - 1;
              i2_average = ii + 1;
              j1_average = jj - 1;
              j2_average = jj + 1;

              if(ii == 1)
                i1_average = 1;
              end

              if(ii == imax)
                i2_average = imax;
              end

              if(jj == 1)
                j1_average = 1;
              end

              if(jj == jmax)
                j2_average = jmax;
              end

              if(any(any(~isnan(data_wrf_DMS_2(i1_average:i2_average, j1_average:j2_average)))))
                data_wrf_DMS(ii, jj, 1) = nan_mean(nan_mean(data_wrf_DMS_2(i1_average:i2_average, j1_average:j2_average)));
                is_filling = 1;
              end

            end % if (isnan(data_wrf_DMS(ii, jj, 1)))
          end % for jj
        end % for ii
      end % while is_filling

      data_wrf_DMS(isnan(data_wrf_DMS)) = 0;

      % -------- Loop over days and times, create wrfchemi, write DMS data to wrfchemi -------- 
      days_in_month = days_list_d(months_list_d == month_i & years_list_d == year_i);

      for day_index = 1:length(days_in_month)
        day_i = days_in_month(day_index);

        for ihour = 0:0 % 1 emission file per day. Change to ihour = 0:23 to get hourly files
          % Date
          date_wrfchemi = datenum(year_i, month_i, day_i, ihour, 0, 0); % Date in matlab format

          % Create & open the wrfchemi netcdf
          date_wrf = datestr(date_wrfchemi, 'yyyy-mm-dd_HH:MM:SS'); % Date in WRF format
          disp(['    Create wrfchemi for date ' datestr(date_wrfchemi, 'yyyy/mm/dd-HH:MM:SS')])
          wrfchemi_filename = ['wrfchemi_', domain, '_', date_wrf, '_DMS_OC'];
          create_netcdffile_mozartmosaic_dms(wrfchemi_filename, SPECNAMES_WRFCHEMI, date_wrf, zdim, imax, jmax, ...
                                   dx, dy, ref_lat, ref_lon, truelat1, truelat2, moad_cen_lat, map_proj, wrf_proj.mminlu,...
                                   'Created by Louis Marelle for WRF V4.0');
          ncid_wrfchemi = netcdf.open(wrfchemi_filename, 'NC_WRITE');

          data_wrf = data_wrf_DMS;

          % Variable name in the wrfchemi file
          wrfchemi_varname = ['E_', SPECNAMES_WRFCHEMI{1}];

          % Write regridded emissions (data_wrf) to the wrfchemi netcdf variable (E_...)
          ncwrite(ncid_wrfchemi, data_wrf, wrfchemi_varname);

          % Close netcdf
          netcdf.close(ncid_wrfchemi);

          % Using the NCO function ncks to add the E_DMS_OC variable to the already existing wrfchemi_<domain>* files
          wrfchemi_dms_filename = wrfchemi_filename;
          wrfchemi_destination_filename = [WRFCHEMI_PATH, '/wrfchemi_', domain, '_', date_wrf];
          disp(['    Adding ', wrfchemi_varname, ' to wrfchemi file ', wrfchemi_destination_filename])
          system(['ncks -v ', wrfchemi_varname, ' ', wrfchemi_dms_filename, ' -A ', wrfchemi_destination_filename, ';']);
          system(['rm -f ', wrfchemi_dms_filename]);

        end % for ihour (hourly loop)
      end % for day_index (daily loop)
    end % for month_index (monthly loop)
  end % for idomain (WRF domain loop)

  disp(' ')
  disp('-------- prep_wrfchemi_oceanic_dms done, everything OK --------')
  disp(' ')

% end % end main matlab function


%-------- FUNCTIONS --------

function [data_var] = ncread(ncid, VARNAME)
  %-------- Read netcdf variable VARNAME --------
  varid = netcdf.inqVarID(ncid, VARNAME);
  data_var = netcdf.getVar(ncid, varid);
end


function [] = ncwrite(ncid, data_var, VARNAME)
  %-------- Write data_var to netcdf variable VARNAME --------
  varid = netcdf.inqVarID(ncid, VARNAME);
  netcdf.putVar(ncid, varid, data_var);
end



