function [] = prep_anthro_emissions_mozcart_eclipse_rcp(startdate, enddate)

  %-------- Create WRF-Chem emissions (wrfchemi_d0*) for the MOZCART mechanism --------------%
  %
  % Louis Marelle, 2022/04/20
  %
  % Purpose:
  %   This MATLAB routine (function) creates daily wrfchemi_dXX input emission
  %   files for WRFChem using:
  %   - Anthropogenic emissions from the ECLIPSEv5 inventory
  %   - Shipping emissions from RCP8.5 (just change the directories/filenames 
  %     to use another RCP inventory)
  %
  %
  % Inputs: 
  %  startdate: first date at which the emissions must be created, in MATLAB datenum format
  %  enddate: last date at which the emissions must be created, in MATLAB datenum format
  %
  % Outputs:
  %   Creates daily wrfchemi files in the path where this routine was run
  %
  % This routine can be modified easily to create hourly files, by replacing the 
  % "for ihour = 0:0" loop by "for ihour = 0:23", and commenting the line:
  % "hourly_factors_now = 1". There is no daily/hourly variation for shipping 
  % emisssions.
  %
  % The routine writes wrfchemi files in the path where the routine is located/run. 
  % These files can take a lot of disk space for big runs
  %
  % Always check the emissions before running WRFChem to make sure that the routine
  % produced reasonable results (I usually check one VOC, one trace gas and one 
  % aerosol). In the future I might print some diagnosis at the end of this routine
  % to make sure that emission mass is conserved and the regridding or netcdf writing
  % goes well.
  %------------------------------------------------------------------------------------------%


  %-------- Input --------
  % This needs to contain the wrfinput file, can be different from the WRF run
  RUN_DIRECTORY = './';
  max_domains = 1;


  %-------- Parameters --------
  % Emission datasets directories
  ECLIPSE_DIRECTORY = '/data/marelle/EMISSIONS/ECLIPSE_v5/';
  RCP_DIRECTORY = '/data/marelle/EMISSIONS/ships/RCP/';

  % Years present in the ECLIPSE and RCP inventories
  years_avail_ECLIPSE = [2000 2005 2010 2020 2030 2040 2050 2060 2070 2080 2090 2100];
  years_avail_RCP = [2005 2010 2030 2050];
  
  %---- Species
  % Mechanism species, written in wrfchemi file
  SPECNAMES_MOZCART = {'CO', 'NH3', 'SO2', 'NO', 'NO2', 'BIGALK', 'BIGENE', 'C2H4', 'C2H5OH', 'C2H6',...
                       'C3H6', 'C3H8', 'CH2O', 'CH3CHO', 'CH3COCH3', 'CH3OH', 'MEK', 'TOLUENE', 'OC', 'BC'};
  % Emission datasets species
  SPECNAMES_ECLIPSE = {'CO', 'NH3', 'SO2', 'NOx', 'NOx', 'VOC', 'VOC', 'VOC', 'VOC', 'VOC', 'VOC', 'VOC', 'VOC', 'VOC', 'VOC', 'VOC', 'VOC', 'VOC', 'OC', 'BC'};
  ECLIPSE_SECTORS = {'ene', 'ind', 'dom', 'tra', 'agr', 'wst', 'slv'};
  SPECNAMES_RCP = {'CO', 'NH3', 'SO2', 'NO', 'NO', 'NMVOC', 'NMVOC', 'NMVOC', 'NMVOC', 'NMVOC', 'NMVOC', 'NMVOC', 'NMVOC', 'NMVOC', 'NMVOC', 'NMVOC', 'NMVOC', 'NMVOC', 'OC', 'BC'};

  % Molarweights for mechanism species, NOx = NO2 in inventory
  MOLARWEIGHTS_ECLIPSE = [28, 17, 64, 46, 46, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN];
  MOLARWEIGHTS_RCP = [28, 17, 64, 30, 30, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN];
  
  % NOx proportion as NO
  % Proportion of NOx stored as NO for ship emissions %EPA, 2000, von glasow et
  % al for full ref
  nox_as_no_ships = 0.94;
  % Proportion of NOx stored as NO for anthropogenic emissions % ~Kaynak et al
  % 2009, Finlayson-Pitts and Pitts, 1999, google "nox "90% no" "10% no2" for
  % more references
  nox_as_no_anthro = 0.90;

  % Number of vertical levels in the wrfchemi file
  ZDIM = 1;

  %---- Time variation factors
  % Daily cycle - hourly factors
  % Factors to apply a daily emission cycle: 1h to 24h, from TNO-MACC (Denier van
  % der Gon et al., full ref in Marelle et al., ACP, 2016)
  hourly_factors = [0.79, 0.72, 0.72, 0.71, 0.74, 0.80, 0.92, 1.08, 1.19, 1.22, 1.21, 1.21,...
                    1.17, 1.15, 1.14, 1.13, 1.10, 1.07, 1.04, 1.02, 1.02, 1.01, 0.96, 0.88;... ENERGY
                    0.75, 0.75, 0.78, 0.82, 0.88, 0.95, 1.02, 1.09, 1.16, 1.22, 1.28, 1.30,...
                    1.22, 1.24, 1.25, 1.16, 1.08, 1.01, 0.95, 0.90, 0.85, 0.81, 0.78, 0.75;... INDUSTRY
                    0.40, 0.40, 0.40, 0.40, 0.40, 0.50, 1.20, 1.50, 1.60, 1.60, 1.40, 1.20,...
                    1.10, 1.10, 1.00, 1.00, 1.00, 1.10, 1.40, 1.50, 1.40, 1.40, 1.00, 0.40;... RESIDENTIAL
                    0.19, 0.09, 0.06, 0.05, 0.09, 0.22, 0.86, 1.84, 1.86, 1.41, 1.24, 1.20,...
                    1.32, 1.44, 1.45, 1.59, 2.03, 2.08, 1.51, 1.06, 0.74, 0.62, 0.61, 0.44;... TRANSPORT
                    0.60, 0.60, 0.60, 0.60, 0.60, 0.65, 0.75, 0.90, 1.10, 1.25, 1.45, 1.60,...
                    1.80, 1.75, 1.70, 1.55, 1.35, 1.10, 0.90, 0.75, 0.65, 0.60, 0.60, 0.60;... AGRICULTURE
                    1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,...
                    1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00;... WASTE = 1.00
                    1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,...
                    1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00]; % SOLVENTS = 1.00
  % Weekly cycle - daily factors
  % Daily factors to apply a weekly emission cycle, Sunday to Saturday from TNO-MACC
  daily_factors = [0.85, 1.06, 1.06, 1.06, 1.06, 1.06, 0.85;... ENERGY
                   0.80, 1.08, 1.08, 1.08, 1.08, 1.08, 0.80;... INDUSTRY
                   0.80, 1.08, 1.08, 1.08, 1.08, 1.08, 0.80;... RESIDENTIAL
                   0.79, 1.02, 1.06, 1.08, 1.10, 1.14, 0.81;... TRANSPORT
                   1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00;... AGRICULTURE
                   1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00;... WASTE = 1.00
                   1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00]; % SOLVENTS = 1.00
  

  %-------- Initialization --------
  routine_filename = mfilename;
  disp(' ')
  disp(['-------- ', routine_filename, ' - produce wrfchemi files --------'])
  disp(' ')
  disp('Initializing')

  % Include my MATLAB toolbox (ijll and llij projection routines,
  % create_netcdf and speciate_voc functions)
  addpath('/home/marelle/scripts/matlab/TOOLBOX');

  % Cleanup existing wrfchemi files
  system(['rm -f wrfchemi_*']);

  % Year for the ECLIPSE and RCP inventory
  year_ECLIPSE = year(startdate);
  year_RCP = year(startdate);

  %---- Initialization: WRF emission dates
  wrfchemi_dates = startdate:enddate;
  ndates = length(wrfchemi_dates);
  years_list = year(wrfchemi_dates);
  months_list = month(wrfchemi_dates);
  days_list = day(wrfchemi_dates);
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

  %---- Initialization: emission inventory years
  % Find the closest years to year_ECLIPSE and year_RCP in the ECLIPSE and RCP inventories
  if(year_ECLIPSE < min(years_avail_ECLIPSE) || year_ECLIPSE >= max(years_avail_ECLIPSE))
    error('Error, year_ECLIPSE is out of range')
  end
  if(year_RCP < min(years_avail_RCP) || year_RCP >= max(years_avail_RCP))
    error('Error, year_RCP is out of range')
  end
  % Find the closest year in the list. If the requested year is after the
  % closest year, pick the closest year and the next year in the list for
  % interpolation. If the requested year is before the closest year, pick the
  % closest year and the previous year in the list for interpolation
  [min_val_ECLIPSE, indx_closest_ECLIPSE_year] = min(abs(year_ECLIPSE - years_avail_ECLIPSE));
  [min_val_RCP, indx_closest_RCP_year] = min(abs(year_RCP - years_avail_RCP));
  if(min_val_ECLIPSE > 0)
    % Closest year is before
    years_closest_ECLIPSE = years_avail_ECLIPSE(indx_closest_ECLIPSE_year:indx_closest_ECLIPSE_year + 1);
  elseif(min_val_ECLIPSE < 0)
    % Closest year is after
    years_closest_ECLIPSE = years_avail_ECLIPSE(indx_closest_ECLIPSE_year - 1:indx_closest_ECLIPSE_year);
  else
    % Closest year is in the list
    years_closest_ECLIPSE = [year_ECLIPSE year_ECLIPSE];
  end
  if(min_val_RCP > 0)
    % Closest year is before
    years_closest_RCP = years_avail_RCP(indx_closest_RCP_year:indx_closest_RCP_year + 1);
  elseif(min_val_RCP < 0)
    % Closest year is after
    years_closest_RCP = years_avail_RCP(indx_closest_RCP_year - 1:indx_closest_RCP_year);
  else
    % Closest year is in the list
    years_closest_RCP = [year_RCP year_RCP];
  end
  % Compute the factors for the linear interpolation of emissions to year_ECLIPSE
  if(any(~(years_avail_ECLIPSE - year_ECLIPSE)))
    % If the requested year is available in the list, no need to interpolate
    interp_factor_before_ECLIPSE = 1;
    interp_factor_after_ECLIPSE = 0; 
  else
    % If the year is not in the list, need to interpolate between closest years
    interp_factor_before_ECLIPSE = ((year_ECLIPSE) - (years_closest_ECLIPSE(1))) ...
                          / ((years_closest_ECLIPSE(2)) - (years_closest_ECLIPSE(1)));
    interp_factor_after_ECLIPSE = ((years_closest_ECLIPSE(2)) - (year_ECLIPSE)) ... 
                          / ((years_closest_ECLIPSE(2)) - (years_closest_ECLIPSE(1)));
  end
  % Same for year_RCP
  if(any(~(years_avail_RCP - year_RCP)))
    % If the requested year is available in the list, no need to interpolate
    years_closest_RCP = [year_RCP year_RCP];
    interp_factor_before_RCP = 1;
    interp_factor_after_RCP = 0;
  else
    % If it is not in the list, need to interpolate
    interp_factor_before_RCP = ((year_RCP) - (years_closest_RCP(1))) ... 
                    / ((years_closest_RCP(2)) - (years_closest_RCP(1)));
    interp_factor_after_RCP = ((years_closest_RCP(2)) - (year_RCP)) ...
                    / ((years_closest_RCP(2)) - (years_closest_RCP(1)));
  end
  disp(['  Year is ', num2str(year_ECLIPSE), ', interpolating between closest years for ECLIPSE, ' num2str(years_closest_ECLIPSE)])
  disp(['  Year is ', num2str(year_RCP), ', interpolating between closest years for RCP, ' num2str(years_closest_RCP)])

  %---- Initialization: get emission grid dimensions
  %-- Get ECLIPSE grid dimensions
  disp('  Reading lat & lon from ECLIPSE emission files')
  eclipse_filename = [ECLIPSE_DIRECTORY '/ETP_base_CLE_V5_CO_' num2str(years_closest_ECLIPSE(1)) '.nc'];
  ncid = netcdf.open(eclipse_filename, 'NC_NOWRITE');
  varid = netcdf.inqVarID(ncid, 'lat');
  emissions_lat_ECLIPSE = netcdf.getVar(ncid, varid);
  varid = netcdf.inqVarID(ncid, 'lon');
  emissions_lon_ECLIPSE = netcdf.getVar(ncid, varid);
  netcdf.close(ncid);
  % Dimensions of ECLIPSE grid
  DLAT_ECLIPSE = emissions_lat_ECLIPSE(2) - emissions_lat_ECLIPSE(1);
  DLON_ECLIPSE = emissions_lon_ECLIPSE(2) - emissions_lon_ECLIPSE(1);
  imax_emiss_ECLIPSE = length(emissions_lon_ECLIPSE);
  jmax_emiss_ECLIPSE = length(emissions_lat_ECLIPSE);
  %-- Get RCP grid dimensions
  disp('  Reading lat & lon from RCP emission files')
  ncid = netcdf.open([ RCP_DIRECTORY '/RCP85_CO_ships_' num2str(years_closest_RCP(1)) '_0.5x0.5_v1_21_12_2009.nc'], 'NC_NOWRITE');
  varid = netcdf.inqVarID(ncid, 'lat');
  emissions_lat_RCP = netcdf.getVar(ncid, varid);
  varid = netcdf.inqVarID(ncid, 'lon');
  emissions_lon_RCP = netcdf.getVar(ncid, varid);
  netcdf.close(ncid);
  DLAT_RCP = emissions_lat_RCP(2) - emissions_lat_RCP(1);
  DLON_RCP = emissions_lon_RCP(2) - emissions_lon_RCP(1);
  imax_emiss_RCP = length(emissions_lon_RCP);
  jmax_emiss_RCP = length(emissions_lat_RCP);
  % Convert RCP lon from 0-360 to -180 to 180
  emissions_lon_RCP = mod((emissions_lon_RCP + 180), 360) - 180;

  %---- Initialization: Check dimensions
  nspecies = length(SPECNAMES_MOZCART);
  nsectors_ECLIPSE = length(ECLIPSE_SECTORS);
  if(nspecies ~= length(SPECNAMES_ECLIPSE) | nspecies ~= length(SPECNAMES_RCP) ...
      | nspecies ~= length(MOLARWEIGHTS_ECLIPSE) | nspecies ~= length(MOLARWEIGHTS_RCP))
    error(['Error, nspecies should be the same for all inventories and for the mechanism, ' num2str([nspecies, length(SPECNAMES_ECLIPSE), length(SPECNAMES_RCP), length(MOLARWEIGHTS_ECLIPSE), length(MOLARWEIGHTS_RCP)])])
  end


  %-------- Main prepemis_anthro program  --------
  % WRF domains loop
  for idomain = 1:max_domains
    % Domain name
    domain = ['d0', num2str(idomain)];
    disp(' ')
    disp(['WRF Domain ', domain])
    disp(['  Initializing for WRF domain ', domain])

    %---- Initializing grid for domain idomain
    % Read projection & grid info from wrfinput file
    disp('  Reading WRF domain info')
    wrfinput_file = [RUN_DIRECTORY, '/wrfinput_', domain];
    [moad_cen_lat, truelat1, truelat2, stdlon, imax, jmax, kmax, dx, dy, ref_lat,...
      ref_lon, map_proj, hemi, ref_x, ref_y] = get_WRF_grid(wrfinput_file);
    wrf_proj = get_WRF_proj(wrfinput_file);
    % Find min & max lat & lon from the WRF domain
    disp('  Opening wrfinput file to get WRF lat & lon boundaries')
    ncid = netcdf.open(wrfinput_file, 'NC_NOWRITE');
    data_lat_wrf = ncread(ncid, 'XLAT');
    data_lon_wrf = ncread(ncid, 'XLONG');
    % Adding -1 & + 1 degrees to avoid missing emissions at the domain
    % boundaries
    lat_min = min(min(data_lat_wrf)) - 1;
    lat_max = max(max(data_lat_wrf)) + 1;
    lon_min = min(min(data_lon_wrf)) - 1;
    lon_max = max(max(data_lon_wrf)) + 1;
    % Read the map factors from the wrfinput file
    data_mapfacx_wrf = ncread(ncid, 'MAPFAC_M');
    data_mapfacy_wrf = ncread(ncid, 'MAPFAC_M');
    data_landmask_wrf = ncread(ncid, 'LANDMASK');
    netcdf.close(ncid);
    % Choose ncells, number of minicells per grid length, used for regridding emissions 
    % This was determined by trial and error, and might not work in all cases, check the 
    % emissions for spurious "lines" or "waves" through the data, which indicate a too low
    % minicell number was chosen
    if(min(dx, dy) >= 100)
      ncells_ECLIPSE = 8;
      ncells_RCP = 8;
    elseif(min(dx, dy) >= 50)
      ncells_ECLIPSE = 16;
      ncells_RCP = 16;
    elseif ( min(dx, dy) >= 25)
      ncells_ECLIPSE = 32;
      ncells_RCP = 32;
    elseif ( min(dx, dy) >= 10)
      ncells_ECLIPSE = 64;
      ncells_RCP = 64;
    else
      ncells_ECLIPSE = 128;
      ncells_RCP = 128;
    end % if min(dx, dy)
    disp(['  Minicell number for domain ', domain, ', ncells_ECLIPSE = ' num2str(ncells_ECLIPSE) ', ncells_RCP = ' num2str(ncells_RCP)])

    
    % -------- Emission mapping on wrf grid --------
    % The objective of this step is to create a mapping between the grid of
    % each inventory and the WRF grid
    % This mapping is later used to distribute emissions from each inventory on the WRF grid 
    disp(['  Emission mapping from emission inventory to WRF grid for domain ', domain])
  
    % Emission mapping - Creating the mapping  - ECLIPSE
    tic
    disp('    Creating the mapping, ECLIPSE')
    [mapping_array_ECLIPSE] = map_emissions(imax, jmax, imax_emiss_ECLIPSE, jmax_emiss_ECLIPSE,...
                                            emissions_lon_ECLIPSE, emissions_lat_ECLIPSE,...
                                            DLON_ECLIPSE, DLAT_ECLIPSE, lon_min, lon_max,...
                                            lat_min, lat_max, ncells_ECLIPSE, truelat1,...
                                            truelat2, hemi, stdlon, ref_lat, ref_lon,...
                                            ref_x, ref_y, dx, map_proj);
    t_elapsed = toc;
    disp(['      Elapsed time ' num2str(t_elapsed) 's']);
    
    % Emission mapping - Creating the mapping - RCP
    tic
    disp('    Creating the mapping, RCP')
    [mapping_array_RCP] = map_emissions(imax, jmax, imax_emiss_RCP, jmax_emiss_RCP, emissions_lon_RCP, emissions_lat_RCP, DLON_RCP, DLAT_RCP, lon_min, lon_max, lat_min, lat_max, ncells_RCP, truelat1, truelat2, hemi, stdlon, ref_lat, ref_lon, ref_x, ref_y, dx, map_proj);
    t_elapsed = toc;
    disp(['      Elapsed time ' num2str(t_elapsed) 's']);


    % -------- Loop over months and interpolate monthly emission data -------- 
    for imonth = 1:length(months_list_m) 
      % months_list_m contains all the individual months (e.g., for a run from
      % 2012/12/1 to 2013/3/1, months_list_m = [12 1 2 3])
      month_i = months_list_m(imonth);
      % Years corresponding to this month (month_i)
      year_i = years_list_m(imonth);
      disp(' ')
      disp([ 'Regridding monthly emission data to WRF domain ', domain, ' for ' datestr(datenum(year_i, month_i, 1), 'yyyy/mm')]);
      
      %---- Monthly scaling factors for ECLIPSE
      % Open ECLIPSE monthly factors for current month
      ECLIPSE_MONTHS_FILENAME = [ECLIPSE_DIRECTORY '/ECLIPSEv5_monthly_patterns.nc'];
      ncid = netcdf.open(ECLIPSE_MONTHS_FILENAME, 'NC_NOWRITE');
      eclipse_mfac_dom = ncread(ncid, 'dom');
      eclipse_mfac_dom = eclipse_mfac_dom(:, :, month_i);
      eclipse_mfac_ene = ncread(ncid, 'ene');
      eclipse_mfac_ene = eclipse_mfac_ene(:, :, month_i);
      eclipse_mfac_ind = ncread(ncid, 'ind');
      eclipse_mfac_ind = eclipse_mfac_ind(:, :, month_i);
      eclipse_mfac_tra = ncread(ncid, 'tra');
      eclipse_mfac_tra = eclipse_mfac_tra(:, :, month_i);
      eclipse_mfac_wst = ncread(ncid, 'other');
      eclipse_mfac_wst = eclipse_mfac_wst(:, :, month_i);
      eclipse_mfac_slv = ncread(ncid, 'other');
      eclipse_mfac_slv = eclipse_mfac_slv(:, :, month_i);
      eclipse_mfac_agr = ncread(ncid, 'agr');
      eclipse_mfac_agr = eclipse_mfac_agr(:, :, month_i);
      eclipse_mfac_agr_NH3 = ncread(ncid, 'agr_NH3');
      eclipse_mfac_agr_NH3 = eclipse_mfac_agr_NH3(:, :, month_i);
      netcdf.close(ncid);
      
      
      % -------- Regrid monthly VOCs --------
      disp(['  Regridding monthly NMVOC'])

      %---- Monthly ECLIPSE VOCs
      disp(['    Regridding ECLIPSE NMVOC'])
      % This is done only once per month (ECLIPSE emissions are monthly), for bulk VOCs
      % The speciation from bulk VOCs to specific WRFChem VOCs is done below and is very fast compared to this step

      %---- Monthly ECLIPSE VOCs: regrid for the closest year before year_ECLIPSE
      ECLIPSE_FILENAME = [ECLIPSE_DIRECTORY '/ETP_base_CLE_V5_VOC_' num2str(years_closest_ECLIPSE(1)) '.nc'];
      disp(['    Opening ' ECLIPSE_FILENAME])
      ncid = netcdf.open(ECLIPSE_FILENAME, 'NC_NOWRITE');
      % ECLIPSE VOCs year_before: Loop over sectors, regrid apply monthly factors
      for isector = 1:nsectors_ECLIPSE
        % Check if emissions exist for this sector
        errorID = 1;
        try
          data_ECLIPSE = ncread(ncid, ['emis_' ECLIPSE_SECTORS{isector}]);
          errorID = 0;
        catch exception
          errorID = 1;
        end
        if(~errorID)
          data_ECLIPSE(isnan(data_ECLIPSE)) = 0;
          eval(['eclipse_mfac = eclipse_mfac_' ECLIPSE_SECTORS{isector} ';']);
          data_wrf_ECLIPSE = regrid(data_ECLIPSE .* eclipse_mfac, mapping_array_ECLIPSE, ncells_ECLIPSE);
          data_wrf_ECLIPSE(isnan(data_wrf_ECLIPSE)) = 0;
          eval(['data_wrf_ECLIPSE_VOC_' ECLIPSE_SECTORS{isector} '_before = data_wrf_ECLIPSE;']);
        else
          %disp(['      ', ECLIPSE_FILENAME ' : ' ECLIPSE_SECTORS{isector}   ' does not exist'])
        end % if error
      end % for isector (ECLIPSE emission sectors)
      netcdf.close(ncid); % Close ECLIPSE emission file
      %---- Monthly ECLIPSE VOCs: regrid for the closest year after year_ECLIPSE
      ECLIPSE_FILENAME = [ECLIPSE_DIRECTORY '/ETP_base_CLE_V5_VOC_' num2str(years_closest_ECLIPSE(2)) '.nc'];
      disp(['    Opening ' ECLIPSE_FILENAME])
      ncid = netcdf.open(ECLIPSE_FILENAME, 'NC_NOWRITE');
      % ECLIPSE VOCs year_after: Loop over sectors, regrid apply monthly factors
      for isector = 1:nsectors_ECLIPSE
        % Check if emissions exist for this sector
        errorID = 1;
        try
          data_ECLIPSE = ncread(ncid, ['emis_' ECLIPSE_SECTORS{isector}]);
          errorID = 0;
        catch exception
          errorID = 1;
        end
        if(~errorID)
          data_ECLIPSE(isnan(data_ECLIPSE)) = 0;
          eval(['eclipse_mfac = eclipse_mfac_' ECLIPSE_SECTORS{isector} ';']);
          data_wrf_ECLIPSE = regrid(data_ECLIPSE .* eclipse_mfac, mapping_array_ECLIPSE, ncells_ECLIPSE);
          data_wrf_ECLIPSE(isnan(data_wrf_ECLIPSE)) = 0;
          eval(['data_wrf_ECLIPSE_VOC_' ECLIPSE_SECTORS{isector} '_after = data_wrf_ECLIPSE;']);
        else
          %disp(['      ', ECLIPSE_FILENAME ' : ' ECLIPSE_SECTORS{isector}   ' does not exist'])
        end % if error
      end % for isector (ECLIPSE emission sectors)
      netcdf.close(ncid); % Close ECLIPSE emission file
      %----- Monthly ECLIPSE VOCs: linear interpolation emissions to
      % year_ECLIPSE between year_before and year_after
      for isector = 1:nsectors_ECLIPSE
        if(exist(['data_wrf_ECLIPSE_VOC_', ECLIPSE_SECTORS{isector}, '_before'], 'var')...
            & exist(['data_wrf_ECLIPSE_VOC_', ECLIPSE_SECTORS{isector}, '_after'], 'var'))
          eval(['data_wrf_ECLIPSE_VOC_' ECLIPSE_SECTORS{isector} ' = data_wrf_ECLIPSE_VOC_' ECLIPSE_SECTORS{isector} ...
                '_after * interp_factor_after_ECLIPSE  + data_wrf_ECLIPSE_VOC_' ECLIPSE_SECTORS{isector} ...
                '_before * interp_factor_before_ECLIPSE;']);
        end % if exist (ECLIPSE regridded variables exist)
      end % for isector (ECLIPSE emission sectors)

      
      %---- Monthly RCP VOCs
      disp(['  Regridding RCP shipping NMVOCs'])
      % This is done only once per month (RCP emissions are monthly), for bulk VOCs
      % The speciation from bulk VOCs to specific WRFChem VOCs is done below and is very fast compared to this step
      % Monthly RCP VOCs: regrid for the closest year before year_RCP
      RCP_FILENAME = [ RCP_DIRECTORY '/RCP85_NMVOC_ships_' num2str(years_closest_RCP(1)) '_0.5x0.5_v1_21_12_2009.nc'];
      disp(['    Opening ' RCP_FILENAME])
      ncid = netcdf.open(RCP_FILENAME, 'NC_NOWRITE');
      data_RCP = ncread(ncid, 'emiss_shp');
      data_RCP(isnan(data_RCP)) = 0;
      netcdf.close(ncid);
      data_wrf_RCP_NMVOC = regrid_area(data_RCP(:, :, month_i), mapping_array_RCP);
      data_wrf_RCP_NMVOC(isnan(data_wrf_RCP_NMVOC)) = 0;
      data_wrf_RCP_NMVOC_before = data_wrf_RCP_NMVOC;
      % Monthly RCP VOCs: regrid for the closest year after year_RCP
      RCP_FILENAME = [RCP_DIRECTORY '/RCP85_NMVOC_ships_' num2str(years_closest_RCP(2)) '_0.5x0.5_v1_21_12_2009.nc'];
      disp(['    Opening ' RCP_FILENAME])
      ncid = netcdf.open(RCP_FILENAME, 'NC_NOWRITE');
      data_RCP = ncread(ncid, 'emiss_shp');
      data_RCP(isnan(data_RCP)) = 0;
      netcdf.close(ncid);
      data_wrf_RCP_NMVOC = regrid_area(data_RCP(:, :, month_i), mapping_array_RCP);
      data_wrf_RCP_NMVOC(isnan(data_wrf_RCP_NMVOC)) = 0;
      data_wrf_RCP_NMVOC_after = data_wrf_RCP_NMVOC;
      % Monthly RCP VOCs: linear interpolation emissions to
      data_wrf_RCP_NMVOC = data_wrf_RCP_NMVOC_before * interp_factor_before_RCP + data_wrf_RCP_NMVOC_after * interp_factor_after_RCP;

      
      %-------- Regrid remaining monthly (non-VOC) emissions --------
      disp(['  Regridding remaining ECLIPSE and RCP species'])
      for ispecies = 1:nspecies

        %-------- Monthly ECLIPSE non-VOC emissions regridding
        specname_ECLIPSE = SPECNAMES_ECLIPSE{ispecies};
        if(specname_ECLIPSE & ~strcmp(specname_ECLIPSE,'VOC'))
          disp(['    Regridding ECLIPSE emissions for ' specname_ECLIPSE ' : '])
          %---- Monthly ECLIPSE non-VOC: regrid species for the closest year before year_ECLIPSE
          ECLIPSE_FILENAME = [ECLIPSE_DIRECTORY '/ETP_base_CLE_V5_' specname_ECLIPSE '_' num2str(years_closest_ECLIPSE(1)) '.nc'];
          disp(['      Opening ' ECLIPSE_FILENAME])
          ncid = netcdf.open(ECLIPSE_FILENAME, 'NC_NOWRITE');
          % Loop over ECLIPSE sectors
          for isector = 1:nsectors_ECLIPSE
            errorID = 1;
            try
              data_ECLIPSE = ncread(ncid, ['emis_' ECLIPSE_SECTORS{isector}]);
              errorID = 0;
            catch exception
              errorID = 1;
            end
            if(~errorID)
              data_ECLIPSE(isnan(data_ECLIPSE)) = 0;
              eval(['eclipse_mfac = eclipse_mfac_' ECLIPSE_SECTORS{isector} ';']);
              if(strcmp(specname_ECLIPSE, 'NH3') && strcmp(ECLIPSE_SECTORS{isector}, 'agr'))
                eclipse_mfac = eclipse_mfac_agr_NH3;
              end
              data_wrf_ECLIPSE = regrid(data_ECLIPSE .* eclipse_mfac, mapping_array_ECLIPSE, ncells_ECLIPSE);
              data_wrf_ECLIPSE(isnan(data_wrf_ECLIPSE)) = 0;
              eval(['data_wrf_ECLIPSE_' specname_ECLIPSE '_' ECLIPSE_SECTORS{isector} '_before = data_wrf_ECLIPSE;']);
            else
              % disp(['        ', ECLIPSE_FILENAME ' : ' ECLIPSE_SECTORS{isector}   ' does not exist'])
            end % if error
          end % for isector (ECLIPSE sectors)
          netcdf.close(ncid); % Close ECLIPSE emission file
          %---- Monthly ECLIPSE non-VOC: regrid species for the closest year after year_ECLIPSE
          ECLIPSE_FILENAME = [ECLIPSE_DIRECTORY '/ETP_base_CLE_V5_' specname_ECLIPSE '_' num2str(years_closest_ECLIPSE(2)) '.nc'];
          disp(['      Opening ' ECLIPSE_FILENAME])
          ncid = netcdf.open(ECLIPSE_FILENAME, 'NC_NOWRITE');
          for isector = 1:nsectors_ECLIPSE
            errorID = 1;
            try
              data_ECLIPSE = ncread(ncid, ['emis_' ECLIPSE_SECTORS{isector}]);
              errorID = 0;
            catch exception
              errorID = 1;
            end
            if(~errorID)
              data_ECLIPSE(isnan(data_ECLIPSE)) = 0;
              eval(['eclipse_mfac = eclipse_mfac_' ECLIPSE_SECTORS{isector} ';']);
              if(strcmp(specname_ECLIPSE, 'NH3') && strcmp(ECLIPSE_SECTORS{isector}, 'agr'))
                  eclipse_mfac = eclipse_mfac_agr_NH3;
              end
              data_wrf_ECLIPSE = regrid(data_ECLIPSE .* eclipse_mfac, mapping_array_ECLIPSE, ncells_ECLIPSE);
              data_wrf_ECLIPSE(isnan(data_wrf_ECLIPSE)) = 0;
              eval(['data_wrf_ECLIPSE_' specname_ECLIPSE '_' ECLIPSE_SECTORS{isector} '_after = data_wrf_ECLIPSE;']);
            else
              %disp(['        ', ECLIPSE_FILENAME ' : ' ECLIPSE_SECTORS{isector}   ' does not exist'])
            end % if error
          end % for isector (ECLIPSE emission sectors)
          netcdf.close(ncid); % Close ECLIPSE emission file
          %----- Monthly ECLIPSE non-VOC: Linear time interpolation of species ispecies in ECLIPSE to year_ECLIPSE
          for isector = 1:nsectors_ECLIPSE
            if(exist(['data_wrf_ECLIPSE_' specname_ECLIPSE '_', ECLIPSE_SECTORS{isector}, '_before'], 'var')...
                & exist(['data_wrf_ECLIPSE_' specname_ECLIPSE '_', ECLIPSE_SECTORS{isector}, '_after'], 'var'))
              % Linear interpolation of ECLIPSE emission data to year_ECLIPSE
              eval(['data_wrf_ECLIPSE_' specname_ECLIPSE '_' ECLIPSE_SECTORS{isector}...
                    ' = data_wrf_ECLIPSE_' specname_ECLIPSE '_' ECLIPSE_SECTORS{isector} ...
                    '_after * interp_factor_after_ECLIPSE  + data_wrf_ECLIPSE_' ...
                    specname_ECLIPSE '_' ECLIPSE_SECTORS{isector} ...
                    '_before * interp_factor_before_ECLIPSE;']);
            end % if exist (ECLIPSE regridded variables exist)
          end % for isector (ECLIPSE emission sectors)
        end % if(specname_ECLIPSE & ~strcmp(specname_ECLIPSE,'VOC'))

        %---- Monthly RCP non-VOC emissions regridding
        specname_RCP = SPECNAMES_RCP{ispecies};
        if(specname_RCP & ~strcmp(specname_RCP, 'NMVOC'))
          disp(['    Regridding RCP emissions for ' specname_RCP ' : '])
          % Monthly RCP non-VOC: regrid species for the closest year before year_RCP
          RCP_FILENAME = [ RCP_DIRECTORY '/RCP85_' specname_RCP '_ships_' num2str(years_closest_RCP(1)) '_0.5x0.5_v1_21_12_2009.nc'];
          disp(['      Opening ' RCP_FILENAME])
          ncid = netcdf.open(RCP_FILENAME, 'NC_NOWRITE');
          data_RCP = ncread(ncid, 'emiss_shp');
          data_RCP(isnan(data_RCP)) = 0;
          data_wrf_RCP = regrid_area(data_RCP(:, :, month_i), mapping_array_RCP);
          eval(['data_wrf_RCP_' specname_RCP '_before = data_wrf_RCP;']);
          netcdf.close(ncid);
          % Monthly RCP non-VOC: regrid species for the closest year after year_RCP
          RCP_FILENAME = [ RCP_DIRECTORY '/RCP85_' specname_RCP '_ships_' num2str(years_closest_RCP(2)) '_0.5x0.5_v1_21_12_2009.nc'];
          disp(['      Opening ' RCP_FILENAME])
          ncid = netcdf.open(RCP_FILENAME, 'NC_NOWRITE');
          data_RCP = ncread(ncid, 'emiss_shp');
          data_RCP(isnan(data_RCP)) = 0;
          data_wrf_RCP = regrid_area(data_RCP(:, :, month_i), mapping_array_RCP);
          eval(['data_wrf_RCP_' specname_RCP '_after = data_wrf_RCP;']);
          netcdf.close(ncid);
          % Monthly RCP non-VOC: Linear time interpolation of species ispecies in RCP to year_RCP
          eval(['data_wrf_RCP_' specname_RCP ' = data_wrf_RCP_' specname_RCP...
                '_before * interp_factor_before_RCP +  data_wrf_RCP_' specname_RCP...
                '_after * interp_factor_after_RCP;']);
        end % if(specname_RCP & ~strcmp(specname_RCP, 'NMVOC'))

      end % for ispecies
      
      
      %-------- Create wrfchemi, assemble emissions, write to file --------
      days_in_month = days_list_d(months_list_d == month_i & years_list_d == year_i); 
      for iday = 1:length(days_in_month)
        day_i = days_in_month(iday);

        % Hourly loop, 1 emission file per day. Change to ihour = 0:23 to get
        % hourly files, also comment hourly_factors_now = 1.
        for ihour = 0:0
          %---- Set exact time
          % Date in matlab format
          wrfchemi_date = datenum(year_i, month_i, day_i, ihour, 0, 0);
          % Calculate day of week for applying daily emission variation factors
          dayofweek_now = weekday(wrfchemi_date);
          dayofweek_now(dayofweek_now == 0) = 1;
          % Calculate the local "solar" time in order to add hourly variations if needed
          local_time = round(mod(ihour + data_lon_wrf / 180 * 12, 24));
          local_time(local_time == 0) = 24;
          % Date in WRF format
          date_wrf = datestr(wrfchemi_date, 'yyyy-mm-dd_HH:MM:SS');

          %---- Create wrfchemi and open it
          disp(['  Create wrfchemi_', domain, ' for date ', datestr(wrfchemi_date, 'yyyy/mm/dd-HH:MM:SS')])
          wrfchemi_filename = ['wrfchemi_', domain, '_', date_wrf];
          create_netcdffile_mozcart(wrfchemi_filename, SPECNAMES_MOZCART,...
                                    date_wrf, ZDIM, imax, jmax, ...
                                    dx, dy, ref_lat, ref_lon, truelat1,...
                                    truelat2, moad_cen_lat, map_proj, wrf_proj.mminlu,...
                                    ['Created by Louis Marelle for WRF V4.0 ',...
                                    'using ECLIPSE_v5 emission files']);
          ncid_wrfchemi = netcdf.open(wrfchemi_filename, 'NC_WRITE');


          %-------- Final emission assembly --------- 
          % Emissions speciation, time variation and unit conversion ; writing
          % to wrfchemi
          for ispecies = 1:nspecies
            specname_mechanism = SPECNAMES_MOZCART{ispecies};
            disp(['    Retrieve data and write to wrfchemi, ', specname_mechanism])

            %---- Assemble emissions: ECLIPSE emissions
            specname_ECLIPSE = SPECNAMES_ECLIPSE{ispecies};
            data_wrf_ECLIPSE = zeros(imax, jmax);
            if(specname_ECLIPSE)
              if(strcmp(specname_ECLIPSE, 'VOC'))
                %-- Assemble emissions: ECLIPSE VOCs
                % If species is a VOC, retrieve ECLIPSE bulk VOCs for
                % isector and speciate to WRFChem VOCs
                for isector = 1:nsectors_ECLIPSE
                  if(exist(['data_wrf_ECLIPSE_VOC_' ECLIPSE_SECTORS{isector}], 'var'))
                    eval(['data_wrf_ECLIPSE_VOC = data_wrf_ECLIPSE_VOC_' ECLIPSE_SECTORS{isector} ';']);
                    % Hour of day emission variation
                    hourly_factors_now = reshape(hourly_factors(isector, local_time), size(local_time));
                    %turn off the hourly factors for daily emissions
                    hourly_factors_now = 1;
                    % Day of week emission variation
                    daily_factor_now = daily_factors(isector, dayofweek_now);
                    % Convert from g(tot VOC)/cell/month to moles(spec VOC)/cell/month
                    data_wrf_ECLIPSE_VOC = voc_speciation_mozart(data_wrf_ECLIPSE_VOC * 1E9,...
                                             specname_mechanism, ECLIPSE_SECTORS{isector});
                    data_wrf_ECLIPSE = data_wrf_ECLIPSE + data_wrf_ECLIPSE_VOC...
                                       .* hourly_factors_now * daily_factor_now;
                  end
                end % for isector
              else 
                %-- Assemble emissions: ECLIPSE non-VOCs
                % If species is not a VOC, just retrieve the already regridded
                % emissions for this species and apply the time variation
                % factors
                for isector = 1:nsectors_ECLIPSE
                  if(exist(['data_wrf_ECLIPSE_' specname_ECLIPSE '_' ECLIPSE_SECTORS{isector}], 'var'))
                    eval(['data_wrf_ECLIPSE_SPECIES = data_wrf_ECLIPSE_' specname_ECLIPSE '_' ECLIPSE_SECTORS{isector} ';']);
                    % Hour of day emission variation
                    hourly_factors_now = reshape(hourly_factors(isector, local_time), size(local_time));
                    % Turn off the hourly factors for daily emissions
                    hourly_factors_now = 1;
                    % Day of week emission variation
                    daily_factor_now = daily_factors(isector, dayofweek_now);
                    data_wrf_ECLIPSE = data_wrf_ECLIPSE + data_wrf_ECLIPSE_SPECIES .* hourly_factors_now * daily_factor_now;
                  end
                end % for isector
              end % if(strcmp(specname_ECLIPSE, 'VOC'))
              %-- Assemble emissions: ECLIPSE speciation and unit conversion
              cell_area_km2 = dx ./ data_mapfacx_wrf .* dy ./ data_mapfacy_wrf;
              cell_area_m2 = cell_area_km2 * 1.0E6;
              hours_in_month = (datenum(0, month_i + 1, 0) - datenum(0, month_i, 0)) * 24;
              seconds_in_month = hours_in_month * 3600;
              if(strcmp(specname_ECLIPSE, 'BC'))
                % kT/month to ug/m2/s
                data_wrf_ECLIPSE = data_wrf_ECLIPSE * 1E15 ./ cell_area_m2 / seconds_in_month;
              elseif(strcmp(specname_ECLIPSE, 'OC'))
                % kT/month OC to ug/m2/s OC
                data_wrf_ECLIPSE = data_wrf_ECLIPSE * 1E15 ./ cell_area_m2 / seconds_in_month;
              elseif(strcmp(specname_ECLIPSE, 'VOC'))
                % moles/month to moles/km2/hour
                data_wrf_ECLIPSE = data_wrf_ECLIPSE ./ cell_area_km2 / hours_in_month;
              else
                % kT/month to moles/km2/hour
                data_wrf_ECLIPSE = data_wrf_ECLIPSE * 1E9 ./ cell_area_km2 / hours_in_month / MOLARWEIGHTS_ECLIPSE(ispecies);
              end
              %-- Assemble emissions: ECLIPSE NOx speciation to NO and NO2
              %TODO should this be specname_mechanism everywhere?
              if (strcmp(specname_mechanism, 'NO'))
                data_wrf_ECLIPSE = data_wrf_ECLIPSE * nox_as_no_anthro;
              elseif (strcmp(specname_mechanism, 'NO2'))
                data_wrf_ECLIPSE = data_wrf_ECLIPSE * (1 - nox_as_no_anthro);
              end
            end % if(specname_ECLIPSE)

            %---- Assemble emissions: RCP emissions
            specname_RCP = SPECNAMES_RCP{ispecies};
            data_wrf_RCP = zeros(imax, jmax);
            if(specname_RCP)
              if(strcmp(specname_RCP, 'NMVOC'))
                %-- Assemble emissions: RCP VOCs
                % If species ispecies is a VOC, retrieve RCP bulk VOCs for
                % isector and speciate to WRFChem VOCs
                data_wrf_RCP = data_wrf_RCP + data_wrf_RCP_NMVOC;
                data_wrf_RCP = voc_speciation_mozart(data_wrf_RCP * 1000, specname_mechanism, 'shp'); %g/m2/s to moles/m2/s
              else
                % If species ispecies is not a VOC, just retrieve the already regridded emissions 
                eval(['data_wrf_RCP_SPECIES = data_wrf_RCP_' specname_RCP ';']);
                data_wrf_RCP = data_wrf_RCP + data_wrf_RCP_SPECIES;
              end
              %-- Assemble emissions: RCP unit conversion (aerosols, VOCs, non-VOCs)
              if(strcmp(specname_RCP, 'BC'))
                % kg/m2/s to ug/m2/s
                data_wrf_RCP = data_wrf_RCP * 1E9;
              elseif(strcmp(specname_RCP, 'OC'))
                % kg/m2/s OC to ug/m2/s OC
                data_wrf_RCP = data_wrf_RCP * 1E9;
              elseif(strcmp(specname_RCP, 'NMVOC'))
                % moles/m2/s to moles/km2/hour
                data_wrf_RCP = data_wrf_RCP * 1E6 * 3600;
              else
                % kg/m2/s to moles/km2/hour
                data_wrf_RCP = data_wrf_RCP * 1000 / MOLARWEIGHTS_RCP(ispecies) * 1E6 * 3600;
              end
              %-- Assemble emissions: RCP NOx speciation to NO and NO2
              if (strcmp(specname_mechanism, 'NO'))
                  data_wrf_RCP = data_wrf_RCP * nox_as_no_ships;
              elseif (strcmp(specname_mechanism, 'NO2'))
                  data_wrf_RCP = data_wrf_RCP * (1 - nox_as_no_ships);
              end
            end % if(specname_RCP)
              
            %---- Assemble emissions: Write emissions to wfchemi file
            % Sum all the regridded emissions in data_wrf
            data_wrf = data_wrf_ECLIPSE + data_wrf_RCP;
            
            % Variable name in the wrfchemi file
            wrfchemi_varname = ['E_' specname_mechanism];
            % Write regridded emissions (data_wrf) to the wrfchemi netcdf
            % variable (E_...)
            ncwrite(ncid_wrfchemi, data_wrf, wrfchemi_varname);

          end % for ispecies

          % Close wrfchemi NetCDF
          netcdf.close(ncid_wrfchemi);

        end % for ihour
      end % for iday
    end % for imonth
  end % for idomain

  disp(' ')
  disp(['-------- ', routine_filename, ' - done --------'])
  disp(' ')

end % function



% -------- Additional sub-functions --------


function [] = ncwrite(ncid, data_var, VARNAME)
  % Write data_var to netcdf variable VARNAME
  varid = netcdf.inqVarID(ncid, VARNAME);
  netcdf.putVar(ncid, varid, data_var);
end


function [data_var] = ncread(ncid, VARNAME)
  % Read netcdf variable VARNAME
  varid = netcdf.inqVarID(ncid, VARNAME);
  data_var = netcdf.getVar(ncid, varid);
end


function [mapping_array] = map_emissions(imax, jmax, imax_emiss, jmax_emiss,...
                             emissions_lon, emissions_lat, DLON_emiss, DLAT_emiss,...
                             lon_min, lon_max, lat_min, lat_max, ncells, truelat1,...
                             truelat2, hemi, stdlon, ref_lat, ref_lon, ref_x,...
                             ref_y, dx, map_proj)
  % Create a "mapping" linking each point of the WRF grid to points of the grid
  % of the emission inventory
  %
  % To regrid emissions, we cut emission inventory cells into minicells
  % We then build for each inventory a table of correspondance linking each
  % inventory minicell to a WRF cell. This table (or mapping) has
  % the horizontal dimensions of the destination WRF grid (xmax, ymax). each
  % value mapping_array{x, y} contains a list
  % (i_cell_1, j_cell_1, i_cell_2, j_cell_2...) where i_cell_n, j_cell_n are
  % the coordinates of the emission inventory cells where the minicell n has
  % been cut

  mapping_array = cell(imax, jmax);
  for i = 1:imax_emiss
    for j = 1:jmax_emiss
      if(emissions_lat(j) > lat_min && emissions_lat(j) < lat_max && emissions_lon(i) > lon_min && emissions_lon(i) < lon_max)
        lat_cells = repmat(linspace(emissions_lat(j) - (DLAT_emiss - DLAT_emiss / ncells) / 2,...
          emissions_lat(j) + (DLAT_emiss - DLAT_emiss / ncells) / 2, ncells), 1, ncells);
        lon_cells = reshape(repmat(linspace(emissions_lon(i) - (DLON_emiss - DLON_emiss / ncells) / 2,...
          emissions_lon(i) + (DLON_emiss - DLON_emiss / ncells) / 2, ncells), ncells, 1), 1, ncells * ncells);
        [x_cells, y_cells] = llij(lat_cells, lon_cells, truelat1, truelat2, hemi, stdlon, ref_lat, ref_lon, ref_x, ref_y, dx, map_proj);
        x_cells = round(x_cells);
        y_cells = round(y_cells);
        for c = 1:ncells * ncells
          if(x_cells(c) > 0 && x_cells(c) < imax + 1 && y_cells(c) > 0 && y_cells(c) < jmax + 1)
            mapping_array{x_cells(c), y_cells(c)} = [mapping_array{x_cells(c), y_cells(c)}, i, j];
          end
        end
      end
    end
  end
end


function [data_wrf_emissions] = regrid_area(data_emissions, mapping_array)
  % Regrid emissions for an inventory with emissions in mass or moles /surface/time

  [imax, jmax] = size(mapping_array);
  data_wrf_emissions = zeros(imax, jmax);
  if(~isempty(data_emissions))
    data_wrf_n_emissions = zeros(imax, jmax);
    for i = 1:imax
      for j = 1:jmax
        if(~isempty(mapping_array{i, j}))
          a = mapping_array{i, j};
          for w = 1:2:length(a)
            % Surface emission data, summed (this sum has no physical meaning without the following average)
            data_wrf_emissions(i, j) = data_wrf_emissions(i, j) + data_emissions(a(w), a(w + 1));
            % Counter for each cell, used to calculate average surface emissions
            data_wrf_n_emissions(i, j) = data_wrf_n_emissions(i, j) + 1;
          end
        end
      end
    end
    % Divide emissions by the counter to calculate average emissions for each wrf grid
    data_wrf_n_emissions(data_wrf_n_emissions == 0) = 1;
    data_wrf_emissions = data_wrf_emissions ./ data_wrf_n_emissions;
  end
end


function [data_wrf_emissions] = regrid(data_emissions, mapping_array, ncells)
  % Regrid emissions for the inventoires with emissions in mass/time

  [imax, jmax] = size(mapping_array);
  data_wrf_emissions = zeros(imax, jmax);
  if(~isempty(data_emissions))
    for i = 1:imax
      for j = 1:jmax
        if(~isempty(mapping_array{i, j}))
          a = mapping_array{i, j};
          for w = 1:2:length(a)
            data_wrf_emissions(i, j) = data_wrf_emissions(i, j) + data_emissions(a(w), a(w + 1)) / (ncells * ncells);
          end
        end
      end
    end
  end
end

