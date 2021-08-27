import numpy as np
import os
import datetime
import pandas as pd
import netCDF4 as nc
from sklearn.linear_model import LinearRegression, HuberRegressor
import sympy as sym
import geopy.distance
import seaborn as sns
import h5py
import scipy.io

Workspace = '/Users/monicazhu/Documents/MATLAB/PAN_Data/Workspaces/SiteData'
wrf_path = os.path.join(Workspace, 'ML-WRF')
omi_path = os.path.join(Workspace, 'ML-OMI')
pasadena_wrf_path = os.path.join(Workspace, 'ML-WRF-Pasadena')
small_radius_wrf_path = os.path.join(Workspace, 'ML-WRF-Small-Radius')
small_radius_omi_path = os.path.join(Workspace, 'ML-OMI-Small-Radius')
obs_path = './Inputs'
output_path = './Outputs'
iso_output_path = os.path.join(output_path, 'ISO')
def_output_path = os.path.join(output_path, 'DEF')
sat_output_path = os.path.join(output_path, 'SAT')
sat_aks_output_path = os.path.join(output_path, 'SAT_AKS')
sat_aks_filter_output_path = os.path.join(output_path, 'SAT_AKS_FILTER')
sat_surr_aks_output_path = os.path.join(output_path, 'SAT_SURR_AKS')
sat_downwind_aks_output_path = os.path.join(output_path, 'SAT_DW_AKS')
oh_index_output_path = os.path.join(output_path, 'OH_Index')
oh_index_sr_output_path = os.path.join(output_path, 'OH_Index_SR')
sat_aks_sr_output_path = os.path.join(output_path, 'SAT_AKS_SR')
sat_iso_aks_sr_output_path = os.path.join(output_path, 'SAT_ISO_AKS_SR')
oh_index_mlg_output_path = os.path.join(output_path, 'OH_Index_MLG')

satellite_features = ['no2_vcd', 'hcho_vcd', 'o3', 'no2_vcd_aks', 'hcho_vcd_aks']

satellite_small_radius_features = ['no2_vcd', 'hcho_vcd', 'no2_vcd_aks', 'hcho_vcd_aks']

wrf_rename_satellite_features = ['wrf_no2_vcd', 'wrf_hcho_vcd', 'wrf_o3', 'wrf_no2_vcd_aks', 'wrf_hcho_vcd_aks']

full_feature_sets = ['no2', 'no', 'o3', 'hcho', 'pan', 'iso', 'gly',
                'act', 'temperature', 'PHOTR_O31D', 'PHOTR_CH2OR', 'o3_photolysis', 'hcho_photolysis',
             'QVAPOR', 'pressure','WindDir', 'WindSpeed', 'COSZEN', 'no2_vcd', 'no2_no', 'iso_vcd', 'hcho_vcd']

iso_feature_sets = ['o3', 'temperature',  'PHOTR_O31D', 'PHOTR_CH2OR',
                    'QVAPOR', 'pressure', 'COSZEN', 'no2_vcd', 'no2_no', 'hcho_vcd', 'iso_vcd', 'o3_photolysis']

oh_index_full_feature_sets = ['o3', 'temperature',  'PHOTR_O31D', 'PHOTR_CH2OR', 'jo1d_h2o',
                    'QVAPOR', 'pressure', 'COSZEN', 'no2', 'no2_no', 'hcho', 'iso', 'o3_photolysis', 'o3_jo1d']


oh_index_forward_feature_sets = ['PHOTR_O31D', 'QVAPOR', 'no2', 'iso', 'pressure', 'temperature', 'co', 'hcho',
                                 'PHOTR_CH2OR', 'o3', 'ch4', 'PHOTR_NO2', 'ndens']

oh_index_feature_sets = ['PHOTR_O31D', 'QVAPOR', 'no2', 'iso', 'pressure', 'temperature']


oh_index_mlg_feature_sets = ['PHOTR_O31D', 'QVAPOR', 'no2', 'iso', 'pressure', 'temperature']


oh_index_fillgap_feature_sets = ['PHOTR_O31D', 'QVAPOR', 'no2', 'iso', 'pressure', 'temperature', 'co', 'ch4']


opt_feature_sets = ['PHOTR_O31D', 'QVAPOR', 'no2_vcd', 'hcho_vcd', 'pressure', 'temperature', 'o3', 'PHOTR_CH2OR',
                    'COSZEN']

sat_feature_sets = ['PHOTR_O31D', 'QVAPOR', 'no2_vcd', 'hcho_vcd', 'pressure', 'temperature']


sat_aks_feature_sets = ['PHOTR_O31D', 'QVAPOR', 'no2_vcd', 'hcho_vcd_aks', 'pressure', 'temperature']


sat_iso_feature_sets = ['PHOTR_O31D', 'QVAPOR', 'no2_vcd', 'iso_vcd', 'pressure', 'temperature']


sat_surr_aks_feature_sets = ['PHOTR_O31D', 'QVAPOR', 'no2_vcd', 'no2_surr_vcd', 'hcho_vcd_aks', 'pressure',
                             'temperature']

sat_downwind_aks_feature_sets = ['PHOTR_O31D', 'QVAPOR', 'no2_vcd', 'no2_upwind_vcd', 'hcho_vcd_aks', 'pressure',
                             'temperature']

geo_mat = scipy.io.loadmat(os.path.join(Workspace, '../WRFData/Summer_WRF_avg_2005_no2.mat'))
wrf_xlon = geo_mat['xlon']
wrf_xlat = geo_mat['xlat']
f = h5py.File(os.path.join(Workspace, 'site_wrf_data_geos.mat'), 'r')
base_locs = f.get('base_locs')


def collect_info_for_loc(loc_ind):
    ref = base_locs['BoxSize'][0][loc_ind-1]
    radius = np.mean(np.array(base_locs[ref])[2:])
    longitude = np.array(base_locs[base_locs['Longitude'][0][loc_ind-1]])[0][0]
    latitude = np.array(base_locs[base_locs['Latitude'][0][loc_ind-1]])[0][0]
    return longitude, latitude, radius


def find_grids_index_from_data(loc_ind, site_lons, site_lats):
    ref = base_locs['BoxSize'][0][loc_ind-1]
    radius = np.mean(np.array(base_locs[ref])[2:])
    longitude = np.array(base_locs[base_locs['Longitude'][0][loc_ind-1]])[0][0]
    latitude = np.array(base_locs[base_locs['Latitude'][0][loc_ind-1]])[0][0]
    indx_x = []
    indx_y = []
    for site_lon, site_lat in zip(site_lons, site_lats):
        dis = (site_lon - wrf_xlon)**2 + (site_lat - wrf_xlat)**2
        this_indx = np.unravel_index(dis.argmin(), dis.shape)
        indx_x.append(this_indx[0])
        indx_y.append(this_indx[1])
    return indx_x, indx_y, longitude, latitude, radius


def create_wrf_values_2d(loc_ind, site_lons, site_lats, values):
    indx_x, indx_y, longitude, latitude, radius = find_grids_index_from_data(loc_ind, site_lons, site_lats)
    wrf_values = np.empty(wrf_xlon.shape)
    wrf_values[:] = np.nan
    for x,y,value in zip(indx_x, indx_y, values):
        wrf_values[x, y] = value
    return wrf_values, longitude, latitude, radius


def read_wrf_data(city_inds):
    """
    Read the wrf chem data from csv files given the city indexes.

    :param city_inds: the array provide city indexes
    :type city_inds: py:class:`numpy.ndarray`
    :return: a panda dataframe describing wrf chem data.

    """
    data = pd.DataFrame()
    for loc_ind in city_inds:
        filename = 'wrf_data_merged_{:02}.csv'.format(loc_ind)
        this_data = pd.read_csv(os.path.join(wrf_path, filename))
        this_data.loc[:, 'city'] = loc_ind
        data = pd.concat([data, this_data])
    return data


def read_wrf_small_radius_data(city_inds):
    """
    Read the wrf chem data from csv files given the city indexes.

    :param city_inds: the array provide city indexes
    :type city_inds: py:class:`numpy.ndarray`
    :return: a panda dataframe describing wrf chem data.

    """
    data = pd.DataFrame()
    for loc_ind in city_inds:
        filename = 'wrf_data_merged_{:02}.csv'.format(loc_ind)
        print('Reading {}'.format(filename))
        this_data = pd.read_csv(os.path.join(small_radius_wrf_path, filename))
        this_data.loc[:, 'city'] = loc_ind
        data = pd.concat([data, this_data])
    return data


def read_pasadena_wrf_data():
    """
    Read the wrf chem data from csv files given the city indexes.

    :param city_inds: the array provide city indexes
    :type city_inds: py:class:`numpy.ndarray`
    :return: a panda dataframe describing wrf chem data.

    """
    data = pd.DataFrame()
    filename = 'wrf_data_merged_pasadena.csv'
    this_data = pd.read_csv(os.path.join(pasadena_wrf_path, filename))
    this_data.loc[:, 'city'] = 23
    data = pd.concat([data, this_data])
    return data


def read_omi_data(city_inds):
    """
    Read the omi data from csv files given the city indexes, including hcho vcd, no2 vcd and surface o3.

    :param city_inds: the array provide city indexes
    :type city_inds: py:class:`numpy.ndarray`
    :return: a panda dataframe describing omi data.

    """
    data = pd.DataFrame()
    for loc_ind in city_inds:
        filename = 'omi_data_merged_{:02}.csv'.format(loc_ind)
        this_data = pd.read_csv(os.path.join(omi_path, filename))
        this_data.loc[:, 'city'] = loc_ind
        data = pd.concat([data, this_data])
    matlab_date = data.loc[:, 'dvec']
    pdvec = []
    for i in range(data.shape[0]):
        pdvec.append(datetime.datetime.fromordinal(int(matlab_date.iloc[i])) + datetime.timedelta(
            days=int(matlab_date.iloc[i]) % 1) - datetime.timedelta(days=366))
    data.loc[:, 'pdvec'] = np.array(pdvec)
    data.loc[:, 'month'] = pd.DatetimeIndex(data.loc[:, 'pdvec']).month
    data.loc[:, 'year'] = pd.DatetimeIndex(data.loc[:, 'pdvec']).year
    return data


def read_omi_small_radius_data(city_inds):
    """
    Read the omi data from csv files given the city indexes, including hcho vcd, no2 vcd and surface o3.

    :param city_inds: the array provide city indexes
    :type city_inds: py:class:`numpy.ndarray`
    :return: a panda dataframe describing omi data.

    """
    data = pd.DataFrame()
    for loc_ind in city_inds:
        filename = 'omi_data_merged_{:02}.csv'.format(loc_ind)
        this_data = pd.read_csv(os.path.join(small_radius_omi_path, filename))
        this_data.loc[:, 'city'] = loc_ind
        data = pd.concat([data, this_data])
    data['pdvec'] = data['dvec'].apply(
        lambda matlab_datenum: datetime.datetime.fromordinal(int(matlab_datenum)) + datetime.timedelta(
            days=int(matlab_datenum) % 1) - datetime.timedelta(days=366))
    data.loc[:, 'month'] = pd.DatetimeIndex(data.loc[:, 'pdvec']).month
    data.loc[:, 'year'] = pd.DatetimeIndex(data.loc[:, 'pdvec']).year
    return data


def read_obs_data(field_campaign):
    if field_campaign == 'calnex':
        data = pd.read_csv(os.path.join(obs_path, 'merge_calnex.csv'))
    if field_campaign == 'soas':
        data = pd.read_csv(os.path.join(obs_path, 'merge_soas.csv'))
    matlab_date = data.loc[:, 'time']
    pdvec = []
    for i in range(data.shape[0]):
        days = matlab_date.iloc[i] % 1
        hours = days % 1 * 24
        minutes = hours % 1 * 60
        seconds = minutes % 1 * 60
        pdvec.append(datetime.datetime.fromordinal(int(matlab_date.iloc[i])) + datetime.timedelta(
            days=int(days)) + datetime.timedelta(hours=int(hours)) + datetime.timedelta(minutes=int(minutes))
                     + datetime.timedelta(seconds=int(seconds)) - datetime.timedelta(days=366))
    data.loc[:, 'ptime'] = np.array(pdvec)
    data.loc[:, 'pdvec'] = pd.DatetimeIndex(pdvec).date
    data.loc[:, 'hour'] = pd.DatetimeIndex(pdvec).hour
    data.loc[:, 'year'] = pd.DatetimeIndex(pdvec).year
    return data


def wrf_features_process(data, lightning_filter):
    """
    Precoss original wrf inputs, including:
    1) Convert units of chemicals from ppm to molec cm^{-3}
    2) Convert matlab datenum to python datetime, create two separate columns 'year' and 'month'
    3) Calculate Ozone photolysis rate and HCHO phtotolysis rate
    4) Calculate NO2 to NO ratio

    :param data: Panda dataframe containing original wrf data
    :return: Panda dataframes X contains all possible features and y contains true OH values
    """
    #chemicals_lists = ['ho', 'pan', 'no', 'no2', 'acd', 'iso', 'act', 'gly', 'o3', 'hcho']
    chemicals_lists = ['ho']
    chemicals = [chemical for chemical in chemicals_lists if chemical in data.columns]
    for chemical in chemicals:
        data.loc[:, chemical] = data.loc[:, chemical] * data.loc[:, 'ndens'] * 1e-12
    data.loc[:, 'no2_no'] = data.loc[:, 'no2'] / data.loc[:, 'no']
    data.loc[:, 'tau'] = (data.loc[:, 'no'] + data.loc[:, 'no2']) / \
                             (data.loc[:, 'LNOXA']+data.loc[:, 'LNOXHNO3']) / 3600
    if 'hcho_vcd_aks' in data.columns:
        indx = data.loc[:, 'hcho_vcd_aks'] <= 0
        data.loc[indx, 'hcho_vcd_aks'] = np.nan
    else:
        data.loc[:, 'hcho_vcd_aks'] = data.loc[:, 'hcho_vcd']
    data.loc[:, 'flashcounts'] = data.loc[:, 'IC_FLASHCOUNT'] + data.loc[:, 'CG_FLASHCOUNT']
    data['pdvec'] = data['dvec'].apply(
        lambda matlab_datenum: datetime.datetime.fromordinal(int(matlab_datenum)) + datetime.timedelta(
            days=int(matlab_datenum) % 1) - datetime.timedelta(days=366))
    data.loc[:, 'month'] = pd.DatetimeIndex(data.loc[:, 'pdvec']).month
    data.loc[:, 'year'] = pd.DatetimeIndex(data.loc[:, 'pdvec']).year
    nan_rows = data[data.isna().any(axis=1)].shape[0]
    print('There are {} observations containing nan values'.format(nan_rows))
    data = data.dropna()
    if lightning_filter:
        data = wrf_lightning_filter(data)
    y = data.loc[:, ['ho', 'tau']]
    X = data.drop(columns=['ho', 'tau'])
    print('The number of observations is {}'.format(X.shape[0]))
    return X, y


def wrf_features_process_aks(data, lightning_filter):
    """
    Precoss original wrf inputs, including:
    1) Convert units of chemicals from ppm to molec cm^{-3}
    2) Convert matlab datenum to python datetime, create two separate columns 'year' and 'month'
    3) Calculate Ozone photolysis rate and HCHO phtotolysis rate
    4) Calculate NO2 to NO ratio

    :param data: Panda dataframe containing original wrf data
    :return: Panda dataframes X contains all possible features and y contains true OH values
    """
    #chemicals_lists = ['ho', 'pan', 'no', 'no2', 'acd', 'iso', 'act', 'gly', 'o3', 'hcho']
    chemicals_lists = ['ho']
    chemicals = [chemical for chemical in chemicals_lists if chemical in data.columns]
    for chemical in chemicals:
        data.loc[:, chemical] = data.loc[:, chemical] * data.loc[:, 'ndens'] * 1e-12
    data.loc[:, 'no2_no'] = data.loc[:, 'no2'] / data.loc[:, 'no']
    data.loc[:, 'tau'] = (data.loc[:, 'no'] + data.loc[:, 'no2']) / \
                             (data.loc[:, 'LNOXA']+data.loc[:, 'LNOXHNO3']) / 3600

    data.loc[:, 'hcho_vcd_aks'] = data.loc[:, 'qhcho_vcd_aks']
    indx = data.loc[:, 'hcho_vcd_aks'] <= 0
    data.loc[indx, 'hcho_vcd_aks'] = np.nan

    data.loc[:, 'flashcounts'] = data.loc[:, 'IC_FLASHCOUNT'] + data.loc[:, 'CG_FLASHCOUNT']
    data['pdvec'] = data['dvec'].apply(
        lambda matlab_datenum: datetime.datetime.fromordinal(int(matlab_datenum)) + datetime.timedelta(
            days=int(matlab_datenum) % 1) - datetime.timedelta(days=366))
    data.loc[:, 'month'] = pd.DatetimeIndex(data.loc[:, 'pdvec']).month
    data.loc[:, 'year'] = pd.DatetimeIndex(data.loc[:, 'pdvec']).year
    nan_rows = data[data.isna().any(axis=1)].shape[0]
    print('There are {} observations containing nan values'.format(nan_rows))
    data = data.dropna()
    if lightning_filter:
        data = wrf_lightning_filter(data)
    y = data.loc[:, ['ho', 'tau']]
    X = data.drop(columns=['ho', 'tau'])
    print('The number of observations is {}'.format(X.shape[0]))
    return X, y


def wrf_lightning_filter(data):
    data_lightning_filter = data.groupby(['city','pdvec'])['flashcounts'].agg('sum').reset_index()
    indx = data_lightning_filter['flashcounts'] > 0
    data_lightning_filter.loc[indx, 'flashcounts'] = np.nan
    data_lightning_filter = data_lightning_filter.rename(columns={'flashcounts':'total_flashcounts'})
    data = data.merge(data_lightning_filter, on=['city', 'pdvec'])
    nan_rows = data[data.isna().any(axis=1)].shape[0]
    print('There are {} observations containing nan values due to lightning'.format(nan_rows))
    data = data.dropna()
    return data


def wrf_features_process_backup(data):
    """
    Precoss original wrf inputs, including:
    1) Convert units of chemicals from ppm to molec cm^{-3}
    2) Convert matlab datenum to python datetime, create two separate columns 'year' and 'month'
    3) Calculate Ozone photolysis rate and HCHO phtotolysis rate
    4) Calculate NO2 to NO ratio

    :param data: Panda dataframe containing original wrf data
    :return: Panda dataframes X contains all possible features and y contains true OH values
    """
    #chemicals_lists = ['ho', 'pan', 'no', 'no2', 'acd', 'iso', 'act', 'gly', 'o3', 'hcho']
    chemicals_lists = ['ho']
    chemicals = [chemical for chemical in chemicals_lists if chemical in data.columns]
    for chemical in chemicals:
        data.loc[:, chemical] = data.loc[:, chemical] * data.loc[:, 'ndens'] * 1e-12
    data.loc[:, 'no2_no'] = data.loc[:, 'no2'] / data.loc[:, 'no']
    data.loc[:, 'jo1d_h2o'] = data.loc[:, 'QVAPOR'] * data.loc[:, 'PHOTR_O31D']
    data.loc[:, 'o3_jo1d'] = data.loc[:, 'o3'] * data.loc[:, 'PHOTR_O31D']
    data.loc[:, 'o3_photolysis'] = data.loc[:, 'o3'] * data.loc[:, 'PHOTR_O31D'] * data.loc[:, 'QVAPOR']
    data.loc[:, 'hcho_photolysis'] = data.loc[:, 'hcho'] * data.loc[:, 'PHOTR_CH2OR']
    data.loc[:, 'wrf_tau_hno3'] = (data.loc[:, 'no'] + data.loc[:, 'no2'])/data.loc[:, 'LNOXHNO3'] / 3600
    data.loc[:, 'wrf_tau_ans'] = (data.loc[:, 'no'] + data.loc[:, 'no2']) / data.loc[:, 'LNOXA'] / 3600
    data.loc[:, 'tau'] = (data.loc[:, 'no'] + data.loc[:, 'no2']) / \
                             (data.loc[:, 'LNOXA']+data.loc[:, 'LNOXHNO3']) / 3600
    poh, pho2 = cal_poh_pho2(data, data['o3'], data['hcho'])
    data.loc[:, 'poh'] = poh
    data.loc[:, 'pho2'] = pho2
    data.loc[:, 'phox'] = poh + pho2
    data.loc[:, 'pho2_phox'] = pho2 / (pho2 + poh)
    if 'hcho_vcd_aks' in data.columns:
        indx = data.loc[:, 'hcho_vcd_aks'] <= 0
        data.loc[indx, 'hcho_vcd_aks'] = np.nan
    else:
        data.loc[:, 'hcho_vcd_aks'] = data.loc[:, 'hcho_vcd']
    data.loc[:, 'flashcounts'] = data.loc[:, 'IC_FLASHCOUNT'] + data.loc[:, 'CG_FLASHCOUNT']
    data['pdvec'] = data['dvec'].apply(
        lambda matlab_datenum: datetime.datetime.fromordinal(int(matlab_datenum)) + datetime.timedelta(
            days=int(matlab_datenum) % 1) - datetime.timedelta(days=366))
    data.loc[:, 'month'] = pd.DatetimeIndex(data.loc[:, 'pdvec']).month
    data.loc[:, 'year'] = pd.DatetimeIndex(data.loc[:, 'pdvec']).year
    nan_rows = data[data.isna().any(axis=1)].shape[0]
    print('There are {} observations containing nan values'.format(nan_rows))
    data = data.dropna()
    y = data.loc[:, ['ho', 'tau']]
    X = data.drop(columns=['ho', 'tau'])
    print('The number of observations is {}'.format(X.shape[0]))
    return X, y


def omi_features_process(data, wrf_data):
    for sat_feature, rename_sat_feature in zip(satellite_features, wrf_rename_satellite_features):
        wrf_data = wrf_data.rename(columns={sat_feature : rename_sat_feature})
    print('Merging WRF data after renaming satellite features {}'.format(satellite_features))
    data = data.merge(wrf_data, on=['pdvec', 'city'], how='left')
    data = data.rename(columns={"behr_no2": "no2_vcd", "qa4ecv_hcho": "hcho_vcd", "omi_o3": "o3"})
    data.loc[:, 'no2_vcd_aks'] = data.loc[:, 'no2_vcd']
    data.loc[:, 'hcho_vcd_aks'] = data.loc[:, 'hcho_vcd']
    for column in satellite_features:
        indx = data.loc[:, column] < 0
        data.loc[indx, column] = np.nan
    indx = data.loc[:, 'hcho_vcd_aks'] > 1e17
    data.loc[indx, column] = np.nan

    used_data = data.loc[:, sat_aks_feature_sets]
    nan_rows = used_data[used_data.isna().any(axis=1)].shape[0]
    print('There are {} observations containing nan values'.format(nan_rows))
    data = data.dropna(subset=wrf_rename_satellite_features + sat_aks_feature_sets)
    print('The number of observations is {}'.format(data.shape[0]))
    data.loc[:, 'wrf_oh'] = data.loc[:, 'oh']*1e6
    return data


def omi_small_radius_features_process(data, wrf_data):
    for sat_feature, rename_sat_feature in zip(satellite_features, wrf_rename_satellite_features):
        wrf_data = wrf_data.rename(columns={sat_feature:rename_sat_feature})
    print('Merging WRF data after renaming satellite features {}'.format(satellite_features))
    data = data.merge(wrf_data, on=['pdvec', 'city', 'xlon', 'xlat'], how='left')
    data = data.rename(columns={"behr_no2": "no2_vcd", "qa4ecv_hcho": "hcho_vcd"})
    data.loc[:, 'no2_vcd_aks'] = data.loc[:, 'no2_vcd']
    data.loc[:, 'hcho_vcd_aks'] = data.loc[:, 'hcho_vcd']
    for column in satellite_small_radius_features:
        indx = data.loc[:, column] < 0
        data.loc[indx, column] = np.nan
    indx = data.loc[:, 'hcho_vcd_aks'] > 1e17
    data.loc[indx, column] = np.nan
    used_data = data.loc[:, sat_aks_feature_sets]
    nan_rows = used_data[used_data.isna().any(axis=1)].shape[0]
    print('There are {} observations containing nan values'.format(nan_rows))
    data = data.dropna(subset=sat_aks_feature_sets)
    print('The number of observations is {}'.format(data.shape[0]))
    return data


def omi_small_radius_nasa_features_process(data, wrf_data):
    for sat_feature, rename_sat_feature in zip(satellite_features, wrf_rename_satellite_features):
        wrf_data = wrf_data.rename(columns={sat_feature:rename_sat_feature})
    print('Merging WRF data after renaming satellite features {}'.format(satellite_features))
    data = data.merge(wrf_data, on=['pdvec', 'city', 'xlon', 'xlat'], how='left')
    data = data.rename(columns={"behr_no2": "no2_vcd", "nasa_hcho": "hcho_vcd"})
    data.loc[:, 'no2_vcd_aks'] = data.loc[:, 'no2_vcd']
    data.loc[:, 'hcho_vcd_aks'] = data.loc[:, 'hcho_vcd']
    for column in satellite_small_radius_features:
        indx = data.loc[:, column] < 0
        data.loc[indx, column] = np.nan
    indx = data.loc[:, 'hcho_vcd_aks'] > 1e17
    data.loc[indx, column] = np.nan
    used_data = data.loc[:, sat_aks_feature_sets]
    nan_rows = used_data[used_data.isna().any(axis=1)].shape[0]
    print('There are {} observations containing nan values'.format(nan_rows))
    data = data.dropna(subset=sat_aks_feature_sets)
    print('The number of observations is {}'.format(data.shape[0]))
    return data


def omi_small_radius_qa4ecv_l3_features_process(data, wrf_data):
    for sat_feature, rename_sat_feature in zip(satellite_features, wrf_rename_satellite_features):
        wrf_data = wrf_data.rename(columns={sat_feature:rename_sat_feature})
    print('Merging WRF data after renaming satellite features {}'.format(satellite_features))
    data = data.merge(wrf_data, on=['pdvec', 'city', 'xlon', 'xlat'], how='left')
    data = data.rename(columns={"behr_no2": "no2_vcd", "qa4ecv_l3_hcho": "hcho_vcd"})
    data.loc[:, 'no2_vcd_aks'] = data.loc[:, 'no2_vcd']
    data.loc[:, 'hcho_vcd_aks'] = data.loc[:, 'hcho_vcd']
    for column in satellite_small_radius_features:
        indx = data.loc[:, column] < 0
        data.loc[indx, column] = np.nan
    indx = data.loc[:, 'hcho_vcd_aks'] > 1e17
    data.loc[indx, column] = np.nan
    used_data = data.loc[:, sat_aks_feature_sets]
    nan_rows = used_data[used_data.isna().any(axis=1)].shape[0]
    print('There are {} observations containing nan values'.format(nan_rows))
    data = data.dropna(subset=sat_aks_feature_sets)
    print('The number of observations is {}'.format(data.shape[0]))
    return data


def obs_features_process(data, wrf_data, field_campaign):
    if field_campaign == 'calnex':
        wrf_data = wrf_data.query('city == 23')
        wrf_data['pdvec'] = pd.DatetimeIndex(wrf_data['pdvec']).date
        data = data.merge(wrf_data, on=['pdvec'], how='left')
    return data


def cal_poh_pho2(data, o3, hcho):
    """
    :param data:
    :param o3: ppm
    :param hcho: ppm
    :return:
    """
    ndens = data.loc[:, 'ndens']
    qvapor = data.loc[:, 'QVAPOR']
    jo1d = data.loc[:, 'PHOTR_O31D']
    jhcho = data.loc[:, 'PHOTR_CH2OR']
    temp = data.loc[:, 'temperature']
    o3 = o3 * ndens * 1e-6
    hcho = hcho * ndens * 1e-6
    h2o_rate = 2.2e-10
    density = ndens / 6.022e23 * 28.97 / 1000
    h2o = qvapor * density / 18.02e-3 * 6.022e23
    j_o1d = jo1d / 60
    j_hcho = jhcho / 60
    k_o1d = 0.78084 * ndens * 1.8e-11 * np.exp(107 / temp) + 0.20946 * ndens * 3.2e-11 * np.exp(67 / temp)
    o1d = (j_o1d * o3) / (h2o_rate * h2o + k_o1d)
    phox_hcho = 2 * hcho * j_hcho
    phox_o1d = 2 * o1d * h2o * h2o_rate
    return phox_o1d, phox_hcho


def cal_k_rates(temp, ndens):
    def cal_k_oh_no2(temp, ndens):
        k0_300K = 1.49e-30
        n = 1.8
        kinf_300k = 2.58e-11
        m = 0
        zt_help = 300.0 / temp
        k0_T = k0_300K * zt_help**(n) * ndens
        kinf_T = kinf_300k * zt_help**(m)
        k_ratio = k0_T / kinf_T
        k_ohno2 = k0_T / (1.0 + k_ratio) * 0.6 ** (1.0 / (1.0 + np.log10(k_ratio)**2))
        return k_ohno2
    krates = dict()
    krates['k_ho2no'] = 3.5e-12 * np.exp(250/temp)
    krates['k_ohh2co'] = 5.5e-12 * np.exp(125/temp)
    krates['k_ohno2'] = cal_k_oh_no2(temp, ndens)
    krates['k_ro2ho2'] = 8e-12
    krates['k_ro2no'] = 8e-12
    krates['k_ro2ro2'] = 6.8e-14
    krates['k_ho2ho2'] = 6.8e-14
    return krates


def prepare_box_model_inputs_for_single_city(loc_ind, X_omi, o3_var, hcho_var):
    this_omi = X_omi[X_omi['city'] == loc_ind]
    [this_omi.loc[:, 'infer_no2'], this_omi.loc[:, 'no2_huber_res'], this_omi.loc[:, 'no2_linear_res']] = \
        infer_surface_mixingratio_from_vcd(this_omi, 'no2')
    [this_omi.loc[:, 'infer_hcho'], this_omi.loc[:, 'hcho_huber_res'], this_omi.loc[:, 'hcho_linear_res']] = \
        infer_surface_mixingratio_from_vcd(this_omi, 'hcho')
    this_omi.loc[:, 'wrf_hcho'] = this_omi.loc[:, 'hcho'] * this_omi['ndens'] * 1e-6
    this_omi = cal_poh_pho2(this_omi, this_omi.loc[:, o3_var], this_omi.loc[:, hcho_var])
    this_omi.loc[:, 'wrf_tau_hno3'] = (this_omi.loc[:, 'no'] + this_omi.loc[:, 'no2']) / \
                                 (this_omi.loc[:, 'LNOXHNO3']) / 3600
    this_omi.loc[:, 'wrf_tau_ans'] = (this_omi.loc[:, 'no'] + this_omi.loc[:, 'no2']) / \
                                      (this_omi.loc[:, 'LNOXA']) / 3600
    this_omi.loc[:, 'wrf_tau'] = (this_omi.loc[:, 'no'] + this_omi.loc[:, 'no2']) / \
                                 (this_omi.loc[:, 'LNOXHNO3'] + this_omi.loc[:, 'LNOXA']) / 3600
    this_omi.loc[:, 'o3'] = this_omi.loc[:, 'o3'] * this_omi['ndens'] * 1e-6
    this_omi.loc[:, 'infer_nox'] = this_omi.loc[:, 'infer_no2'] * (1 + this_omi.loc[:, 'no2_no']) / this_omi.loc[:,
                                                                                        'no2_no']
    this_omi.loc[:, 'wrf_nox'] = this_omi.loc[:, 'no2'] * this_omi['ndens'] *1e-6 * \
                                 (1 + this_omi.loc[:, 'no2_no']) / this_omi.loc[:, 'no2_no']

    this_omi.loc[:, 'phox'] = this_omi.loc[:, 'poh'] + this_omi.loc[:, 'pho2']
    this_omi.loc[:, 'pho2_phox'] = this_omi.loc[:, 'pho2'] / this_omi.loc[:, 'phox']
    this_omi_annual = this_omi.groupby(['year_y']).agg('mean')
    this_omi_annual = this_omi_annual.loc[:,
                      ['phox', 'pho2_phox', 'o3', 'infer_hcho', 'infer_no2', 'infer_nox', 'ndens', 'temperature',
                       'pressure', 'no2_no', 'oh_pred']]
    this_omi_annual.loc[:, 'city'] = loc_ind
    return this_omi, this_omi_annual.reset_index()


def steady_state_cal(data, oh_varname, nox_varname, hcho_varname, alphaeff):
    krates = cal_k_rates(data['temperature'], data['ndens'])
    poh = data['phox'] * (1 - data['pho2_phox'])
    pho2 = data['phox'] * data['pho2_phox']
    no = data[nox_varname] / (1 + data['no2_no'])
    no2 = data[nox_varname] - no
    tau_hno3 = data[nox_varname] / (krates['k_ohno2'] * data[oh_varname] * no2 * 3600)
    phox = data['phox']
    pho2_phox = data['pho2_phox']
    oh = data[oh_varname]
    hcho = data[hcho_varname]

    phox_total = (phox * pho2_phox / 2 + krates['k_ohh2co'] * oh * hcho) / alphaeff
    ho2 = (phox_total - poh) / (krates['k_ho2no'] * no)
    rate_1st = 2 * krates['k_ro2ro2']
    rate_2nd = krates['k_ro2no'] * no + krates['k_ro2ho2'] * ho2
    ro2 = (-rate_2nd + np.sqrt(rate_2nd**2 + 4*rate_1st*phox_total))/(2*rate_1st)
    alpha = 1 - (ho2 * krates['k_ho2no'] * no + 2 * krates['k_ho2ho2'] * ho2 ** 2 +
                 krates['k_ro2ho2'] * ro2 * ho2 - pho2) /(krates['k_ro2no'] * ro2 * no)
    vocr = phox_total / data[oh_varname]
    tau_ans = data[nox_varname] / (alpha * krates['k_ro2no'] * ro2 * no * 3600)
    tau_ans_fix = data[nox_varname] / (0.04 * krates['k_ro2no'] * ro2 * no * 3600)
    tau = 1 / (1 / tau_hno3 + 1 / tau_ans)
    tau_fix = 1 / (1 / tau_hno3 + 1 / tau_ans_fix)
    data.loc[:, 'vocr'] = vocr
    data.loc[:, 'alpha'] = alpha
    data.loc[:, 'tau_hno3_ss'] = tau_hno3
    data.loc[:, 'tau_ans_ss'] = tau_ans
    data.loc[:, 'tau_ans_fix_ss'] = tau_ans_fix
    data.loc[:, 'tau_ss'] = tau
    data.loc[:, 'tau_fix_ss'] = tau_fix
    data.loc[:, 'ho2'] = ho2
    data.loc[:, 'ro2'] = ro2
    #indx = (data.loc[:, 'alpha'] < 0) | (data.loc[:, 'alpha'] > 0.3)
    #data.loc[indx, 'alpha'] = np.nan
    #data.loc[indx, 'tau'] = np.nan
    #data.loc[indx, 'tau_ans'] = np.nan
    return data


def solve_steady_state_equations(temperature, phox, poh, pho2, pho2_phox, oh, hcho, no, no2):
    alphaeffs = np.arange(0.1, 0.4, 0.02)
    alphas = []
    for alphaeff in alphaeffs:
        result = solve_steady_state_equations_fixedalphaeff(temperature, phox, poh, pho2,
                                                            pho2_phox, oh, hcho, no, no2, alphaeff)
        alphas.append(result[0]/100)


def solve_steady_state_equations_fixedalphaeff(temperature, phox, poh, pho2, pho2_phox, oh, hcho, no, no2, alphaeff):
    krates = cal_k_rates(temperature)
    alpha, phox_total, ho2, ro2 = sym.symbols('alpha, phox_total, ho2, ro2')
    eq1 = sym.Eq(alphaeff * phox_total, (phox * pho2_phox / 2 +
                                   krates['k_ohh2co'] * oh * hcho)/1e8)
    eq2 = sym.Eq(phox_total - ho2 * krates['k_ho2no'] * no, poh/1e8)
    eq3 = sym.Eq(ho2 * krates['k_ho2no'] * no + 2 * krates['k_ho2ho2'] * ho2 ** 2 * 1e8 +
                 krates['k_ro2ho2'] * ro2 * ho2 * 1e8,
           krates['k_ro2no'] * ro2 * no * (1 - alpha/100) + pho2 / 1e8)
    eq4 = sym.Eq(phox_total, ro2 * (krates['k_ro2no'] * no +
                                                krates['k_ro2ho2'] * ho2 * 1e8 +
                                                2 * krates['k_ro2ro2'] * ro2 * 1e8))
    eq5 = sym.Eq(phox/1e8, alpha / 100 * krates['k_ro2no'] * ro2 * no + krates['k_ohno2'] * oh * no2 / 1e8 +
           2 * krates['k_ro2ho2'] * ro2 * ho2 * 1e8 + 2 * krates['k_ro2ro2'] * ro2 ** 2 * 1e8 +
           2 * krates['k_ho2ho2'] * ho2 ** 2 * 1e8)
    result = sym.solve([eq1, eq2, eq3, eq4], (alpha, phox_total, ho2, ro2))

    return result[0]


def infer_surface_mixingratio_from_vcd(data, species):
    if species == 'no2':
        suf_var = 'no2'
        vcd_var = 'wrf_no2_vcd'
        sat_var = 'no2_vcd'
    elif species == 'hcho':
        suf_var = 'hcho'
        vcd_var = 'wrf_hcho_vcd_aks'
        sat_var = 'hcho_vcd_aks'
    n_obs = data.shape[0]
    x = np.column_stack([data[suf_var].values * data['ndens'] / 1e16, np.ones(n_obs)])
    y = data[vcd_var] / 1e15
    linear = LinearRegression().fit(x, y)
    linear_res = ('y = {:.2f}x + {:.2f}, r2:{:.2f}, n:{:d}'.format(linear.coef_[0], linear.intercept_,
                                                                   linear.score(x, y), n_obs))
    huber = HuberRegressor(epsilon=2).fit(x, y)
    huber_res = ('y = {:.2f}x + {:.2f}, r2:{:.2f}, n:{:d}'.format(huber.coef_[0], huber.intercept_,
                                                                  huber.score(x[~huber.outliers_, :],
                                                                              y[~huber.outliers_]),
                                                                  n_obs - np.sum(huber.outliers_)))
    infer_conc = (data[sat_var]/1e15 - huber.intercept_)*1e16/(huber.coef_[0]*data['ndens'])
    return infer_conc, huber_res, linear_res


def find_features_merge_from_wrf(feature_opt):
    if feature_opt == 'default':
        feature_sets = opt_feature_sets
    elif feature_opt == 'isoprene':
        feature_sets = iso_feature_sets
    elif feature_opt == 'oh_index':
        feature_sets = oh_index_feature_sets
    else:
        feature_sets = full_feature_sets

    merge_features = [feature for feature in feature_sets if feature not in satellite_features]
    return merge_features


def read_site_location_ncdf():
    data = nc.Dataset(os.path.join(Workspace, 'trend_locations.nc'))
    lon = data['Longitude'][:]
    lat = data['Latitude'][:]
    ShortName = data['ShortName']
    names = [nc.chartostring(ShortName[:, i]) for i in range(ShortName.shape[1])]
    Radius = data['Radius'][:]
    city_index = np.arange(1, lon.shape[0]+1)
    return names, lon, lat, Radius, city_index


def read_site_location():
    data = pd.read_csv(os.path.join(output_path, 'City_locations.csv'))
    data.loc[:, 'city_index'] = data.index+1
    indx = data.loc[:, 'city_index'] < 50
    shortname = data.loc[indx, 'ShortName'].values
    lon = data.loc[indx, 'Longitude'].values
    lat = data.loc[indx, 'Latitude'].values
    return shortname, lon, lat


def find_save_path(feature_opt):
    if feature_opt == 'default':
        save_path = def_output_path
    elif feature_opt == 'isoprene':
        save_path = iso_output_path
    elif feature_opt == 'oh_index':
        save_path = oh_index_output_path
    elif feature_opt == 'satellite_aks':
        save_path = sat_aks_output_path
    else:
        save_path = output_path
    return save_path


if __name__ == "__main__":
    read_obs_data('soas')
    city_index = [1]
    feature_opt = 'satellite_aks'
    wrf_data = read_wrf_data(city_index)
    omi_data = read_omi_data(city_index)
    X_omi, y_omi = omi_features_process(omi_data, wrf_data)

    feature_opt = 'oh_index'
    wrf_data = read_wrf_data([23])
    X_wrf, y_wrf = wrf_features_process(wrf_data)
    obs_data = read_obs_data('calnex')
    X_obs = obs_features_process(obs_data, wrf_data, 'calnex')
    read_obs_data('calnex')
    shortname, lon, lat = read_site_location()