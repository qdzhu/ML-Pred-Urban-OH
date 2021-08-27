import RandomForest
import numpy as np
from Utils import *
from sklearn.metrics import mean_squared_error
from sklearn.metrics import r2_score
from sklearn.cluster import KMeans
from scipy.spatial.distance import cdist
import matplotlib.pyplot as plt
import os
import seaborn as sns
from RunRF import iso_output_path
from RunRF import def_output_path
from RunRF import output_path
import datashader as ds
import matplotlib.pyplot as plt
import csv

class PlotProcedures:
    @staticmethod
    def plot_cluster_feature_importance(save_path, feature_opt):
        if feature_opt == 'default':
            features = RandomForest.opt_feature_sets
        elif feature_opt == 'isoprene':
            features = RandomForest.iso_feature_sets
        else:
            print('Have implemented for this feature set')
            exit()
        importance_cluster = np.load(os.path.join(save_path, 'feature_importance_cluster.npy'))
        plt.figure(figsize=(10, 5))
        for i, importance in enumerate(importance_cluster):
            plt.subplot(1, 4, i + 1)
            indices = np.argsort(-1 * importance)
            sns.barplot(importance[indices], features)
            plt.yticks(range(len(indices)), [features[i] for i in indices])
            plt.xlabel('Relative Importance')
        plt.rcParams.update({'font.size': 50})
        plt.tight_layout()
        plt.savefig(os.path.join(save_path, 'cluster_feature_importance.png'))

    @staticmethod
    def plot_city_maps():
        city_index = np.arange(1, 50)
        ShortName, lon, lat, Radius, city_index_full = read_site_location()

        lon_lim = [-125, 65]
        lat_lim = [25, 50]

    @staticmethod
    def plot_model_sensitivity_to_noise(feature_opt):
        uncerts = np.linspace(0, 1.05, 22)
        save_path = find_save_path(feature_opt)
        param_names = ['no2_vcd', 'hcho_vcd', 'o3']
        rmses = dict()
        for param_name in param_names:
            save_name = 'rmse_' + param_name + '_noise.npy'
            this_rmse = np.load(os.path.join(save_path, save_name))
            rmses[param_name] = this_rmse
        plt.figure(figsize=(10,5))
        for key, value in rmses.items():
            sns.lineplot(uncerts, value, label=key)
        plt.show()

    @staticmethod
    def plot_iso_jo1d_h2o_variaions_cities():
        feature_opt = 'isoprene'
        save_path = find_save_path(feature_opt)
        cluster_name = os.path.join(save_path, 'city_clusters.npy')
        city_cluster = np.load(cluster_name)
        #ShortName, lon, lat, Radius, city_index_full = read_site_location()
        #city_cluster_names = [name for (i, name) in enumerate(ShortName) if city_index_full[i] in city_locs]
        city_index = np.arange(1, 50)
        data = read_wrf_data(city_index)
        X, y = wrf_features_process(data)
        X.loc[:, 'ho'] = y.values
        X = X.set_index('city')
        X.loc[:, 'PHOTR_O31D'] = X.loc[:, 'PHOTR_O31D']/60*1e5
        features = ['PHOTR_O31D', 'iso', 'QVAPOR']
        for i_cluster, city_c in enumerate(city_cluster):
            city = [city for city in city_c if city > 0]
            X.loc[city, 'cluster'] = i_cluster
        plt.figure(figsize=(5, 5))
        sns.scatterplot(x='PHOTR_O31D', y='ho', data=X.reset_index(), s=1)
        sns.kdeplot(x='PHOTR_O31D', y='ho', data=X.reset_index(), fill=True, alpha=0.5)
        plt.savefig(os.path.join(save_path, 'oh_jo1d_wrf_urban.png'))

        plt.figure(figsize=(10, 8))
        plt.subplot(1, 3, 1)
        sns.scatterplot(x='iso', y='ho', hue='city', data=X.reset_index().query('cluster == 0'))
        plt.subplot(1, 3, 2)
        sns.scatterplot(x='QVAPOR', y='ho', hue='PHOTR_O31D', data=X.query('cluster == 1'))
        plt.subplot(1, 3, 3)
        sns.scatterplot(x='PHOTR_O31D', y='ho', hue='QVAPOR', data=X.query('cluster == 2'))
        plt.show()

        x_city = X.groupby(['city'])[features].agg(['mean', 'std'])
        for i_cluster, city_c in enumerate(city_cluster):
            city = [city for city in city_c if city > 0]
            x_city.loc[city, 'cluster'] = i_cluster
        plt.figure(figsize=(10, 8))
        plt.subplot(1, 3, 1)
        sns.scatterplot(x_city['iso']['mean'], x_city['iso']['std'], hue=x_city['cluster'])
        plt.subplot(1, 3, 2)
        sns.scatterplot(x_city['PHOTR_O31D']['mean'], x_city['PHOTR_O31D']['std'], hue=x_city['cluster'])
        plt.subplot(1, 3, 3)
        sns.scatterplot(x_city['QVAPOR']['mean'], x_city['QVAPOR']['std'], hue=x_city['cluster'])
        plt.tight_layout()
        plt.show()




if __name__ == "__main__":
    feature_opt = 'isoprene'
    if feature_opt == 'default':
        save_path = def_output_path
    elif feature_opt == 'isoprene':
        save_path = iso_output_path
    else:
        save_path = output_path

    PlotProcedures.plot_cluster_feature_importance(save_path, feature_opt)
    PlotProcedures.plot_iso_jo1d_h2o_variaions_cities()
    PlotProcedures.plot_model_sensitivity_to_noise(feature_opt)
    PlotProcedures.plot_city_maps()
