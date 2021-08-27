from sklearn.metrics import mean_squared_error
from sklearn.metrics import r2_score
import XGBoostModel
import numpy as np
import xgboost as xgb
from Utils import *
from sklearn.cluster import KMeans
from RandomForest import train_test_split
from RandomForest import feature_selection
from scipy.spatial.distance import cdist
import matplotlib.pyplot as plt
import random

class RunXGB:
    @staticmethod
    def make_forward_feature_selection(feature_opts, save_path, lightning_filter):
        city_index = np.arange(1, 50)
        data = read_wrf_small_radius_data(city_index)
        X, y = wrf_features_process(data, lightning_filter)
        [x_train_unfilter, x_test_unfilter, y_train, y_test] = train_test_split(
            X, y, test_size=0.25, random_state=50)
        x_train, x_test, features = feature_selection(x_train_unfilter, x_test_unfilter, feature_opts)
        r2_col = []
        rmse_col = []
        for i in range(len(x_train.columns)):
            this_x_train = x_train.iloc[:, :i + 1]
            this_x_test = x_test.iloc[:, :i + 1]
            model = XGBoostModel.train_xgb_ml_model(this_x_train, y_train['ho'])
            pred = XGBoostModel.make_predictions(model, this_x_test)
            rmse_col.append(mean_squared_error(y_test['ho'], pred))
            r2_col.append(r2_score(y_test['ho'], pred))
        eval = np.concatenate((rmse_col, r2_col), axis=0)
        np.save(os.path.join(save_path, 'forward_model_eval'), eval)

    def make_forward_feature_selection_shuffle(lightning_filter):
        save_path = './Outputs/SAT_AKS_SR_Shuffle'
        city_index = np.arange(1, 50)
        data = read_wrf_small_radius_data(city_index)
        X, y = wrf_features_process(data, lightning_filter)
        [x_train_unfilter, x_test_unfilter, y_train, y_test] = train_test_split(
            X, y, test_size=0.25, random_state=50)
        for k in range(10):
            x_train, x_test, features = feature_selection(x_train_unfilter, x_test_unfilter, 'shuffle')
            r2_col = []
            rmse_col = []
            for i in range(len(x_train.columns)):
                this_x_train = x_train.iloc[:, :i + 1]
                this_x_test = x_test.iloc[:, :i + 1]
                model = XGBoostModel.train_xgb_ml_model(this_x_train, y_train['ho'])
                pred = XGBoostModel.make_predictions(model, this_x_test)
                rmse_col.append(mean_squared_error(y_test['ho'], pred))
                r2_col.append(r2_score(y_test['ho'], pred))
            eval = np.concatenate((features, rmse_col, r2_col), axis=0)
            np.save(os.path.join(save_path, 'forward_model_eval_{}'.format(k)), eval)

    def make_forward_feature_selection_improve(save_path, lightning_filter):
        city_index = np.arange(1, 50)
        data = read_wrf_small_radius_data(city_index)
        X, y = wrf_features_process(data, lightning_filter)
        [x_train_unfilter, x_test_unfilter, y_train, y_test] = train_test_split(
            X, y, test_size=0.25, random_state=50)
        original_feature_sets = sat_aks_feature_sets
        added_feature_sets = ['iso', 'iso_vcd', 'pan', 'o3', 'PHOTR_CH2OR', 'PBLH', 'no', 'PHOTR_NO2', 'no2', 'hcho', 'PO3']
        r2_col = []
        rmse_col = []
        feature_sets = original_feature_sets.copy()
        print(feature_sets)
        this_x_train = x_train_unfilter.loc[:, feature_sets]
        this_x_test = x_test_unfilter.loc[:, feature_sets]
        model = XGBoostModel.train_xgb_ml_model(this_x_train, y_train['ho'])
        pred = XGBoostModel.make_predictions(model, this_x_test)
        rmse_col.append(mean_squared_error(y_test['ho'], pred))
        r2_col.append(r2_score(y_test['ho'], pred))
        feature_sets.clear()
        print(rmse_col)
        print(r2_col)
        for feature in added_feature_sets:
            feature_sets = original_feature_sets.copy()
            feature_sets.append(feature)
            print(feature_sets)
            this_x_train = x_train_unfilter.loc[:, feature_sets]
            this_x_test = x_test_unfilter.loc[:, feature_sets]
            model = XGBoostModel.train_xgb_ml_model(this_x_train, y_train['ho'])
            pred = XGBoostModel.make_predictions(model, this_x_test)
            rmse_col.append(mean_squared_error(y_test['ho'], pred))
            r2_col.append(r2_score(y_test['ho'], pred))
            feature_sets.clear()
        eval = np.concatenate((rmse_col, r2_col), axis=0)
        np.save(os.path.join(save_path, 'forward_model_eval'), eval)

    @staticmethod
    def make_cross_validation(feature_opts, save_path, lightning_filter, model_type):
        city_index = np.arange(1, 50)
        data = read_wrf_small_radius_data(city_index)
        X, y = wrf_features_process(data)
        print('number of obs:{}'.format(X.shape[0]))
        if lightning_filter:
            indx = X.loc[:, 'flashcounts'] == 0
            X = X.loc[indx, :]
            y = y.loc[indx]
            print('number of obs after lightning filtering:{}'.format(X.shape[0]))
        [x_train_unfilter, x_test_unfilter, y_train, y_test] = train_test_split(
            X, y, test_size=0.25, random_state=50)
        x_train, x_test, features = feature_selection(x_train_unfilter, x_test_unfilter, feature_opts)
        best_params = XGBoostModel.param_tuning_xgb(x_train, y_train.loc[:, model_type])

    @staticmethod
    def oh_index_general_model():
        city_index = np.arange(1, 50)
        data = read_wrf_small_radius_data(city_index)
        X, y = wrf_features_process(data)
        feature_opt = 'oh_index'
        model, train, test, features = XGBoostModel.make_model(X, y, feature_opt)
        importance = model.get_score(importance_type='gain')
        rmse = mean_squared_error(test.loc[:, 'y_true'], test.loc[:, 'y_pred'])
        r2 = r2_score(test.loc[:, 'y_true'], test.loc[:, 'y_pred'])
        importance = np.append(importance, rmse)
        importance = np.append(importance, r2)
        np.save(os.path.join(oh_index_sr_output_path, 'general_model_eval'), importance)

    @staticmethod
    def sat_aks_general_model(lightning_filter):
        city_index = np.arange(1, 50)
        data = read_wrf_small_radius_data(city_index)
        X_wrf, y_wrf = wrf_features_process(data, lightning_filter)
        feature_opt = 'satellite_aks'
        model, train, test, features = XGBoostModel.make_model(X_wrf, y_wrf['ho'], feature_opt)
        model_name = os.path.join(sat_aks_sr_output_path, 'oh_model_filer_validation.model')
        model.save_model(model_name)
        train.to_pickle(os.path.join(sat_aks_sr_output_path, 'ml_wrf_train.pkl'))
        test.to_pickle(os.path.join(sat_aks_sr_output_path, 'ml_wrf_test.pkl'))

    @staticmethod
    def sat_aks_general_model_aks_update(lightning_filter):
        city_index = np.arange(1, 50)
        data = read_wrf_small_radius_data(city_index)
        X_wrf, y_wrf = wrf_features_process_aks(data, lightning_filter)
        feature_opt = 'satellite_aks'
        model, train, test, features = XGBoostModel.make_model(X_wrf, y_wrf['ho'], feature_opt)
        model_name = os.path.join(sat_aks_sr_output_path, 'oh_model_filer_validation_aks_update.model')
        model.save_model(model_name)
        train.to_pickle(os.path.join(sat_aks_sr_output_path, 'ml_wrf_train_aks_update.pkl'))
        test.to_pickle(os.path.join(sat_aks_sr_output_path, 'ml_wrf_test_aks_update.pkl'))

    @staticmethod
    def sat_aks_noisy_general_model(lightning_filter):
        city_index = np.arange(1, 50)
        data = read_wrf_small_radius_data(city_index)
        X_wrf, y_wrf = wrf_features_process(data, lightning_filter)
        feature_opt = 'satellite_aks'
        param_noise = dict({'no2_vcd': 0.3, 'hcho_vcd_aks': 0.6})
        train, test, noise_model= XGBoostModel.make_noise_model(X_wrf, y_wrf['ho'], feature_opt, param_noise)
        model_name = os.path.join(sat_aks_sr_output_path, 'oh_model_filer_validation_noise.model')
        noise_model.save_model(model_name)
        train.to_pickle(os.path.join(sat_aks_sr_output_path, 'ml_wrf_train_noise.pkl'))
        test.to_pickle(os.path.join(sat_aks_sr_output_path, 'ml_wrf_test_noise.pkl'))

    @staticmethod
    def make_final_datasets(lightning_filter):
        city_index = np.arange(1, 50)
        data = read_wrf_small_radius_data(city_index)
        X_wrf, y_wrf = wrf_features_process(data, lightning_filter)
        feature_opt = 'satellite_aks'
        model, train = XGBoostModel.make_final_model(X_wrf, y_wrf['ho'], feature_opt)
        model_name = os.path.join(sat_aks_sr_output_path, 'oh_model_filer.model')
        model.save_model(model_name)
        X_wrf.loc[:, 'wrf_pred'] = train.loc[:, 'y_pred']
        X_wrf.loc[:, 'wrf_oh'] = y_wrf['ho']
        X_wrf.to_pickle(os.path.join(sat_aks_sr_output_path, 'ml_wrf_filer.pkl'))
        omi_data = read_omi_small_radius_data(city_index)
        X_omi = omi_small_radius_features_process(omi_data, X_wrf)
        x, x, features = feature_selection(X_omi, X_omi, feature_opt)
        y_pred = XGBoostModel.make_predictions(model, x)
        X_omi.loc[:, 'omi_pred'] = y_pred
        X_omi.to_pickle(os.path.join(sat_aks_sr_output_path, 'ml_wrf_omi_filer.pkl'))

    @staticmethod
    def make_final_datasets_aks_update(lightning_filter):
        city_index = np.arange(1, 50)
        data = read_wrf_small_radius_data(city_index)
        X_wrf, y_wrf = wrf_features_process_aks(data, lightning_filter)
        feature_opt = 'satellite_aks'
        model, train = XGBoostModel.make_final_model(X_wrf, y_wrf['ho'], feature_opt)
        model_name = os.path.join(sat_aks_sr_output_path, 'oh_model_filer_aks_update.model')
        model.save_model(model_name)
        X_wrf.loc[:, 'wrf_pred'] = train.loc[:, 'y_pred']
        X_wrf.loc[:, 'wrf_oh'] = y_wrf['ho']
        X_wrf.to_pickle(os.path.join(sat_aks_sr_output_path, 'ml_wrf_filer_aks_update.pkl'))
        omi_data = read_omi_small_radius_data(city_index)
        X_omi = omi_small_radius_qa4ecv_l3_features_process(omi_data, X_wrf)
        x, x, features = feature_selection(X_omi, X_omi, feature_opt)
        y_pred = XGBoostModel.make_predictions(model, x)
        X_omi.loc[:, 'omi_pred'] = y_pred
        X_omi.to_pickle(os.path.join(sat_aks_sr_output_path, 'ml_wrf_omi_filer_aks_update.pkl'))

    @staticmethod
    def make_final_datasets_nasahcho(lightning_filter):
        city_index = np.arange(1, 50)
        data = read_wrf_small_radius_data(city_index)
        X_wrf, y_wrf = wrf_features_process(data, lightning_filter)
        feature_opt = 'satellite_aks'
        model, train = XGBoostModel.make_final_model(X_wrf, y_wrf['ho'], feature_opt)
        X_wrf.loc[:, 'wrf_pred'] = train.loc[:, 'y_pred']
        X_wrf.loc[:, 'wrf_oh'] = y_wrf['ho']
        omi_data = read_omi_small_radius_data(city_index)
        X_omi = omi_small_radius_nasa_features_process(omi_data, X_wrf)
        x, x, features = feature_selection(X_omi, X_omi, feature_opt)
        y_pred = XGBoostModel.make_predictions(model, x)
        X_omi.loc[:, 'omi_pred'] = y_pred
        X_omi.to_pickle(os.path.join(sat_aks_sr_output_path, 'ml_wrf_omi_filer_nasa.pkl'))

    @staticmethod
    def make_final_datasets_qa4ecv_l3_hcho(lightning_filter):
        city_index = np.arange(1, 50)
        data = read_wrf_small_radius_data(city_index)
        X_wrf, y_wrf = wrf_features_process(data, lightning_filter)
        feature_opt = 'satellite_aks'
        model, train = XGBoostModel.make_final_model(X_wrf, y_wrf['ho'], feature_opt)
        X_wrf.loc[:, 'wrf_pred'] = train.loc[:, 'y_pred']
        X_wrf.loc[:, 'wrf_oh'] = y_wrf['ho']
        omi_data = read_omi_small_radius_data(city_index)
        X_omi = omi_small_radius_qa4ecv_l3_features_process(omi_data, X_wrf)
        x, x, features = feature_selection(X_omi, X_omi, feature_opt)
        y_pred = XGBoostModel.make_predictions(model, x)
        X_omi.loc[:, 'omi_pred'] = y_pred
        X_omi.to_pickle(os.path.join(sat_aks_sr_output_path, 'ml_wrf_omi_filer_qa4ecv_l3.pkl'))

    @staticmethod
    def make_final_datasets_citeis(lightning_filter):
        city_index = np.arange(1, 50)
        for city_ind in city_index:
            data = read_wrf_small_radius_data([city_ind])
            X_wrf, y_wrf = wrf_features_process(data, lightning_filter)
            feature_opt = 'satellite_aks'
            model, train = XGBoostModel.make_final_model(X_wrf, y_wrf['ho'], feature_opt)
            X_wrf.loc[:, 'wrf_pred'] = train.loc[:, 'y_pred']
            X_wrf.loc[:, 'wrf_oh'] = y_wrf['ho']
            omi_data = read_omi_small_radius_data([city_ind])
            X_omi = omi_small_radius_features_process(omi_data, X_wrf)
            x, x, features = feature_selection(X_omi, X_omi, feature_opt)
            y_pred = XGBoostModel.make_predictions(model, x)
            X_omi.loc[:, 'omi_pred'] = y_pred
            save_filename = 'ml_wrf_omi_filer_{:02}.pkl'.format(city_ind)
            save_pathway = os.path.join(sat_aks_sr_output_path, 'cities')
            X_omi.to_pickle(os.path.join(save_pathway, save_filename))


    @staticmethod
    def collect_feature_importance():
        """
        Train model for all 49 cities and
        :param save_path: the save path
        :return:
        """
        save_path = sat_aks_sr_output_path
        lightning_filter = True
        city_index = np.arange(1, 50)
        importance = []
        rmse = []
        r2 = []
        for city_indx in city_index:
            print('Reading wrf data for city {:02}.'.format(city_indx))
            data = read_wrf_small_radius_data([city_indx])
            X_wrf, y_wrf = wrf_features_process(data, lightning_filter)
            feature_opt = 'satellite_aks'
            print('Training model for city {:02}.'.format(city_indx))
            model, train, test, features = XGBoostModel.make_model(X_wrf, y_wrf['ho'], feature_opt)
            this_importance = model.get_score(importance_type='gain')

            this_importance_value = [this_importance[k] for k in sat_aks_feature_sets]
            this_importance_value = np.array(this_importance_value)/sum(this_importance_value)
            importance.append(this_importance_value)
            rmse.append(mean_squared_error(test.loc[:, 'y_true'], test.loc[:, 'y_pred']))
            r2.append(r2_score(test.loc[:, 'y_true'], test.loc[:, 'y_pred']))
        np.save(os.path.join(save_path, 'feature_importance_cities'), np.array(importance))
        np.save(os.path.join(save_path, 'rmse_cities'), rmse)
        np.save(os.path.join(save_path, 'r2_cities'), r2)

    @staticmethod
    def find_best_k_for_cluster_analysis():
        """
        K-Mean cluster analysis to group cities with difference feature importance
        :param importances: importance array
        :return: kmeans
        """
        save_path = sat_aks_sr_output_path
        features = sat_aks_feature_sets
        importances = np.load(os.path.join(save_path, 'feature_importance_cities.npy'))
        kmeans = range(1, 20)
        distortions = []
        inertias = []
        for kmean in kmeans:
            kmeanmodel = KMeans(n_clusters=kmean).fit(importances)
            distortions.append(sum(np.min(cdist(importances, kmeanmodel.cluster_centers_,
                                                'euclidean'), axis=1)) / importances.shape[0])
            inertias.append(kmeanmodel.inertia_)
        plt.subplot(1, 2, 1)
        plt.plot(kmeans, distortions, 'bx-')
        plt.xlabel('Values of K')
        plt.ylabel('Distortion')
        plt.title('The Elbow Method using Distortion')
        plt.subplot(1, 2, 2)
        plt.plot(kmeans, inertias, 'bx-')
        plt.xlabel('Values of K')
        plt.ylabel('Inertia')
        plt.title('The Elbow Method using Inertia')
        plt.savefig(os.path.join(save_path, 'kmean_elbow.png'))

    @staticmethod
    def cluster_analysis_feature_importance():
        save_path = sat_aks_sr_output_path
        features = sat_aks_feature_sets
        lightning_filter = True
        city_index = np.arange(1, 50)
        importances = np.load(os.path.join(save_path, 'feature_importance_cities.npy'))
        kmeanmodel = KMeans(n_clusters=3).fit(importances)
        labels = kmeanmodel.labels_
        groups = []
        for label in np.unique(labels):
            groups.append(city_index[labels == label])

        importance_cluster = []
        city_cluster = []
        for label in np.unique(labels):
            city_locs = city_index[labels == label]
            data = read_wrf_small_radius_data(city_locs)
            X_wrf, y_wrf = wrf_features_process(data, lightning_filter)
            feature_opt = 'satellite_aks'
            print('Reading wrf data for label {:02}.'.format(label))
            model, train, test, features = XGBoostModel.make_model(X_wrf, y_wrf['ho'], feature_opt)
            print('Training model for label {:02}.'.format(label))
            this_importance = model.get_score(importance_type='gain')
            this_importance_value = [this_importance[k] for k in sat_aks_feature_sets]
            this_importance_value = np.array(this_importance_value)/sum(this_importance_value)

            importance_cluster.append(this_importance_value)
            city_cluster.append(city_locs[:])
        # np.save(os.path.join(save_path, 'city_clusters'), np.array(city_cluster))
        np.save(os.path.join(save_path, 'feature_importance_cluster_num3'), importance_cluster)
        transformed_city_cluster = np.zeros([len(city_cluster), len(max(city_cluster, key=lambda x: len(x)))])
        for i, j in enumerate(city_cluster):
            transformed_city_cluster[i][0:len(j)] = j
        np.save(os.path.join(save_path, 'city_clusters_num3'), transformed_city_cluster)

    @staticmethod
    def cluster_analysis_no2_importance():
        save_path = sat_aks_sr_output_path
        features = sat_aks_feature_sets
        lightning_filter = True
        city_index = np.arange(1, 50)
        importances = np.load(os.path.join(save_path, 'feature_importance_cities.npy'))
        no2_importances = []
        labels = []
        for importance in importances:
            no2_importances.append(importance[2])
            if importance[2] < 0.3:
                labels.append(0)
            elif importance[2] > 0.5:
                labels.append(2)
            else:
                labels.append(1)
        groups = []
        for label in np.unique(labels):
            groups.append(city_index[labels == label])

        importance_cluster = []
        city_cluster = []
        for label in np.unique(labels):
            city_locs = city_index[labels == label]
            data = read_wrf_small_radius_data(city_locs)
            X_wrf, y_wrf = wrf_features_process(data, lightning_filter)
            feature_opt = 'satellite_aks'
            print('Reading wrf data for label {:02}.'.format(label))
            model, train, test, features = XGBoostModel.make_model(X_wrf, y_wrf['ho'], feature_opt)
            print('Training model for label {:02}.'.format(label))
            this_importance = model.get_score(importance_type='gain')
            this_importance_value = [this_importance[k] for k in sat_aks_feature_sets]
            this_importance_value = np.array(this_importance_value) / sum(this_importance_value)

            importance_cluster.append(this_importance_value)
            city_cluster.append(city_locs[:])
        # np.save(os.path.join(save_path, 'city_clusters'), np.array(city_cluster))
        np.save(os.path.join(save_path, 'no2_importance_cluster_num3'), importance_cluster)
        transformed_city_cluster = np.zeros([len(city_cluster), len(max(city_cluster, key=lambda x: len(x)))])
        for i, j in enumerate(city_cluster):
            transformed_city_cluster[i][0:len(j)] = j
        np.save(os.path.join(save_path, 'city_clusters_no2_num3'), transformed_city_cluster)

    @staticmethod
    def cluster_analysis_no2_hcho_importance():
        save_path = sat_aks_sr_output_path
        features = sat_aks_feature_sets
        lightning_filter = True
        city_index = np.arange(1, 50)
        importances = np.load(os.path.join(save_path, 'feature_importance_cities.npy'))
        no2_importances = []
        hcho_importances = []
        labels = []
        for importance in importances:
            no2_importances.append(importance[2])
            hcho_importances.append(importance[3])
            if importance[2] >= 0.65:
                labels.append(0)
            elif importance[3] > 0.18:
                labels.append(1)
            elif importance[2] + importance[3] < 0.45:
                labels.append(2)
            else:
                labels.append(3)
        groups = []
        for label in np.unique(labels):
            groups.append(city_index[labels == label])

        importance_cluster = []
        city_cluster = []
        for label in np.unique(labels):
            city_locs = city_index[labels == label]
            data = read_wrf_small_radius_data(city_locs)
            X_wrf, y_wrf = wrf_features_process(data, lightning_filter)
            feature_opt = 'satellite_aks'
            print('Reading wrf data for label {:02}.'.format(label))
            model, train, test, features = XGBoostModel.make_model(X_wrf, y_wrf['ho'], feature_opt)
            print('Training model for label {:02}.'.format(label))
            this_importance = model.get_score(importance_type='gain')
            this_importance_value = [this_importance[k] for k in sat_aks_feature_sets]
            this_importance_value = np.array(this_importance_value) / sum(this_importance_value)

            importance_cluster.append(this_importance_value)
            city_cluster.append(city_locs[:])
        # np.save(os.path.join(save_path, 'city_clusters'), np.array(city_cluster))
        np.save(os.path.join(save_path, 'no2_hcho_importance_cluster_num4'), importance_cluster)
        transformed_city_cluster = np.zeros([len(city_cluster), len(max(city_cluster, key=lambda x: len(x)))])
        for i, j in enumerate(city_cluster):
            transformed_city_cluster[i][0:len(j)] = j
        np.save(os.path.join(save_path, 'city_clusters_no2_hcho_num4'), transformed_city_cluster)

    @staticmethod
    def cal_steady_state_inputs_daily(alphaeff, oh_var, nox_var, hcho_var, o3_var, model_type):
        city_index = np.arange(1, 50)
        feature_opt = 'satellite_aks'
        wrf_data = read_wrf_small_radius_data(city_index)
        X_wrf, y_wrf = wrf_features_process(wrf_data)
        y_wrf = y_wrf.loc[:, model_type]
        model, wrf = XGBoostModel.make_final_model(X_wrf, y_wrf, feature_opt)
        importance = model.feature_importances_
        X_wrf.loc[:, 'y_pred'] = wrf.loc[:, 'y_pred']
        X_wrf.loc[:, 'oh'] = y_wrf
        omi_data = read_omi_data(city_index)
        X_omi = omi_features_process(omi_data, X_wrf)
        y_pred = XGBoostModel.make_predictions(model, X_omi, feature_opt)
        X_omi.loc[:, 'oh_pred'] = y_pred * 1e6
        X_omi_p = pd.DataFrame()
        for city in city_index:
            this_omi, this_omi_annual = prepare_box_model_inputs_for_single_city(city, X_omi, o3_var, hcho_var)
            X_omi_p = pd.concat([X_omi_p, this_omi])
        X_omi_p = steady_state_cal(X_omi_p, oh_var, nox_var, hcho_var, alphaeff)
        save_name = 'box_model_output_01_to_49_annual_weekday_alpha_{}_{}_{}_{}_{:.2f}.csv'.format(
            oh_var, nox_var, hcho_var, o3_var, alphaeff)
        X_omi_p.to_csv(os.path.join('./Outputs/Boxmodel/', save_name))


if __name__ == '__main__':
    RunXGB.make_forward_feature_selection_shuffle(True)
    RunXGB.make_final_datasets_aks_update(True)
    RunXGB.sat_aks_general_model_aks_update(True)
    RunXGB.make_final_datasets_qa4ecv_l3_hcho(True)
    RunXGB.make_final_datasets_nasahcho(True)
    RunXGB.make_final_datasets(True)
    RunXGB.make_final_datasets_citeis(True)
    RunXGB.make_forward_feature_selection_improve(sat_iso_aks_sr_output_path, True)
    RunXGB.make_forward_feature_selection('satellite_iso_aks', sat_iso_aks_sr_output_path, True)
    RunXGB.sat_aks_noisy_general_model(True)
    RunXGB.cluster_analysis_no2_hcho_importance()
    RunXGB.cluster_analysis_no2_importance()
    RunXGB.cluster_analysis_feature_importance()
    RunXGB.find_best_k_for_cluster_analysis()
    RunXGB.collect_feature_importance()

    RunXGB.make_final_datasets(True)
    RunXGB.make_forward_feature_selection('satellite_aks', sat_aks_sr_output_path, True)
    RunXGB.make_cross_validation('satellite_aks', sat_aks_sr_output_path, True, 'ho')

    RunXGB.make_forward_feature_selection('oh_index', oh_index_sr_output_path, False)
    RunXGB.make_forward_feature_selection('oh_index', oh_index_sr_output_path, False)
    RunXGB.oh_index_general_model()
