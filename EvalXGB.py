import XGBoostModel
import numpy as np
from Utils import *
from sklearn.metrics import mean_squared_error
from sklearn.metrics import r2_score
from sklearn.cluster import KMeans
from scipy.spatial.distance import cdist
import matplotlib.pyplot as plt
import os



class EvalXGB:
    @staticmethod
    def noise_rf(feature_opt, param_noise, iterations, X, y):
        rmse = []
        for iter in range(iterations):
            print('Iteration: {}'.format(iter))
            train, test = XGBoostModel.make_noise_model(X, y['ho'], feature_opt, param_noise)
            rmse.append(mean_squared_error(test.loc[:, 'y_true'], test.loc[:, 'y_pred_noise']))
        return rmse

    @staticmethod
    def run_noise_test_rf_a_param(feature_opt, param_name, iterations):
        lightning_filter = True
        city_index = np.arange(1, 50)
        data = read_wrf_small_radius_data(city_index)
        X, y = wrf_features_process(data, lightning_filter)

        uncerts = np.linspace(0, 1.05, 11)
        rmses = []
        for uncert in uncerts:
            param_noise = dict({param_name: uncert})
            print(
                'Training a generalization model with {} noise in {}.'.format(param_noise.keys(), param_noise.values()))
            rmse = EvalXGB.noise_rf(feature_opt, param_noise, iterations, X, y)
            rmses.append(rmse)
        save_path = sat_aks_sr_output_path
        save_name = 'rmse_' + param_name + '_noise'
        np.save(os.path.join(save_path, save_name), rmses)

    @staticmethod
    def run_noise_test_rf_two_params(feature_opt, param_names):
        lightning_filter = True
        city_index = np.arange(1, 50)
        data = read_wrf_small_radius_data(city_index)
        X, y = wrf_features_process(data, lightning_filter)

        uncerts = np.linspace(0, 1.05, 22)
        rmses = np.zeros([np.size(uncerts), np.size(uncerts)])
        iterations = 1
        for i_p1, uncert_p1 in enumerate(uncerts):
            for i_p2, uncert_p2 in enumerate(uncerts):
                param_noise = dict({param_names[0]: uncert_p1, param_names[1]: uncert_p2})
                rmse = EvalXGB.noise_rf(feature_opt, param_noise, iterations, X, y)
                rmses[i_p1, i_p2] = rmse[0]
        save_path = sat_aks_sr_output_path
        save_name = 'rmse_noise_params_{}_{}'.format(*param_names)
        np.save(os.path.join(save_path, save_name), rmses)

    @staticmethod
    def run_a_noise_test(feature_opt, param_noise, iterations):
        rmse = EvalXGB.noise_rf(feature_opt, param_noise, iterations)
        std = np.std(rmse)


def main():
    feature_opt = 'satellite_aks'
    param_names = ['no2_vcd', 'hcho_vcd_aks', 'PHOTR_O31D', 'QVAPOR', 'pressure', 'temperature']
    iterations = 1
    for param_name in param_names:
        EvalXGB.run_noise_test_rf_a_param(feature_opt, param_name, iterations)

    param_names = ['no2_vcd', 'hcho_vcd_aks']
    EvalXGB.run_noise_test_rf_two_params(feature_opt, param_names)



if __name__ == '__main__':
    main()
    #EvalXGB.run_a_noise_test('default', dict({'hcho_vcd': 0.2}), 10)