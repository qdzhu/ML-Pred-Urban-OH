import xgboost as xgb
from RandomForest import feature_selection
from RandomForest import train_test_split
import numpy as np

def train_xgb_ml_model(x_train, y_train):
    """
    Train a random forest model, execute backward features eliminations if feature_opts is "full", otherwise
    train a normal random forest model with existing features
    :param x_train: x train dataset
    :param y_train: y train dataset
    :param feature_opts: options of feature sets
    :return: random forest model
    """
    dtrain = xgb.DMatrix(x_train, label=y_train)
    param = {'verbosity': 1, 'max_depth': 9, 'min_child_weight': 2, 'eta': 0.3}
    evallist = [(dtrain, 'train')]
    num_round = 50
    model = xgb.train(param, dtrain, num_round, evallist, early_stopping_rounds=5)
    return model


def param_tuning_xgb(x_train, y_train):
    max_depth_range = range(6, 10)
    min_child_weight_range = range(1, 5)
    eta_range = [0.1, 0.2, 0.3, 0.4]
    gridsearch_params = [(max_depth, min_child_weight, eta)
                         for max_depth in max_depth_range
                         for min_child_weight in min_child_weight_range
                         for eta in eta_range]
    dtrain = xgb.DMatrix(x_train, label=y_train)
    min_rmse = float('Inf')
    for max_depth, min_child_weight, eta in gridsearch_params:
        params = {'max_depth': max_depth, 'min_child_weight': min_child_weight,
                  'eta': eta, 'verbosity': 1}
        print("CV with max_depth={}, min_child_weight={}, eta={}".format(
            max_depth, min_child_weight, eta))
        num_boost_round = 50
        cv_results = xgb.cv(
            params,
            dtrain,
            num_boost_round=num_boost_round,
            seed=42,
            nfold=5,
            metrics='rmse',
            early_stopping_rounds=3
        )
        mean_rmse = cv_results['test-rmse-mean'].min()
        boost_rounds = cv_results['test-rmse-mean'].argmin()
        print("RMSE {} for {} rounds".format(mean_rmse, boost_rounds))
        if mean_rmse < min_rmse:
            min_mae = mean_rmse
            best_params = (max_depth, min_child_weight, eta)
    print(best_params)
    return best_params


def make_predictions(model, x):
    """
    Make oh predictions using the random forest model
    :param model: A trained random forest model
    :param x: features
    :return: model prediction
    """
    dtest = xgb.DMatrix(x)
    pred = model.predict(dtest)
    return pred


def make_model(X, y, feature_opts):
    """
    Make a generalization model
    :param cities: city indexed
    :param feature_opts: feature options
    :return: model, train data and test data
    """
    [x_train_unfilter, x_test_unfilter, y_train, y_test] = train_test_split(X, y, test_size=0.25, random_state=50)
    x_train, x_test, features = feature_selection(x_train_unfilter, x_test_unfilter, feature_opts)
    model = train_xgb_ml_model(x_train, y_train)
    x_train_unfilter.loc[:, 'y_true'] = y_train.values
    x_test_unfilter.loc[:, 'y_true'] = y_test.values
    x_train_unfilter.loc[:, 'y_pred'] = make_predictions(model, x_train)
    x_test_unfilter.loc[:, 'y_pred'] = make_predictions(model, x_test)
    return model, x_train_unfilter, x_test_unfilter, features


def make_final_model(X, y, feature_opts):
    """
        Make final model for prediction
        :param cities: city indexed
        :param feature_opts: feature options
        :return: model, train data and test data
        """
    x_train, x_test, features = feature_selection(X, X, feature_opts)
    model = train_xgb_ml_model(x_train, y)
    x_train.loc[:, 'y_true'] = y.values
    x_train.loc[:, 'y_pred'] = make_predictions(model, x_test)
    return model, x_train

def add_noise(x, param_noise):
    """
    Add noise term to the corresponding features
    :param x
    :param param_noise: noise dict
    :return:
    """
    for param, noise in param_noise.items():
        orig_values = x.loc[:, param].values
        gen_noise = np.zeros(np.size(orig_values))
        for i, value in enumerate(orig_values):
            gen_noise[i] = np.random.normal(0, value*noise, 1)[:]
        noise_values = orig_values + np.array(gen_noise)
        x.loc[:, param] = noise_values
    return x

def make_noise_model(X, y, feature_opts, para_noise):
    """
    Make a generalization model while noises are adding to the model
    :param X: features sets
    :param y: y
    :param feature_opts: feature options
    :param para_noise: a dict describing the parameter to add noise and the according uncertainty
    :return: model, train data and test data
    """
    [x_train_unfilter, x_test_unfilter, y_train, y_test] = train_test_split(X, y, test_size=0.25, random_state=50)
    x_train, x_test, features = feature_selection(x_train_unfilter, x_test_unfilter, feature_opts)
    model = train_xgb_ml_model(x_train, y_train)
    x_train_unfilter.loc[:, 'y_true'] = y_train.values
    x_test_unfilter.loc[:, 'y_true'] = y_test.values
    x_train_unfilter.loc[:, 'y_pred'] = make_predictions(model, x_train)
    x_test_unfilter.loc[:, 'y_pred'] = make_predictions(model, x_test)
    x_train = add_noise(x_train, para_noise)
    x_test = add_noise(x_test, para_noise)
    noise_model = train_xgb_ml_model(x_train, y_train)
    x_train_unfilter.loc[:, 'y_pred_noise'] = make_predictions(noise_model, x_train)
    x_test_unfilter.loc[:, 'y_pred_noise'] = make_predictions(noise_model, x_test)
    for param, noise in para_noise.items():
        noise_column = param+'_noise'
        x_train_unfilter.loc[:, noise_column] = x_train.loc[:, param].values
        x_test_unfilter.loc[:, noise_column] = x_test.loc[:, param].values
    return x_train_unfilter, x_test_unfilter, noise_model