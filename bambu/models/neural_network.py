from flaml.model import SKLearnEstimator
from flaml import tune
from sklearn.neural_network import MLPClassifier
import numpy as np
import itertools
import shutil
import time

class NeuralNetworkEstimator(SKLearnEstimator):
    """The class for tuning a Neural Network Classifier"""

    @classmethod
    def search_space(cls, data_size, **params):

        NETWORK_CONFIGURATIONS = []
        MIN_NUMBER_OF_HIDDEN_LAYERS = 2
        MIN_NUMBER_OF_HIDDEN_LAYERS = 4
        MIN_NUMBER_OF_NEURONS_IN_HIDDEN_LAYERS = 1
        MAX_NUMBER_OF_NEURONS_IN_HIDDEN_LAYERS = 100

        for hidden_layers in range(MIN_NUMBER_OF_HIDDEN_LAYERS, MIN_NUMBER_OF_HIDDEN_LAYERS+1):
            for configuration in itertools.permutations(range(MIN_NUMBER_OF_NEURONS_IN_HIDDEN_LAYERS, MAX_NUMBER_OF_NEURONS_IN_HIDDEN_LAYERS), hidden_layers):
                NETWORK_CONFIGURATIONS.append(configuration)

        search_space = {
            "hidden_layer_sizes": {
                "domain": tune.choice(NETWORK_CONFIGURATIONS),
                "init_value": NETWORK_CONFIGURATIONS[0]
            },
        
            "activation": {
                "domain": tune.choice(["identity", "logistic", "tanh", "relu"]),
                "init_value": "relu"
            },
            "solver": {
                "domain": tune.choice(["lbfgs", "sgd", "adam"]),
                "init_value": "adam"
            },
            "alpha": {
                "domain": tune.uniform(lower=1e-5, upper=1e-2),
                "init_value": 1e-5
            },
            "batch_size": {
                "domain": tune.randint(lower=1, upper=data_size[0]),
                "init_value": 1
            },
            "learning_rate": {
                "domain": tune.choice(["constant", "ïnvscaling", "adaptative"]),
                "init_value": "constant"
            },
            "learning_rate_init": {
                "domain": tune.uniform(lower=1e-5, upper=1e-2),
                "init_value": 1e-5
            },
            "momentum": {
                "domain": tune.uniform(lower=1e-5, upper=1 - 1e-5),
                "init_value": 1e-5
            }
        }

        cls._hyperameters = search_space.keys()
        return search_space
    
    @classmethod
    def size(cls, config):
        return 1.0

    def config2params(self, config: dict) -> dict:
        params = config.copy()
        return params

    def __init__(
        self,
        task="classification",
        **config,
    ):
        super().__init__(task, **config)

    def fit(self, X_train, y_train, budget=None, **kwargs):
        hyperparameters = self.params.copy()
        for param in self.params:
            if param not in self._hyperameters:
                del hyperparameters[param]
        nn = MLPClassifier(**hyperparameters)
        start_time = time.time()
        deadline = start_time + budget if budget else np.inf
        nn.fit(X_train, y_train)        
        train_time = time.time() - start_time
        self._model = nn
        return train_time

    def predict(self, X_test):
        return super().predict(X_test)

    def predict_proba(self, X_test):
        return super().predict_proba(X_test)


