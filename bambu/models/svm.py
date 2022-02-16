from flaml.model import SKLearnEstimator
from flaml import tune
from sklearn.svm import SVC
import numpy as np
import shutil
import time

class SvmEstimator(SKLearnEstimator):
    """The class for tuning a Support Vector Machine Classifier"""

    @classmethod
    def search_space(cls, data_size, **params):

        search_space = {
            "C": {
                "domain": tune.uniform(lower=1e-5, upper=1e5),
                "init_value": 1e-5
            },
            "kernel": {
                "domain": tune.choice(["linear", "poly", "rbf", "sigmoid"]),
                "init_value": "rbf"
            },
            "degree": {
                "domain": tune.randint(lower=2, upper=10),
                "init_value": 2
            },
            "gamma": {
                "domain": tune.choice(["scale", "auto"]),
                "init_value": "auto"
            },
            "coef0": {
                "domain": tune.uniform(lower=1e-5, upper=1e3),
                "init_value": 1e-5
            },
            "tol": {
                "domain": tune.uniform(lower=1e-5, upper=1e5),
                "init_value": 1e-5
            },
            "shrinking": {
                "domain": tune.choice([False, True]),
                "init_value": False
            },
            "cache_size": {
                "domain": tune.uniform(lower=1e-1, upper=5e2),
                "init_value": 1e-1                
            },
            "max_iter": {
                "domain": tune.randint(lower=-1, upper=10),
                "init_value": -1               
            },
            "break_ties": {
                "domain": tune.choice([False, True]),
                "init_value": False                  
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
        svm = SVC(**hyperparameters, probability=True)
        start_time = time.time()
        deadline = start_time + budget if budget else np.inf
        svm.fit(X_train, y_train)        
        train_time = time.time() - start_time
        self._model = svm
        return train_time

    def predict(self, X_test):
        return super().predict(X_test)

    def predict_proba(self, X_test):
        return super().predict_proba(X_test)


