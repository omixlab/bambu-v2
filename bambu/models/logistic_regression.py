from flaml.model import SKLearnEstimator
from flaml import tune
from sklearn.linear_model import LogisticRegression
import numpy as np
import shutil
import time

class LogisticRegressionEstimator(SKLearnEstimator):
    """The class for tuning a Logistic Regression Classifier"""

    @classmethod
    def search_space(cls, data_size, **params):

        search_space = {
            "dual": {
                "domain": tune.choice([True, False]),
                "init_value": False
            },
        
            "tol": {
                "domain": tune.uniform(lower=1e-5, upper=1e5),
                "init_value": 1e-5
            },
            "C": {
                "domain": tune.uniform(lower=1e-5, upper=1e5),
                "init_value": 1e-5
            },
            "fit_intercept": {
                "domain": tune.choice([True, False]),
                "init_value": False
            },
            "intercept_scaling": {
                "domain": tune.uniform(lower=1e-5, upper=1e5),
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
        lr = LogisticRegression(**hyperparameters)
        start_time = time.time()
        deadline = start_time + budget if budget else np.inf
        lr.fit(X_train, y_train)        
        train_time = time.time() - start_time
        self._model = lr
        return train_time

    def predict(self, X_test):
        return super().predict(X_test)

    def predict_proba(self, X_test):
        return super().predict_proba(X_test)


