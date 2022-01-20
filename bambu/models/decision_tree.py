from flaml.model import SKLearnEstimator
from flaml import tune
from sklearn.tree import DecisionTreeClassifier
import numpy as np
import time

class DecisionTreeEstimator(SKLearnEstimator):
    """The class for tuning a Decision Tree Classifier"""

    @classmethod
    def search_space(cls, data_size, **params):
       
        search_space = {
            "criterion": {
                "domain": tune.choice(["gini", "entropy"]),
                "init_value": "gini"
            },
            "max_depth": {
                "domain": tune.randint(lower=4, upper=200),
                "init_value": 4
            },
            "splitter": {
                "domain": tune.choice(["best", "random"]),
                "init_value": "best"
            },
            "min_samples_split": {
                "domain": tune.randint(lower=2, upper=data_size[0]),
                "init_value": 2
            },
            "min_samples_leaf": {
                "domain": tune.randint(lower=2, upper=data_size[0]),
                "init_value": 2
            },
            "min_weight_fraction_leaf": {
                "domain": tune.uniform(lower=1e-5, upper=1),
                "init_value": 1e-5
            },
            "max_features": {
                "domain": tune.choice([None, "auto", "sqrt", "log2"]),
                "init_value": None                
            },
            "max_leaf_nodes": {
                "domain": tune.randint(lower=2, upper=data_size[0]),
                "init_value": 2             
            },
            "min_impurity_decrease": {
                "domain": tune.uniform(lower=1e-5, upper=1.0),
                "init_value": 1e-5  
            },
            "ccp_alpha": {
                "domain": tune.uniform(lower=1e-5, upper=1e3),
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
        dtc = DecisionTreeClassifier(**hyperparameters)
        start_time = time.time()
        deadline = start_time + budget if budget else np.inf
        dtc.fit(X_train, y_train)        
        train_time = time.time() - start_time
        self._model = dtc
        return train_time

    def predict(self, X_test):
        return super().predict(X_test)

    def predict_proba(self, X_test):
        return super().predict_proba(X_test)


