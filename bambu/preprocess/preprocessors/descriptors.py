from mordred import Calculator, descriptors
from rdkit.Chem import Descriptors

class DescriptorsPreprocessor:

    def __init__(self):
        self.calculator = Calculator(descriptors)
        self.features = None

    def compute_features(self, mol):

        df = self.calculator.pandas([mol], quiet=True)
        if self.features is None:
            self.features = list(df.columns)
        descriptors = {column:df.iloc[0][column] for column in df.columns}
        for descriptor in descriptors:
           try:
               descriptors[descriptor] = float(descriptors[descriptor])
           except:
               descriptors[descriptor] = None
        return descriptors
