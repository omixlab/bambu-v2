from mordred import Calculator, descriptors
from rdkit.Chem import Descriptors

class DescriptorsPreprocessor:

    def __init__(self):
        self.calculator = Calculator(descriptors)
        self.features = [descriptor for descriptor in self.calculator.descriptors]

    def compute_features(self, mol):

        df = self.calculator.pandas([mol], quiet=True)
        descriptors = {column:df.iloc[0][column] for column in df.columns}
        return descriptors