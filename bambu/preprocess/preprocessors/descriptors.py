from rdkit.Chem import Descriptors

class DescriptorsPreprocessor:

    def __init__(self):
        self.features = list(dict(Descriptors._descList).keys())

    def compute_features(self, mol):
        descriptors = {}
        for descriptor_name, descriptor_function in Descriptors._descList:
            descriptors[descriptor_name] = descriptor_function(mol)
        return descriptors