from mol2vec.features import mol2alt_sentence, MolSentence, sentences2vec
from gensim.models import word2vec
import numpy as np

class Mol2VecPreprocessor:

    def __init__(self, pretrained_model):
        self.pretrained_model = pretrained_model
        self.model = word2vec.Word2Vec.load(pretrained_model)
        self.features = [f"vector_{v}" for v in range(self.model.vector_size)]

    def compute_features(self,mol):
        sentence = MolSentence(mol2alt_sentence(mol, 1))
        vectors  = sentences2vec([sentence], self.model, unseen='UNK')[0]
        return {f"vector_{v}":value for v, value in enumerate(vectors)}