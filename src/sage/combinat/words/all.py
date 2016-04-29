r"""
Word features that are imported by default in the interpreter namespace
"""
from alphabet import Alphabet, build_alphabet
from morphism import WordMorphism
from paths import WordPaths
from word import Word
from word_options import WordOptions
from word_generators import words
from words import Words, FiniteWords, InfiniteWords
from free_group_automorphism import FreeGroupMorphism
from free_group_automorphism import FreeGroupAutomorphism
from free_group_automorphism import free_group_automorphisms
from free_group_word import FreeGroupWord
from inverse_alphabet import AlphabetWithInverses
from inverse_graph import GraphWithInverses
from inverse_graph import MetricGraph
from free_group import FreeGroup
from marked_graph import MarkedGraph
from marked_graph import MarkedMetricGraph
from graph_map import GraphMap
from graph_self_map import GraphSelfMap
from train_track_map import TrainTrackMap
from convex_core import ConvexCore