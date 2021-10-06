from higgs_dna.workflows.dystudies import DYStudiesProcessor
from higgs_dna.workflows.taggers import taggers

workflows = {}

workflows["dystudies"] = DYStudiesProcessor

__all__ = ["workflows", "taggers", "DYStudiesProcessor"]
