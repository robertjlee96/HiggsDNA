from higgs_dna.workflows.dystudies import (
    DYStudiesProcessor,
    TagAndProbeProcessor,
)
from higgs_dna.workflows.taggers import taggers
from higgs_dna.workflows.Zmmy import ZmmyProcessor, ZmmyHist, ZmmyZptHist

workflows = {}

workflows["dystudies"] = DYStudiesProcessor
workflows["tagandprobe"] = TagAndProbeProcessor
workflows["zmmy"] = ZmmyProcessor
workflows["zmmyHist"] = ZmmyHist
workflows["zmmyZptHist"] = ZmmyZptHist

__all__ = ["workflows", "taggers", "DYStudiesProcessor"]
