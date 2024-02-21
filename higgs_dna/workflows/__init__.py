from higgs_dna.workflows.dystudies import (
    DYStudiesProcessor,
    TagAndProbeProcessor,
)
from higgs_dna.workflows.taggers import taggers
from higgs_dna.workflows.Zmmy import ZmmyProcessor, ZmmyHist, ZmmyZptHist
from higgs_dna.workflows.top import TopProcessor

workflows = {}

workflows["dystudies"] = DYStudiesProcessor
workflows["tagandprobe"] = TagAndProbeProcessor
workflows["zmmy"] = ZmmyProcessor
workflows["zmmyHist"] = ZmmyHist
workflows["zmmyZptHist"] = ZmmyZptHist
workflows["top"] = TopProcessor

__all__ = ["workflows", "taggers", "DYStudiesProcessor"]
