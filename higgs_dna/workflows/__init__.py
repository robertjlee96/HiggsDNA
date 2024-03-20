from higgs_dna.workflows.dystudies import (
    DYStudiesProcessor,
    TagAndProbeProcessor,
)
from higgs_dna.workflows.taggers import taggers
from higgs_dna.workflows.HHbbgg import HHbbggProcessor
from higgs_dna.workflows.particleLevel import ParticleLevelProcessor
from higgs_dna.workflows.top import TopProcessor
from higgs_dna.workflows.Zmmy import ZmmyProcessor, ZmmyHist, ZmmyZptHist

workflows = {}

workflows["dystudies"] = DYStudiesProcessor
workflows["tagandprobe"] = TagAndProbeProcessor
workflows["HHbbgg"] = HHbbggProcessor
workflows["particleLevel"] = ParticleLevelProcessor
workflows["top"] = TopProcessor
workflows["zmmy"] = ZmmyProcessor
workflows["zmmyHist"] = ZmmyHist
workflows["zmmyZptHist"] = ZmmyZptHist

__all__ = ["workflows", "taggers", "DYStudiesProcessor"]
