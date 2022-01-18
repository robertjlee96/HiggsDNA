from higgs_dna.workflows.dystudies import (
    DYStudiesProcessor,
    TagAndProbeProcessor,
)
from higgs_dna.workflows.taggers import taggers

workflows = {}

workflows["dystudies"] = DYStudiesProcessor
workflows["tagandprobe"] = TagAndProbeProcessor

__all__ = ["workflows", "taggers", "DYStudiesProcessor"]
