import re
from typing import Dict


def get_truth_info_dict(filename: str) -> Dict[str, float]:

    # Get the dictonary of parameters and their values searching in the filename

    param_values = {}

    patterns = {
        "kl": r"kl[-_](m?[\d]+p?[\d]*)",
        "kt": r"kt[-_](m?[\d]+p?[\d]*)",
        "c2": r"c2[-_](m?[\d]+p?[\d]*)",
        "CV": r"CV[-_](m?[\d]+p?[\d]*)",
        "C2V": r"C2V[-_](m?[\d]+p?[\d]*)",
        "C3": r"C3[-_](m?[\d]+p?[\d]*)",
        "Radion_M": r"RadiontoHHto2B2G_M-([\d]+)",
        "BulkGraviton_M": r"BulkGravitontoHHto2B2G_M-([\d]+)"
    }
    # mY can be added later when the sample requests are submitted

    for key, pattern in patterns.items():
        match = re.search(pattern, filename)
        if match:
            param_values[key] = float(match.group(1).replace("m", "-").replace("p", "."))
    return param_values
