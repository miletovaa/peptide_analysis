import os
from pathlib import Path
import yaml
from dotenv import load_dotenv


class DotDict(dict):
    __getattr__ = dict.get
    __setattr__ = dict.__setitem__


def _to_dotdict(x):
    return DotDict({k: _to_dotdict(v) if isinstance(v, dict) else v for k, v in x.items()})


load_dotenv()

FILE_PATH_INPUT = os.getenv("FILE_PATH_INPUT", "./data")
FILE_PATH_OUTPUT = os.getenv("FILE_PATH_OUTPUT", "./results")
PARAMS_ID = os.getenv("PARAMS_ID", "000")

yml_candidates = list(Path("params").glob(f"{PARAMS_ID}_*.yml"))
if not yml_candidates:
    raise FileNotFoundError(f"No params YAML found for ID {PARAMS_ID} in ./params")
with open(yml_candidates[0], "r", encoding="utf-8") as f:
    yml_cfg = _to_dotdict(yaml.safe_load(f) or {})

cfg = DotDict(
    {
        "FILE_PATH_INPUT": FILE_PATH_INPUT,
        "FILE_PATH_OUTPUT": FILE_PATH_OUTPUT,
        "PARAMS_ID": PARAMS_ID,
        **yml_cfg,
    }
)