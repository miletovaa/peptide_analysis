import os
import sys
import copy
from pathlib import Path
import yaml
from dotenv import load_dotenv


class DotDict(dict):
    __getattr__ = dict.get
    __setattr__ = dict.__setitem__


def _to_dotdict(x):
    return DotDict({k: _to_dotdict(v) if isinstance(v, dict) else v for k, v in x.items()})


def _detect_project_root():
    candidates = []
    try:
        nb_path = Path(__file__).resolve()
        candidates.extend([nb_path.parent, *nb_path.parents])
    except NameError:
        pass
    cwd = Path.cwd()
    candidates.extend([cwd, *cwd.parents])

    visited = set()
    for candidate in candidates:
        if not candidate:
            continue
        resolved = candidate.expanduser().resolve()
        if resolved in visited:
            continue
        visited.add(resolved)
        if (resolved / "params").is_dir() and (resolved / "scripts").is_dir():
            return resolved

    raise RuntimeError("Could not locate project root.")


PROJECT_ROOT = _detect_project_root()
os.chdir(PROJECT_ROOT)

if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

load_dotenv(PROJECT_ROOT / ".env")


def _load_params(params_id: str):
    candidates = sorted(Path("params").glob(f"{params_id}_*.yml"))
    if not candidates:
        raise FileNotFoundError(f"No params YAML found for ID {params_id} in ./params")
    with open(candidates[0], "r", encoding="utf-8") as f:
        return yaml.safe_load(f) or {}


def _deep_merge(base: dict, extra: dict):
    for key, value in extra.items():
        if (
            key in base
            and isinstance(base[key], dict)
            and isinstance(value, dict)
        ):
            _deep_merge(base[key], value)
        else:
            base[key] = copy.deepcopy(value)
    return base


PROJECT_ID = os.getenv("PROJECT_ID", "p02")
PROJECT_NAME = os.getenv("PROJECT_NAME", "peptide_analysis")
PARAMS_ID = os.getenv("PARAMS_ID", "000")
PARAMS_NAME = os.getenv("PARAMS_NAME", "default")
FILE_PATH_INPUT = os.getenv("FILE_PATH_INPUT", "./data")
FILE_PATH_OUTPUT = os.getenv("FILE_PATH_OUTPUT", "./results")

params_prot = _load_params(PARAMS_ID)

cfg = DotDict(
    {
        "FILE_PATH_INPUT": FILE_PATH_INPUT,
        "FILE_PATH_OUTPUT": FILE_PATH_OUTPUT,
        "PROJECT_ID": PROJECT_ID,
        "PROJECT_NAME": PROJECT_NAME,
        "PARAMS_ID": PARAMS_ID,
        "PARAMS_NAME": PARAMS_NAME,
        **_to_dotdict(params_prot),
    }
)
