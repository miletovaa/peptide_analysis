from typing import Literal
import os
import pandas as pd
from datetime import datetime

def create_output_dir(cfg):
    try:
        stamp = datetime.now().strftime("%Y%m%d-%H%M%S")  # avoid ':' for Windows
        OUTPUT_DIR = os.path.join(cfg.get("FILE_PATH_OUTPUT", "./results"), f"run_{cfg.get('PARAMS_ID', '000')}_{stamp}")
        os.makedirs(OUTPUT_DIR, exist_ok=True)
    except Exception as e:
        print(f"Failed to create output directory: {e}")
    return OUTPUT_DIR


OutputType = Literal[
    "raw",
    "prep",
    "norm",
    "stats",
    "ml",
    "viz",
    "fig",
    "output",
    "report",
]

def get_output_file_name(step: str, type: OutputType, description: str, ext: str, cfg, output_dir: str) -> str:
    date = pd.Timestamp.now().strftime("%Y%m%d")
    return os.path.join(output_dir, f"{cfg.get("PROJECT_ID")}_{cfg.get("PROJECT_NAME")}_{cfg.get("PARAMS_ID")}_{cfg.get("PARAMS_NAME")}_{step}_{type}_{description}_{date}.{ext}")


def save_figure(fig, step: str, description: str, cfg, output_dir: str):
    file_name = get_output_file_name(step, "fig", description, "png", cfg, output_dir)
    fig.write_image(file_name)
    fig.write_html(file_name.replace(".png", ".html"))
    print(f"Figure saved to {file_name}")