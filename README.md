# 📓 Universal template for lab project management

In order to remain consistent according to the [SOPs protocol](link_to_be_added) use this repository template

## 💻 Code structure

Following file layout is a default maximum setup *(meaning the majority of the directories and files may not be needed for ceratain projects)*.

```
[project-name]/
    ├─ .vscode/
    │  └─ extensions.json          # recommendations (e.g. Jupyter, GitLens, etc.)
    ├─ data/                       # default directory of input files; NEVER commit
    ├─ devcontainer/               # VS Code Dev Containers
    │  ├─ devcontainer.json
    │  └─ Dockerfile
    ├─ env/
    │  ├─ conda-lock.yml           # exact, reproducible lock
    │  ├─ environment.yml          # main conda env
    │  └─ requirements.txt         # for pip-only contexts
    ├─ logs/
    │  ├─ app.log
    │  └─ errors.log
    ├─ notebooks/                  # Jupiter/Quarto/etc.
    ├─ params/                     # experiment parameters
    │  └─ 000_default.yaml
    ├─ results/                    # default output directory
    │  ├─ reports/                 # .md, .pdf, .html, etc.
    │  └─ [run]/                   # create new run_[parameters-yaml-id]_[date-time] for each run
    │     ├─ docs/                 # .doc, .pdf, .txt, etc.
    │     ├─ tables/               # .csv, .xlsx, etc.
    │     └─ figures/              # .png, .jpg, .svg, etc.
    ├─ scripts/                    # small runnable scripts (shell/py)
    │  ├─ script.py
    │  └─ script.sh
    ├─ src/
    │  └─ project_name/
    │     ├─ utils/                # helpers
    │     ├─ __init__.py
    │     └─ main.py/              # entrypoint
    ├─ tests/
    │  └─ test.py
    ├─ workflows/                  # pipelines
    │  └─ snakemake/               # Snakefiles, rule-level envs
    │     ├─ Snakefile
    │     └─ envs/
    ├─ .gitignore                  # editable template; files that are ignored by Git
    ├─ .gitattributes              # LFS patterns + notebook filters
    ├─ .env.example                # template for paths/keys/etc.
    ├─ .env                        # local only, gitignored
    ├─ DOCS.md                     # template for description, documentation and installation 
    └─ README.md                   # default repository layout guide
```
