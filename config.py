from pathlib import Path

try:
    PROJECT_ROOT = Path(__file__).resolve().parent
except NameError:
    PROJECT_ROOT = Path.cwd()



# from pathlib import Path
# import sys
# import subprocess

# # ── Project Root ──────────────────────────────────────────────
# PROJECT_ROOT = Path(__file__).resolve().parent

# # ── Install dependencies ──────────────────────────────────────
# requirements = PROJECT_ROOT / "python_requirements.txt"
# subprocess.check_call([sys.executable, "-m", "pip", "install", "-r", str(requirements)])
