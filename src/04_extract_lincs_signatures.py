import pandas as pd
  from cmapPy.pandasGEXpress.parse import parse
  import os

  print("Installing cmapPy via pip if missing...")
  import subprocess
  subprocess.call(["pip", "install", "cmapPy"])

  # The gctx file needs to be downloaded since it's 20GB. Let's do it outside this script using curl/wget.
  