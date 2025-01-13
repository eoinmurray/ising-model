import os
import json

data_folder = os.path.join(os.path.dirname(__file__), 'data')

try:
  files = os.listdir(data_folder)
except OSError as e:
  files = []

directories = []
if files:
  for f in files:
    if os.path.isdir(os.path.join(data_folder, f)):
      directories.append({
        'name': f,
        'path': f
      })

print(
  json.dumps(
    {
      'directories': directories
    }, 
    indent=2
  )
)