"""Execute with the current Python, without relying on a global kernelspec."""
import argparse
import os
from pathlib import Path
import sys
import nbformat
from nbclient import NotebookClient
from jupyter_client import KernelManager

parser=argparse.ArgumentParser()
parser.add_argument('--output',type=Path,required=True)
args=parser.parse_args()
root=Path(__file__).resolve().parents[1]
notebook=nbformat.read(root/'analysis.ipynb',as_version=4)
# A fresh subprocess kernel uses exactly this interpreter, including a venv.
km=KernelManager(kernel_name='python3')
km.kernel_spec.argv=[sys.executable,'-m','ipykernel_launcher','-f','{connection_file}']
client=NotebookClient(notebook,km=km,timeout=300,startup_timeout=60,resources={'metadata':{'path':str(root)}})
try:
    notebook=client.execute()
finally:
    if client.kc is not None:
        client.kc.stop_channels()
    if km.has_kernel:
        km.shutdown_kernel(now=True)
        km.cleanup_resources()
args.output.parent.mkdir(parents=True,exist_ok=True)
nbformat.write(notebook,args.output)
print(f'Executed all notebook cells: {args.output}')
