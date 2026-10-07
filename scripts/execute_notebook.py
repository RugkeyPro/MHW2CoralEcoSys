"""Execute every notebook cell in the current validated Python environment."""
import argparse
import json
import sys
import tempfile
from pathlib import Path

import nbformat
from nbclient import NotebookClient
from jupyter_client import KernelManager
from jupyter_client.kernelspec import KernelSpecManager
from ipykernel.kernelspec import write_kernel_spec

ROOT = Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--output', type=Path, default=ROOT / 'outputs/notebook_execution/MHW2CoralEcoSys_demo.executed.ipynb')
parser.add_argument('--inplace', action='store_true', help='Save successful real outputs into the primary notebook')
args = parser.parse_args()
notebook_path = ROOT / 'MHW2CoralEcoSys_demo.ipynb'
notebook = nbformat.read(notebook_path, as_version=4)

with tempfile.TemporaryDirectory(prefix='mhw2coral-kernel-') as directory:
    kernel_path = Path(directory) / 'mhw2coral-demo'
    write_kernel_spec(path=str(kernel_path))
    manager = KernelManager(kernel_name='mhw2coral-demo',
                            kernel_spec_manager=KernelSpecManager(kernel_dirs=[directory]))
    client = NotebookClient(notebook, km=manager, timeout=600, allow_errors=False,
                            resources={'metadata': {'path': str(ROOT)}})
    try:
        executed = client.execute()
    finally:
        if manager.has_kernel:
            manager.shutdown_kernel(now=True)

code_cells = [cell for cell in executed.cells if cell.cell_type == 'code']
if any(cell.execution_count is None for cell in code_cells):
    raise RuntimeError('A code cell was not executed')
if any(output.output_type == 'error' for cell in code_cells for output in cell.outputs):
    raise RuntimeError('Notebook contains an execution error')
checks = json.loads((ROOT / 'outputs/notebook_demo/run_checks.json').read_text(encoding='utf-8'))
if checks.get('status') != 'passed' or not checks.get('full_frozen_design'):
    raise RuntimeError('Notebook did not complete the full accepted scientific design')
if checks['environment']['numpy'] != '2.3.5' or checks['environment']['pandas'] != '2.3.3':
    raise RuntimeError('Notebook did not use the pinned scientific environment')
executed.metadata.kernelspec = {'display_name': 'Python 3 (ipykernel)', 'language': 'python', 'name': 'python3'}
nbformat.validate(executed)
target = notebook_path if args.inplace else args.output
target.parent.mkdir(parents=True, exist_ok=True)
nbformat.write(executed, target)
print(json.dumps({'notebook': target.name, 'executed_code_cells': len(code_cells),
                  'status': 'passed', 'draw_records': checks['draw_rows'],
                  'python_executable': sys.executable}, indent=2))
