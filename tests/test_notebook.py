"""Check the primary notebook is a staged analysis rather than a CLI wrapper."""
import unittest
from pathlib import Path

import nbformat


class NotebookStructureTests(unittest.TestCase):
    def test_staged_scientific_entrypoint(self):
        root = Path(__file__).resolve().parents[1]
        notebook = nbformat.read(root / 'MHW2CoralEcoSys_demo.ipynb', as_version=4)
        nbformat.validate(notebook)
        cells = [cell for cell in notebook.cells if cell.cell_type == 'code']
        self.assertEqual(len(cells), 8)
        source = '\n'.join(cell.source for cell in cells)
        for cell in cells:
            compile(cell.source, '<notebook>', 'exec')
        for operation in ('core.point_estimates(', 'core.bootstrap(', 'core.summarize(', 'compare_frames(', 'render_figures('):
            self.assertIn(operation, source)
        self.assertNotIn('subprocess', source)
        self.assertNotIn('run_demo.py', source)


if __name__ == '__main__':
    unittest.main()
