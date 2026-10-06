export type Metric = { mean: number; min: number; max: number };
export type Future = {
  scenario: string; region: string; quality: Metric; burden: Metric; pressure: Metric;
  historicalPressure: number; futurePressure: number; futureHabitat: number;
  members: { member: string; quality: number; burden: number; pressure: number }[];
};
export type Conservation = { scenario: string; region: string; reallocation: Metric; gain: Metric };
export type History = { region: string; taxon: string; group: number; exposure: number; change: number; min: number; max: number };
export type Source = { file: string; project_relative_source: string; sha256: string; bytes: number; rows?: number; columns?: string[]; kind: string; unchanged_copy: boolean };
export type DemoData = {
  schemaVersion: number; lineageDate: string; curatedDate: string; future: Future[];
  conservation: Conservation[]; historical: History[];
  provenance: { files: Source[]; selection: string; historical_lineage: string; future_lineage: string };
  checks: { sourceHashes: number; numericComparisons: number; ensembleRows: number; conservationRows: number; historicalRows: number };
};
export type MapCell = [number, number, number, number];
