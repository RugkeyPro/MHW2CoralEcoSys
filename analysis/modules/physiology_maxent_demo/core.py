"""Native MaxEnt engine helpers and the accepted Acropora projection constraint."""
import hashlib
import os
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

VARIABLES = ['bathymetry', 'rugosity', 'bottomT', 'thetao', 'si', 'so', 'zos']


def coral_multiplier(temperature, omega):
    """Accepted structural assumptions: asymmetric TPC times carbonate sigmoid."""
    temperature = np.asarray(temperature, dtype=float)
    omega = np.asarray(omega, dtype=float)
    sigma = np.where(temperature <= 27.0, 5.0, 2.0)
    phi = np.exp(-((temperature - 27.0) ** 2) / (2 * sigma ** 2))
    phi *= 1 / (1 + np.exp(-4.0 * (omega - 3.0)))
    return np.clip(phi, 0, 1) ** 0.5


def split_spatial(frame):
    """Deterministic 5-degree blocks, shared between presences and background."""
    longitude_block = np.floor((frame.longitude + 180) / 5).astype(int)
    latitude_block = np.floor((frame.latitude + 90) / 5).astype(int)
    return ((longitude_block * 73 + latitude_block * 37) % 5 == 0).to_numpy()


def verify_sources(root, manifest):
    for record in manifest['packaged_files']:
        digest = hashlib.sha256((root / record['file']).read_bytes()).hexdigest()
        if digest != record['sha256']:
            raise ValueError(f"Input integrity failed: {record['file']}")


def java_executable():
    executable = shutil.which('java')
    if executable:
        return executable
    if os.name == 'nt':
        candidates = sorted(Path('C:/Program Files/Java').glob('*/bin/java.exe'))
        if candidates:
            return str(candidates[-1])
    raise RuntimeError('Install Java 17 or newer and place java on PATH.')


def run_java(jar, arguments, log):
    command = [java_executable(), '-Xmx2g', '-Djava.awt.headless=true', '-cp', str(jar)] + arguments
    result = subprocess.run(command, capture_output=True, text=True, timeout=300)
    log.write_text(result.stdout + '\n' + result.stderr, encoding='utf-8')
    if result.returncode:
        raise RuntimeError(f'MaxEnt failed with code {result.returncode}; inspect {log.name}')
    return result.stdout


def swd(frame, species, destination):
    result = frame[['longitude', 'latitude'] + VARIABLES].copy()
    result.insert(0, 'species', species)
    result.to_csv(destination, index=False, float_format='%.17g')


def auc(presence_scores, background_scores):
    """Presence-background rank AUC, not presence/true-absence accuracy."""
    positive = np.asarray(presence_scores)
    negative = np.asarray(background_scores)
    ranks = pd.Series(np.r_[positive, negative]).rank(method='average').to_numpy()
    return float((ranks[:len(positive)].sum() - len(positive) * (len(positive) + 1) / 2)
                 / (len(positive) * len(negative)))


def project(jar, model, frame, output, label):
    """Use the native engine's SWD projection, with the parent projection flags."""
    valid = frame[VARIABLES].notna().all(axis=1)
    predictors = output / f'{label}_environment.csv'
    swd(frame.loc[valid], 'projection', predictors)
    destination = output / f'{label}_maxent.csv'
    run_java(jar, ['density.Project', str(model), str(predictors), str(destination),
                   'doclamp=true', 'extrapolate=false', 'fadebyclamping=true', 'outputformat=logistic'],
             output / f'{label}_projection.log')
    predictions = pd.read_csv(destination)
    print(f'{label}: native output columns {list(predictions.columns)}')
    return valid, predictions
