import os
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from alsi_to_microsim import (  # noqa: E402
    build_phase_functions,
    find_eutectic,
    load_projected_database,
    mole_fraction_to_wt_percent_si,
    temperature_grid,
    validate_reference_eutectic,
    wt_percent_to_mole_fraction_si,
    TARGET_PHASES,
)


@pytest.mark.parametrize("wt_percent", [0.1, 1.6, 3.5, 7.5, 12.6, 99.9])
def test_weight_mole_round_trip(wt_percent):
    recovered = mole_fraction_to_wt_percent_si(wt_percent_to_mole_fraction_si(wt_percent))
    assert recovered == pytest.approx(wt_percent, abs=1e-12)


def test_temperature_grid_clips_to_eutectic_and_refines():
    values = temperature_grid(850.0, 853.0, 1.0, 850.2)
    assert values.tolist() == pytest.approx([850.2, 850.3, 850.4, 850.7, 851.0, 852.0, 853.0])


@pytest.mark.skipif("ALSI_TDB" not in os.environ, reason="Set ALSI_TDB for private integration test")
def test_private_database_eutectic():
    db = load_projected_database(Path(os.environ["ALSI_TDB"]))
    functions = {phase: build_phase_functions(db, phase) for phase in TARGET_PHASES}
    eutectic, residual = find_eutectic(functions)
    errors, _ = validate_reference_eutectic(eutectic, 850.0, 12.6, 1.6)
    assert residual < 1e-3
    assert max(errors.values()) < 2.0
