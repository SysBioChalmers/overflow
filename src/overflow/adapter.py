"""Model adapter and model loading.

The adapter carries the organism-specific parameters that geckopy reads
(biomass reaction, carbon source, protein content, ...). Its defaults
come from ``model_adapter.toml`` at the repository root; per-condition
values are layered on top.
"""
from __future__ import annotations

from typing import TYPE_CHECKING, Any, Optional

from geckopy import ModelAdapter, load_conventional_gem, load_ec_model

from overflow.config import CONV_GEM, EC_MODEL, ROOT

if TYPE_CHECKING:  # pragma: no cover - typing only
    import cobra
    from geckopy import EcModel

    from overflow.config import Condition


def build_adapter(condition: Optional["Condition"] = None, **overrides: Any) -> ModelAdapter:
    """Load the repository adapter, optionally specialised to a condition.

    ``condition`` supplies the measured total protein content and the
    dilution rate as the reference growth rate. Explicit keyword
    ``overrides`` win over both.
    """
    adapter = ModelAdapter.from_folder(ROOT)
    updates: dict[str, Any] = {}
    if condition is not None:
        updates["p_tot"] = condition.p_tot
        updates["gr_exp"] = condition.d_rate
    updates.update(overrides)
    if updates:
        adapter.params = adapter.params.model_copy(update=updates)
    return adapter


def load_model(adapter: Optional[ModelAdapter] = None, **overrides: Any) -> "EcModel":
    """Load the GECKO 4 ecModel distributed with this repository."""
    if adapter is None:
        adapter = build_adapter(**overrides)
    return load_ec_model(str(EC_MODEL), adapter=adapter)


def load_gem(adapter: Optional[ModelAdapter] = None, **overrides: Any) -> "cobra.Model":
    """Load the conventional yeast-GEM the ecModel was built from.

    The returned model carries the adapter, which the geckopy flux-data
    helpers read even when acting on a plain cobra model.
    """
    if adapter is None:
        adapter = build_adapter(**overrides)
    model = load_conventional_gem(adapter)
    model.adapter = adapter
    return model
