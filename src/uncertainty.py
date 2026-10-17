"""Seed stability and Monte Carlo uncertainty interfaces."""


def run_seed_ensemble(dataset, pysr_config, seeds):
    """Run the shared PySR stability experiment over multiple seeds."""
    raise NotImplementedError


def summarize_seed_ensemble(results):
    """Summarize selected structures, parameters, losses, and diagnostics."""
    raise NotImplementedError


def run_monte_carlo(dataset, fit_result, mc_config):
    """Propagate observational uncertainty through the reconstruction."""
    raise NotImplementedError


def summarize_confidence_bands(samples, levels=(1, 2, 3)):
    """Return central curves and confidence bands."""
    raise NotImplementedError

