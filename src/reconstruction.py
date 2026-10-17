"""Distance derivatives, cosmological reference curves, and w(z) interfaces."""


def evaluate_distance(expression, redshift):
    """Evaluate a fitted distance expression on a redshift grid."""
    raise NotImplementedError


def evaluate_distance_derivatives(expression, redshift):
    """Evaluate the first and second redshift derivatives."""
    raise NotImplementedError


def reconstruct_w(expression, redshift, cosmology):
    """Reconstruct w(z) from D_M/r_d and its derivatives."""
    raise NotImplementedError


def evaluate_reference_curves(redshift, config):
    """Evaluate the configured Lambda-CDM and CPL comparison curves."""
    raise NotImplementedError


def diagnose_reconstruction(expression, redshift, cosmology):
    """Check denominator crossings and other invalid reconstruction regions."""
    raise NotImplementedError

