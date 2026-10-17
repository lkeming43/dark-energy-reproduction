"""Common plotting interfaces for the Figure 2 and Figure 3 layouts."""


def plot_distance_panel(dataset, fit_result, confidence_bands, references, ax=None):
    """Plot the D_M/r_d reconstruction panel."""
    raise NotImplementedError


def plot_w_panel(reconstruction, confidence_bands, references, ax=None):
    """Plot the w(z) reconstruction panel."""
    raise NotImplementedError


def plot_figure_2_or_3(result, config):
    """Create the shared two-panel paper-style figure."""
    raise NotImplementedError

