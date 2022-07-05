import os
import sys
import time
import numpy as np
import pandas as pd


__all__ = ['qsc_order_gene', 'qsc_distribution',
           'qsc_activation_ratios']


def _qsc_probabilities(ngenes, labels, counts, drop_zero=True):
    """
    Complete the spectrum of quantum states for a given distribution of states.
    Parameters
    ----------
    ngenes : int
        The number of genes in the Quantum GRN.
    labels : np.array
        The existing labels in the distribution of states.
    counts : np.array
        The counts for the labels in the distribution of states.
    drop_zeros : bool
        Flag to normalize the distributions. It set the `|0>_n` state to 0 and
        rescale the rest of the distribution.
    Returns
    -------
    probabilities: np.array
        The complete spectrum of quantum states for the QuantumGRN model.
    """
    length = 2 ** ngenes
    probabilities = np.zeros(length)

    nlabels, _ = labels.shape
    idx = np.zeros(nlabels, dtype=int)
    for i in range(ngenes):
        idx = 2 * idx + labels[:, i]

    probabilities[idx] = counts
    if drop_zero:
        probabilities[0] = 0.

    return probabilities / np.sum(probabilities)


def qsc_distribution(dataframe, threshold=0, drop_zero=True):
    """
    Compute the observed distribution for the binarization step for the dataset
    using the threshold values as described in the manuscript.
    Parameters
    ----------
    dataframe : pd.DataFrame
        The dataset as DataFrame object where the columns are genes and
        rows are cells.
    threshold : float
        The threshold value for the binarization.
        Default: 0
    drop_zero : bool
        Flag to normalize the distributions. It set the `|0>_n` state to 0 and
        rescale the rest of the distribution.
    Returns
    -------
    prob : np.array
        The complete spectrum of the observed distribution `p_obs` for a
        given dataset.
    Raises
    ------
    ValueError
        If dataframe is not a pd.DataFrame object.
    """
    if not isinstance(dataframe, pd.DataFrame):
        raise ValueError("The dataset parameter is not a pd.DataFrame "
                         "object")

    scdata = dataframe.to_numpy().T
    ngenes, _ = scdata.shape

    labels, counts = np.unique((scdata > threshold) + 0, axis=1, \
                               return_counts=True)
    labels = np.flip(labels, axis=0).T
    prob = _qsc_probabilities(ngenes, labels, counts, drop_zero)
    info_print("The observed probability `p_obs` is calculated")
    return prob


def qsc_order_gene(dataframe, threshold=0):
    """
    Order the dataset according to the expression level of each gene from
    largest to smallest.
    Parameters
    ----------
    dataframe : pd.DataFrame
        The dataset as DataFrame object where the columns are genes and
        rows are cells.
    threshold : float
        The threshold value for the binarization.
        Default: 0
    Returns
    -------
    dataframe : pd.DataFrame
        The ordered dataset as DataFrame object
    Raises
    ------
    ValueError
        If dataframe is not a pd.DataFrame object.
    """
    if not isinstance(dataframe, pd.DataFrame):
        raise ValueError("The dataset parameter is not a pd.DataFrame "
                         "object")
    mask = dataframe > threshold
    counts = mask.sum(axis=0).sort_values(ascending=False)
    ordered_genes = counts.index.to_list()
    info_print("The dataframe genes are ordered")
    return dataframe[ordered_genes]


def qsc_activation_ratios(dataframe, threshold=0):
    """
    Compute the activation ratios for each gene in the ordered dataset.
    Parameters
    ----------
    dataframe : pd.DataFrame
        The dataset as DataFrame object where the columns are genes and
        rows are cells.
    threshold : float
        The threshold value for the binarization.
        Default: 0
    Returns
    -------
    dataframe : pd.DataFrame
        The ordered dataset as DataFrame object
    Raises
    ------
    ValueError
        If dataframe is not a pd.DataFrame object.
    """
    if not isinstance(dataframe, pd.DataFrame):
        raise ValueError("The dataset parameter is not a pd.DataFrame "
                         "object")
    scdata = dataframe.to_numpy().T
    _, ncells = scdata.shape
    info_print("Activation ratios are computed")
    return np.sum((scdata > threshold) + 0, axis=1) / ncells


def _qsc_labels(ngenes):
    """
    Retrieves a complete set of basis state in the QuantumGRN given the number
    of genes. It is used in histogram visualization.
    Parameters
    ----------
    ngenes : int
        The number of genes for GRN modelling.
    Returns
    -------
    labels : list
        The complete set of basis state in the GRN modelling.
    Raises
    ------
    TypeError
        If the ngene parameters is not an integer.
    """
    if not isinstance(ngenes, int):
        raise TypeError("The number of genes must be an integer value")
    labels = []

    for i in range(2**ngenes):
        label = np.binary_repr(i, width=ngenes)
        labels.append(label)

    return labels


def info_print(msg, level="I"):
    date = time.strftime("%Y-%m-%d %H:%M:%S")
    sys.stdout.write("{date} | {level} | {msg}\n"
                     .format(date=date, msg=msg, level=level))


def _print_msg(message, line_break=True):
    if line_break:
        sys.stdout.write(message + "\n")
    else:
        sys.stdout.write(message)
    sys.stdout.flush()


class Progbar:
    """
    target: Total number of steps expected, None if unknown.
    width: Progress bar width on screen
    interval: Minimum visual progress update interval (in seconds)
    unit_name: Display name for step counts
    re-use from tf.keras.utils.progbar
    """

    def __init__(self, target, width=50, interval=1,
                 unit_name="step"):
        self.target = target
        self.width = width
        self.interval = interval
        self.unit_name = unit_name

        self._dynamic_display = ((hasattr(sys.stdout, 'isatty') and
                                  sys.stdout.isatty()) or
                                 'ipykernel' in sys.modules or
                                 'posix' in sys.modules or
                                 'PYCHARM_HOSTED' in os.environ)
        self._total_width = 0
        self._seen_so_far = 0
        # We use a dict + list to avoid garbage collection
        # issues found in OrderedDict
        self._start = time.time()
        self._last_update = 0
        self._time_at_epoch_start = self._start
        self._time_at_epoch_end = None
        self._time_after_first_step = None


    def update(self, current, finalize=None):
        """
        Update the progress bar.
        """
        if finalize is None:
            if self.target is None:
                finalize = False
            else:
                finalize = current >= self.target

        self._seen_so_far = current

        message = ''
        now = time.time()
        info = ' - %.0fs' % (now - self._start)
        if current == self.target:
            self._time_at_epoch_end = now
        if now - self._last_update < self.interval and not finalize:
            return

        prev_total_width = self._total_width
        if self._dynamic_display:
            message += '\b' * prev_total_width
            message += '\r'
        else:
            message += '\n'

        numdigits = int(np.log10(self.target)) + 1
        bar = ('%' + str(numdigits) + 'd/%d [') % (current, self.target)
        prog = float(current) / self.target
        prog_width = int(self.width * prog)
        if prog_width > 0:
            bar += ('=' * (prog_width - 1))
            if current < self.target:
                bar += '>'
            else:
                bar += '='
        bar += ('.' * (self.width - prog_width))
        bar += ']'

        self._total_width = len(bar)
        message += bar

        time_per_unit = self._estimate_step_duration(current, now)

        if self.target is None or finalize:
            info += self._format_time(time_per_unit, self.unit_name)
        else:
            eta = time_per_unit * (self.target - current)
            eta_format = '%d:%02d' % (eta // 60, eta % 60)
            info = ' - ETA: %s' % eta_format


        self._total_width += len(info)
        if prev_total_width > self._total_width:
            info += (' ' * (prev_total_width - self._total_width))

        if finalize:
            info += '\n'

        message += info
        _print_msg(message, line_break=False)
        message = ''

    def _estimate_step_duration(self, current, now):
        if self._time_after_first_step is not None and current > 1:
            time_per_unit = (now - self._time_after_first_step) / (current - 1)
        else:
            time_per_unit = (now - self._start) / current

        if current == 1:
            self._time_after_first_step = now
        return time_per_unit

    def _format_time(self, time_per_unit, unit_name):
        """format a given duration to display to the user."""
        formatted = ''
        if time_per_unit >= 1 or time_per_unit == 0:
            formatted += ' %.0fs/%s' % (time_per_unit, unit_name)
        elif time_per_unit >= 1e-3:
            formatted += ' %.0fms/%s' % (time_per_unit * 1e3, unit_name)
        else:
            formatted += ' %.0fus/%s' % (time_per_unit * 1e6, unit_name)
        return formatted


