"""Structured logger factory for PDAC research modules."""

import logging
import sys


def get_logger(name: str) -> logging.Logger:
    """Return a named logger with structured-friendly formatting.

    Args:
        name: Module name, typically ``__name__``.

    Returns:
        Configured :class:`logging.Logger` instance.
    """
    logger = logging.getLogger(name)
    if not logger.handlers:
        handler = logging.StreamHandler(sys.stderr)
        handler.setFormatter(
            logging.Formatter("%(asctime)s %(levelname)s %(name)s %(message)s")
        )
        logger.addHandler(handler)
    logger.setLevel(logging.INFO)
    return logger
