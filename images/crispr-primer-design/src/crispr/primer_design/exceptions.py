"""Exceptions raised by primer_design."""


class PrimerDesignError(Exception):
    """Base class for all primer design errors."""

class PrimerDesignIOError(PrimerDesignError):
    """An IO error occurred during primer design."""

class PrimerDesignGenomeError(PrimerDesignError):
    """An error with the reference genome files."""

class PrimerDesignGuideError(PrimerDesignError):
    """An error processing the guides for primer design."""

class PrimerDesignNGSPrimerError(PrimerDesignError):
    """An error while generating NGS primers."""
