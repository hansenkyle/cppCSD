# Copyright (c) 2026, Kyle Hansen (khansen3@ncsu.edu)
#
# Funded by CARRE (https://carre-psaapiv.org/)
#
# Licensed under BSD 3-Clause License; Redistribution and use in source and binary forms, with
# or without modification are permitted provided that the terms of the license are met.

from .results import (
    Convergence,
    Input,
    MaterialData,
    Run,
    cell_average,
    multigroup,
    read,
    spectrum,
)

__all__ = [
    "Convergence",
    "Input",
    "MaterialData",
    "Run",
    "cell_average",
    "multigroup",
    "read",
    "spectrum",
]

# The C++ extension (input deck builder + solvers). Optional, so the reader above still works
# from a plain checkout without a C++ build (e.g. PYTHONPATH=python). Only a missing module is
# tolerated; one that's built but fails to load still raises.
try:
    from ._core import (
        Angle,
        Energy,
        InputDeck,
        Material,
        Mesh,
        Method,
        SecondMoment,
        SourceIteration,
        configure_logging,
    )
except ModuleNotFoundError as e:
    if e.name != f"{__name__}._core":
        raise
else:
    __all__ += [
        "Angle",
        "Energy",
        "InputDeck",
        "Material",
        "Mesh",
        "Method",
        "SecondMoment",
        "SourceIteration",
        "configure_logging",
    ]
