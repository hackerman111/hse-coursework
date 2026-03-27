"""
I/O слой: разбор и форматирование генераторов.
"""

from .parse import (
    parse_generator_spec,
    spec_from_string,
    spec_to_string,
    describe_generator,
    describe_random_generator,
    generator_term_count,
    generator_max_degree,
)

__all__ = [
    "parse_generator_spec",
    "spec_from_string",
    "spec_to_string",
    "describe_generator",
    "describe_random_generator",
    "generator_term_count",
    "generator_max_degree",
]
