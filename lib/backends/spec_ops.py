"""
Backend-shaped operations over DerivationSpec.

These helpers are intentionally separate from `lib.core.spec` so that the core
module remains a DTO layer plus structural helpers.
"""

from __future__ import annotations

from lib.core.spec import DerivationSpec, Scalar


def _poly_diff(poly: dict, i: int) -> dict:
    """Дифференцировать полином по переменной i."""
    result: dict[tuple[int, ...], Scalar] = {}
    for alpha, coeff in poly.items():
        if alpha[i] > 0:
            new_alpha = tuple(a - (1 if k == i else 0) for k, a in enumerate(alpha))
            new_coeff = coeff * alpha[i]  # type: ignore[operator]
            prev = result.get(new_alpha, 0)
            result[new_alpha] = prev + new_coeff  # type: ignore[operator]
    return {alpha: coeff for alpha, coeff in result.items() if abs(float(coeff)) > 1e-15}


def _poly_mul(f: dict, g: dict) -> dict:
    """Перемножить два полинома."""
    result: dict[tuple[int, ...], Scalar] = {}
    for alpha, coeff_f in f.items():
        for beta, coeff_g in g.items():
            gamma = tuple(a + b for a, b in zip(alpha, beta))
            prev = result.get(gamma, 0)
            result[gamma] = prev + coeff_f * coeff_g  # type: ignore[operator]
    return {alpha: coeff for alpha, coeff in result.items() if abs(float(coeff)) > 1e-15}


def _spec_apply(spec: DerivationSpec, poly: dict) -> dict:
    """Применить деривацию к полиному: D(poly) = sum_i f_i * ∂(poly)/∂x_i."""
    result: dict[tuple[int, ...], Scalar] = {}
    for axis, component in spec.terms.items():
        differentiated = _poly_diff(poly, axis)
        if not differentiated:
            continue
        product = _poly_mul(component, differentiated)
        for alpha, coeff in product.items():
            prev = result.get(alpha, 0)
            result[alpha] = prev + coeff  # type: ignore[operator]
    return {alpha: coeff for alpha, coeff in result.items() if abs(float(coeff)) > 1e-15}


def lie_bracket(
    spec1: DerivationSpec,
    spec2: DerivationSpec,
    D_max: int | None = None,
) -> DerivationSpec:
    """
    Вычислить скобку Ли [spec1, spec2].

    На вход принимает два `DerivationSpec` и необязательное усечение `D_max`.
    На выходе возвращает их скобку Ли как новый `DerivationSpec`.

    >>> from lib.core.spec import monomial_spec, partial_spec, scale_spec
    >>> lie_bracket(monomial_spec(2, 0, (1, 0)), partial_spec(2, 0)) == scale_spec(partial_spec(2, 0), -1)
    True
    """
    if spec1.n != spec2.n:
        raise ValueError(
            f"lie_bracket требует одинакового n: {spec1.n} != {spec2.n}"
        )

    terms: dict[int, dict[tuple[int, ...], Scalar]] = {}

    for axis in range(spec1.n):
        second_component = spec2.terms.get(axis, {})
        first_component = spec1.terms.get(axis, {})

        contrib: dict[tuple[int, ...], Scalar] = {}
        for alpha, coeff in _spec_apply(spec1, second_component).items():
            prev = contrib.get(alpha, 0)
            contrib[alpha] = prev + coeff  # type: ignore[operator]
        for alpha, coeff in _spec_apply(spec2, first_component).items():
            prev = contrib.get(alpha, 0)
            contrib[alpha] = prev - coeff  # type: ignore[operator]

        clean = {
            alpha: coeff
            for alpha, coeff in contrib.items()
            if abs(float(coeff)) > 1e-15
        }
        if D_max is not None:
            clean = {
                alpha: coeff
                for alpha, coeff in clean.items()
                if sum(alpha) - 1 <= D_max
            }

        if clean:
            terms[axis] = clean

    return DerivationSpec(n=spec1.n, terms=terms)


def ad(spec: DerivationSpec):
    """
    Присоединённое действие: ad(D)(E) = [D, E].

    На вход принимает фиксированный `DerivationSpec` `D`.
    На выходе возвращает функцию, которая для `E` вычисляет `[D, E]`.

    >>> from lib.core.spec import monomial_spec, partial_spec, scale_spec
    >>> op = ad(monomial_spec(2, 0, (1, 0)))
    >>> op(partial_spec(2, 0)) == scale_spec(partial_spec(2, 0), -1)
    True
    """

    def _ad(other: DerivationSpec, D_max: int | None = None) -> DerivationSpec:
        return lie_bracket(spec, other, D_max=D_max)

    return _ad
