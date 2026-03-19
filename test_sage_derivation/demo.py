from sage.all import PolynomialRing, QQ


def show_basic_derivation_demo():
    ring = PolynomialRing(QQ, "x,y")
    x, y = ring.gens()

    try:
        # Для текущей версии Sage список образов работает стабильнее словаря.
        derivation = ring.derivation([y, 2 * x])
        print(f"Type of D: {type(derivation)}")
        print(f"D(x) = {derivation(x)}")
        print(f"D(x^2 + y) = {derivation(x**2 + y)}")

        derivation_dx = ring.derivation([1, 0])
        print(f"(D+D2)(x) = {(derivation + derivation_dx)(x)}")

        try:
            if hasattr(derivation, "bracket"):
                print("Has bracket method")
            if hasattr(derivation, "commutator"):
                print("Has commutator method")
            print(f"Dir of D: {dir(derivation)}")
        except Exception as exc:
            print(f"Bracket check failed: {exc}")

        try:
            if hasattr(derivation, "leading_term"):
                print(f"Leading term: {derivation.leading_term()}")
            else:
                print("No leading_term method")
        except Exception as exc:
            print(f"Leading term check failed: {exc}")
    except Exception as exc:
        print(f"Creation failed: {exc}")
