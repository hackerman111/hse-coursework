import traceback

from .scenario import test_quotient_circle


def main():
    try:
        test_quotient_circle()
    except ImportError:
        print("SageMath недоступен, пропуск выполнения теста.")
    except Exception as exc:
        print(f"Тест завершился с ошибкой: {exc}")
        traceback.print_exc()


if __name__ == "__main__":
    main()
