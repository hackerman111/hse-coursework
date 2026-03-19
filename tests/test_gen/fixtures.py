"""
Общие фикстуры для тестов генераторов.
"""

from sage.all import QQ


class GeneratorTestMixin:
    def setUp(self):
        self.base_field = QQ
