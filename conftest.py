# conftest.py

import pytest

def pytest_collection_modifyitems(items):
    for item in items:
        if item.parent.name == "test_annotate.py":
            item.add_marker(pytest.mark.conda_env("py37_pysam"))