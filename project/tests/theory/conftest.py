"""Markers for the theory tests."""


def pytest_configure(config):
    config.addinivalue_line("markers", "slow: test that takes longer than ~10 s (deselect with -m 'not slow')")
