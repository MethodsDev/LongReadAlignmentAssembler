def pytest_configure(config):
    config.addinivalue_line(
        "markers",
        "serial: measures real CPU/timing; run outside the parallel (xdist) pass",
    )
