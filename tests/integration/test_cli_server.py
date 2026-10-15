"""Integration tests for beg server CLI subcommand."""

import subprocess
import sys


def test_cli_server_help() -> None:
    """Test 'beg server --help' displays subcommand description and arguments."""
    result = subprocess.run(
        [sys.executable, "beg.py", "server", "--help"],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0
    assert "Launch interactive web browser module" in result.stdout
    assert "--host" in result.stdout
    assert "--port" in result.stdout
    assert "--reloadable" in result.stdout
