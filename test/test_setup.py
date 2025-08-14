import importlib.util
import pathlib
import os
import pytest

def test_setup_py_defines_expected_extension(tmp_path, monkeypatch):
    # Path to setup.py
    setup_path = pathlib.Path(__file__).parent.parent / "src" / "setup.py"

    # Prepare a fake setup() to capture arguments
    captured_args = {}
    def fake_setup(**kwargs):
        captured_args.update(kwargs)

    # Monkeypatch distutils.core.setup
    monkeypatch.setattr("distutils.core.setup", fake_setup)

    # Save current working directory
    cwd = os.getcwd()
    try:
        # Change working directory to src/ so cythonize can find .pyx file
        os.chdir(setup_path.parent)

        # Load and execute setup.py module
        spec = importlib.util.spec_from_file_location("setup_module", setup_path)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
    finally:
        # Restore original working directory
        os.chdir(cwd)

    # Check that 'ext_modules' is defined in setup() args
    assert "ext_modules" in captured_args
    ext_modules = captured_args["ext_modules"]
    # Ensure at least one extension with expected name
    names = [ext.name for ext in ext_modules]
    assert "collision_wrapper" in names

