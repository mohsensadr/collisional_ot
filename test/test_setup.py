import importlib.util
import pathlib
import os
import pytest
import setuptools

def test_setup_py_defines_expected_extension(tmp_path, monkeypatch):
    # Path to setup.py at project root
    setup_path = pathlib.Path(__file__).resolve().parent.parent / "setup.py"

    # Prepare a fake setup() to capture arguments
    captured_args = {}
    def fake_setup(**kwargs):
        captured_args.update(kwargs)

    # Monkeypatch setuptools.setup
    monkeypatch.setattr(setuptools, "setup", fake_setup)

    cwd = os.getcwd()
    try:
        # Change working directory to project root
        os.chdir(setup_path.parent)

        # Load and execute setup.py module
        spec = importlib.util.spec_from_file_location("setup_module", setup_path)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
    finally:
        os.chdir(cwd)

    # Check that 'ext_modules' is defined in setup() args
    assert "ext_modules" in captured_args
    ext_modules = captured_args["ext_modules"]
    names = [ext.name for ext in ext_modules]
    assert "collisional_ot.cython.collision_wrapper" in names

