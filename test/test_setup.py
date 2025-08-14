import importlib.util
import pathlib
import types


def test_setup_py_defines_expected_extension(tmp_path, monkeypatch):
    setup_path = pathlib.Path(__file__).parent.parent / "src" / "setup.py"

    # Prepare a fake setup() to capture arguments
    captured_args = {}
    def fake_setup(**kwargs):
        captured_args.update(kwargs)

    monkeypatch.setattr("distutils.core.setup", fake_setup)

    # Load the setup.py module (will run top-level code)
    spec = importlib.util.spec_from_file_location("setup_module", setup_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    # Check that our fake setup() was called with expected data
    assert "name" in captured_args
    assert captured_args["name"] == "OT_Collision"
    assert "ext_modules" in captured_args
    assert len(captured_args["ext_modules"]) >= 1

    # Check first extension properties
    ext = captured_args["ext_modules"][0]
    assert "collision_wrapper" in ext.name
    assert any("collision.c" in src for src in ext.sources)
    assert any("collision_wrapper.pyx" in src for src in ext.sources)

    # Check macros and compile args
    macros = dict(ext.define_macros)
    assert macros.get("NPY_NO_DEPRECATED_API") == "NPY_1_7_API_VERSION"
    assert "-O3" in ext.extra_compile_args

