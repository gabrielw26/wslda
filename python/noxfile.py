import nox

# import nox_poetry.patch
from nox.sessions import Session


@nox.session(python=["3.6", "3.7", "3.8", "3.9"], reuse_venv=False)  # "2.7", "3.5",
def test(session: Session) -> None:
    """Run the test suite."""
    print()
    session.run("poetry", "env", "use", session.virtualenv.interpreter, external=True)
    # session.install(".")
    # session.run("poetry", "shell", external=True)
    session.run("poetry", "install", external=True)
    session.run("poetry", "run", "pytest", external=True)
    session.run("poetry", "run", "coverage", "xml", "--fail-under=0", external=True)
