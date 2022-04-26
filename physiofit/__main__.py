from subprocess import run
import sys
import pathlib

import physiofit
import physiofit.ui.cli

def main():
    """The main routine"""

    if len(sys.argv) > 1:
        physiofit.ui.cli.main()
    else:
        path_to_app = pathlib.Path(physiofit.__path__[0])
        path_to_app = path_to_app / "ui/gui.py"
        run(["streamlit", "run", str(path_to_app)])

if __name__ == "__main__":
    sys.exit(main())