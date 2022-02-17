from subprocess import run
import sys
import pathlib

import physiofit

def main():
    """The main routine"""

    path_to_app = pathlib.Path(physiofit.__path__[0])
    path_to_app = path_to_app / "gui/gui.py"
    run(["streamlit", "run", str(path_to_app)])

if __name__ == "__main__":
    sys.exit(main())