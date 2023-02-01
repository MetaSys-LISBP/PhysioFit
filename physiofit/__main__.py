from subprocess import run
import sys
from pathlib import Path

from threading import Thread
import requests

import physiofit
import physiofit.ui.cli


def get_last_version():
    """Get last Physiofit version."""
    try:
        pf_path = Path(physiofit.__file__).parent
        # Get the version from pypi
        response = requests.get(f'https://pypi.org/pypi/physiofit/json')
        latest_version = response.json()['info']['version']
        with open(str(Path(pf_path, "last_version.txt")), "w") as f:
            f.write(latest_version)
    except Exception:
        pass

def main():
    """The main routine"""

    if len(sys.argv) > 1:
        physiofit.ui.cli.main()
    else:
        thread = Thread(target=get_last_version)
        thread.start()
        path_to_app = Path(physiofit.__file__).parent
        path_to_app = path_to_app / "ui/gui.py"
        run(["streamlit", "run", str(path_to_app)])

if __name__ == "__main__":
    sys.exit(main())