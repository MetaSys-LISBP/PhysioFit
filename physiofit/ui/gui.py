from ast import literal_eval
import json
import tkinter as tk
from tkinter import filedialog
from pathlib import Path
from copy import copy

import streamlit as st
import pandas as pd

import physiofit
from physiofit.base.io import IoHandler, DEFAULTS


class App:
    """
    Physiofit Graphical User Interface
    """

    def __init__(self):

        self.defaults = copy(DEFAULTS)
        self.select_menu = None
        self.io_handler = None
        self.update_info = None

    def start_app(self):
        """Launch the application"""

        st.set_page_config(page_title=f"PhysioFit (v{physiofit.__version__})")
        st.title(f"Welcome to PhysioFit (v{physiofit.__version__})")
        self.update_info = st.empty()
        self.check_uptodate()
        self.select_menu = st.selectbox(
            "Select a task to execute",
            ["Calculate extracellular fluxes",
             "Calculate degradation constant"]
        )
        if self.select_menu == "Calculate extracellular fluxes":
            self.io_handler = IoHandler()
            self._build_flux_menu()
        else:
            st.header("Implementation in progress...")

    def check_uptodate(self):
        """Compare installed and most recent Physiofit versions."""
        try:
            pf_path = Path(physiofit.__file__).parent
            with open(str(Path(pf_path, "last_version.txt")), "r") as f:
                lastversion = f.read()
            if lastversion != physiofit.__version__:
                # change the next line to streamlit
                self.update_info = st.info(
                    f'New version available ({lastversion}). '
                    f'You can update PhysioFit with: "pip install --upgrade physiofit". '
                    f'Check the documentation for more information.')
        except:
            pass

    def _build_flux_menu(self):
        """Build the starting menu with the data upload button"""

        self.data_file = st.file_uploader(
            "Load a data file or a json configuration file",
            key="data_uploader",
            accept_multiple_files=False
        )

        if self.data_file:
            self._initialize_opt_menu()

    def _initialize_opt_menu(self):
        """
        Initialize all the options. If a json file is given as input, options
        are parsed from it. If the input is tsv, the defaults are used instead.
        """

        file_extension = self.data_file.name.split(".")[1]
        input_values = self.defaults
        if file_extension == "json":
            config = IoHandler.read_json_config(self.data_file)
            input_values.update(config)
        elif file_extension not in ["tsv", "txt"]:
            raise KeyError(
                f"Wrong input file format. Accepted formats are tsv for data "
                f"files or json for configuration files. Detected file: "
                f"{self.data_file.name}")
        else:
            data = pd.read_csv(self.data_file, sep="\t")
            self.io_handler.data = data
            self.io_handler.data = self.io_handler.data.sort_values(
                "time", ignore_index=True
            )
            self.io_handler.get_models()
            try:
                self.defaults["sd"].update({"X": 0.2})
                for col in data.columns[2:]:
                    self.defaults["sd"].update({col: 0.5})
            except Exception:
                raise

        submitted = self._initialize_opt_menu_widgets(
            input_values, file_extension
        )

        if submitted:
            if file_extension in ["tsv", "txt"]:
                self.io_handler.names = self.io_handler.data.columns[
                                        1:].to_list()
                kwargs = self._build_fitter_kwargs()
                self.io_handler.initialize_fitter(kwargs)
            if file_extension == "json":
                st.session_state.submitted = True
                final_json = self._build_internal_json(
                    input_values["path_to_data"]
                )
                self.io_handler.launch_from_json(final_json)
            self.io_handler.fitter.optimize()
            if self.mc:
                self.io_handler.fitter.monte_carlo_analysis()
            self.io_handler.fitter.khi2_test()
            outputs = ["data", "plot", "pdf"]
            self.io_handler.local_out(*outputs)
            st.success(
                f"Run is finished. Check {self.io_handler.res_path} for "
                f"the results."
            )

    def _initialize_opt_menu_widgets(self, input_values, file_extension):

        model_options = [
                model.model_name for model in self.io_handler.models
            ]
        model_options.insert(0, "--")
        model_name = st.selectbox(
            label="Model",
            options=model_options,
            key = "model_selector",
            help="Select the model to use for flux calculation"
        )

        if model_name != "--":
            for model in self.io_handler.models:
                if model.model_name == model_name:
                    self.model = model
                    self.model.get_params()
                    break
            expand_basic_settings = st.expander("Basic settings", expanded=True)
            with expand_basic_settings:
                if self.model.fixed_parameters is not None:
                    for key, val in self.model.fixed_parameters.items():
                        self.model.fixed_parameters[key] = st.text_input(
                            label = key,
                            value = val,
                            help = f"Input value for fixed parameter {key}"
                        )
                self.mc = st.checkbox(
                    "Sensitivity analysis (Monte Carlo)",
                    value=True,
                    help="Determine the precision on estimated fluxes by Monte "
                         "Carlo sensitivity analysis."
                )
                enable_mc = False if self.mc else True
                self.iterations = st.number_input(
                    "Number of iterations",
                    value=input_values["iterations"],
                    help="Number of iterations for the Monte Carlo analysis.",
                    disabled=enable_mc
                )

                if file_extension in ["tsv", "txt"]:

                    # Set up tkinter for directory chooser
                    root = tk.Tk()
                    root.withdraw()

                    # Make folder picker dialog appear on top of other windows
                    root.wm_attributes('-topmost', 1)

                    # Initialize folder picker button and add logic
                    clicked = st.button(
                        "Select output data directory", key="clicker"
                    )
                    if clicked:

                        # Initialize home path from directory selector and add
                        # to session state
                        st.session_state.home_path = Path(st.text_input(
                            "Selected output data directory:",
                            filedialog.askdirectory(master=root)
                        ))
                        if st.session_state.home_path == Path(".") \
                                or not st.session_state.home_path.exists():
                            raise RuntimeError("Please provide a valid path")
                        self.io_handler.home_path = copy(
                            st.session_state.home_path)

                    elif hasattr(st.session_state, "home_path"):

                        self.io_handler.home_path = Path(st.text_input(
                            "Selected output data directory:",
                            st.session_state.home_path
                        ))
                        if self.io_handler.home_path == Path(".") \
                                or not self.io_handler.home_path.exists():
                            raise RuntimeError("Please provide a valid path")

                        # Initialize the result export directory
                        self.io_handler.res_path = self.io_handler.home_path / (
                                    self.io_handler.home_path.name + "_res")
                        if not clicked:
                            if not self.io_handler.res_path.is_dir():
                                self.io_handler.res_path.mkdir()

            # Build the form for advanced parameters
            form = st.form("Parameter_form")
            with form:
                expand_run_params = st.expander("Advanced settings")
                with expand_run_params:
                    for key, val in self.model.initial_values.items():
                        self.model.initial_values[key] = st.text_input(
                            label=f"{key} initial value",
                            value=val,
                            help=f"Input initial value for {key}"
                        )
                        self.model.bounds[key] = st.text_input(
                            label=f"{key} bounds",
                            value= self.model.bounds[key],
                            help=f"Input bounds for {key} as a tuple (min, max)"
                        )
                    self.sd = st.text_input(
                        "Standard deviation on measurements",
                        value=input_values["sd"],
                        help="Standard deviation on the measurements. If empty, default is 0.2 for biomass and"
                             " 0.5 for metabolites"
                    )
                    self.debug_mode = st.checkbox(
                        "Verbose logs",
                        help="Useful in case of trouble. Join it to the issue on github."
                    )
                submitted = st.form_submit_button("Run flux calculation")
            return submitted

    def _build_fitter_kwargs(self):

        kwargs = {
            "vini": float(self.vini),
            "sd": literal_eval(self.sd),
            "model": self.model,
            "mc": self.mc,
            "iterations": self.iterations,
            "debug_mode": self.debug_mode,
        }
        return kwargs

    def _build_internal_json(self, path_to_data):

        final_json = json.dumps({
            "vini": float(self.vini),
            "sd": literal_eval(self.sd),
            "model": self.model.model_name,
            "mc": self.mc,
            "iterations": self.iterations,
            "debug_mode": self.debug_mode,
            "path_to_data": path_to_data
        })
        return final_json


if __name__ == "__main__":
    physiofit_app = App()
    physiofit_app.start_app()
