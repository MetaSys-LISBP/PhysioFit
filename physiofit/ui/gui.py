from ast import literal_eval
import json
import tkinter as tk
from tkinter import filedialog
from pathlib import Path
from copy import copy

import streamlit as st
import pandas as pd

from physiofit.base.io import IoHandler


class App:

    def __init__(self):

        self.select_menu = None
        self.io_handler = None
        self.defaults = {
            "vini": 0.2,
            "weight": {},
            "conc_met_bounds": [1e-06, 50],
            "flux_met_bounds": [-50, 50],
            "conc_biom_bounds": [1e-06, 50],
            "flux_biom_bounds": [-50, 50],
            "t_lag": 0,
            "deg": {},
            "iterations": 100
        }

    def start_app(self):

        st.title("Welcome to Physiofit 2.0")
        self.select_menu = st.selectbox(
            "Select task to execute",
            ["Calculate extracellular fluxes", "Calculate degradation constant"]
        )
        if self.select_menu == "Calculate extracellular fluxes":
            self.io_handler = IoHandler()
            self._build_flux_menu()
        else:
            st.header("Not yet implemented")

    def _build_flux_menu(self):

        self.data_file = st.file_uploader(
            "Upload data file or json configuration file",
            key="data_uploader",
            accept_multiple_files=False
        )

        if self.data_file:
            self._initialize_opt_menu()

    def _initialize_opt_menu(self):

        file_extension = self.data_file.name.split(".")[1]
        input_values = self.defaults
        if file_extension == "json":
            config = IoHandler.read_json_config(self.data_file)
            input_values.update(config)
        elif file_extension != "tsv":
            raise KeyError(f"Wrong input file format. Accepted formats are tsv for the data file or json for config "
                           f"files. Detected file:{self.data_file.name}")
        else:
            data = pd.read_csv(self.data_file, sep="\t")
            try:
                self.defaults["weight"].update({"X": 0.2})
                for col in data.columns[2:]:
                    self.defaults["weight"].update({col: 0.5})
            except Exception:
                raise

        submitted = self._initialize_opt_menu_widgets(input_values, file_extension)

        if submitted:
            if file_extension == "tsv":
                self.io_handler.data = pd.read_csv(self.data_file, sep="\t")
                self.io_handler.data = self.io_handler.data.sort_values("time", ignore_index=True)
                self.io_handler.names = self.io_handler.data.columns[1:].to_list()
                kwargs = self._build_fitter_kwargs()
                self.io_handler.initialize_fitter(kwargs)
            if file_extension == "json":
                st.session_state.submitted = True
                final_json = self._build_internal_json(input_values["path_to_data"])
                self.io_handler.launch_from_json(final_json)
            self.io_handler.fitter.optimize()
            if self.mc:
                self.io_handler.fitter.monte_carlo_analysis()
            self.io_handler.fitter.khi2_test()
            outputs = ["data", "plot", "pdf"]
            self.io_handler.local_out(*outputs)
            st.write(f"Run is finished. Check {self.io_handler.res_path} for the results.")

    def _initialize_opt_menu_widgets(self, input_values, file_extension):

        expand_basic_settings = st.expander("Basic settings", expanded=True)
        with expand_basic_settings:
            self.t_lag = st.checkbox(
                "Lag",
                key="lag_check"
            )
            self.deg_check = st.checkbox(
                "Degradation",
                key="deg_check"
            )
            enable_deg = False if self.deg_check else True
            self.deg = st.text_input(
                "Degradation constants",
                value={},
                help="Dictionary of the different degradation constants for each metabolite. "
                     "Format: {'met1': value, 'met2': value, ... 'metn': value}",
                disabled=enable_deg
            )
            self.mc = st.checkbox(
                "Sensitivity analysis (Monte Carlo)",
                help="Should Monte Carlo statistical analysis be performed"
            )
            enable_mc = False if self.mc else True
            self.iterations = st.number_input(
                "Number of iterations",
                value=input_values["iterations"],
                help="How many iterations should the Monte Carlo analysis perform",
                disabled=enable_mc
            )

            if file_extension == "tsv":

                # Set up tkinter for directory chooser
                root = tk.Tk()
                root.withdraw()

                # Make folder picker dialog appear on top of other windows
                root.wm_attributes('-topmost', 1)

                # Initialize folder picker button and add logic
                clicked = st.button("Select destination folder", key="clicker")
                if clicked:

                    # Initialize home path from directory selector and add to session state
                    st.session_state.home_path = Path(st.text_input(
                        "Selected output folder:", filedialog.askdirectory(master=root)
                    ))
                    if st.session_state.home_path == Path(".") or not st.session_state.home_path.exists():
                        raise RuntimeError("Please input a valid file path")
                    self.io_handler.home_path = copy(st.session_state.home_path)

                elif hasattr(st.session_state, "home_path"):

                    self.io_handler.home_path = Path(st.text_input(
                        "Selected output folder:", st.session_state.home_path
                    ))
                    if self.io_handler.home_path == Path(".") or not self.io_handler.home_path.exists():
                        raise RuntimeError("Please input a valid file path")

                    # Initialize the result export directory
                    self.io_handler.res_path = self.io_handler.home_path / (self.io_handler.home_path.name + "_res")
                    if not clicked:
                        if not self.io_handler.res_path.is_dir():
                            self.io_handler.res_path.mkdir()

        # Build the form for advanced parameters
        form = st.form("Parameter_form")
        with form:
            expand_run_params = st.expander("Advanced settings")
            with expand_run_params:
                self.vini = st.text_input(
                    "Initial flux values",
                    value=input_values["vini"],
                    key="vini",
                    help="Select an initial value for the parameters to estimate. Default=0.2"
                )
                self.weight = st.text_input(
                    "Weights for cost calculation & sensitivity analysis" "Weights",
                    value=input_values["weight"],
                    help="Input weights to apply on the different variables. If empty, default is 0.02 for biomass and"
                         " 0.05 for metabolites"
                )
                self.conc_met_bounds = st.text_input(
                    "Bounds on metabolite concentrations" ,
                    value=input_values["conc_met_bounds"],
                    help="Give the bounds for the initial concentrations of the metabolites (Mi0 value). "
                         "They will limit the range of possibilities during optimization. Defaults: [1e-06, 50]"
                )
                self.flux_met_bounds = st.text_input(
                    "Bounds on fluxes",
                    value=input_values["flux_met_bounds"],
                    help="Give the bounds for the metabolite fluxes (q value). They will limit the range of "
                         "possibilities during optimization. Defaults: [0.01, 50]"
                )
                self.conc_biom_bounds = st.text_input(
                    "Biomass initial concentration bounds",
                    value=input_values["conc_biom_bounds"],
                    help="Give the bounds for the initial concentrations of the metabolites (Mi0 value). "
                         "They will limit the range of possibilities during optimization. Defaults: [1e-06, 50]"
                )
                self.flux_biom_bounds = st.text_input(
                    "Bounds on growth rate",
                    value=input_values["flux_biom_bounds"],
                    help="Give the bounds for the metabolite fluxes (q value). They will limit the range of "
                         "possibilities during optimization. Defaults: [0.01, 50]"
                )
                self.debug_mode = st.checkbox(
                    "Debug Mode",
                    help="Should debug information be written to the log file"
                )
            submitted = st.form_submit_button("Run flux calculation")
        return submitted

    def _build_fitter_kwargs(self):

        kwargs = {
            "vini": float(self.vini),
            "weight": literal_eval(self.weight),
            "conc_met_bounds": tuple(literal_eval(self.conc_met_bounds)),
            "flux_met_bounds": tuple(literal_eval(self.flux_met_bounds)),
            "conc_biom_bounds": tuple(literal_eval(self.conc_biom_bounds)),
            "flux_biom_bounds": tuple(literal_eval(self.flux_biom_bounds)),
            "t_lag": self.t_lag,
            "deg": literal_eval(self.deg),
            "mc": self.mc,
            "iterations": self.iterations,
            "debug_mode": self.debug_mode,
        }
        return kwargs

    def _build_internal_json(self, path_to_data):

        final_json = json.dumps({
            "vini": float(self.vini),
            "weight": literal_eval(self.weight),
            "conc_met_bounds": literal_eval(self.conc_met_bounds),
            "flux_met_bounds": literal_eval(self.flux_met_bounds),
            "conc_biom_bounds": literal_eval(self.conc_biom_bounds),
            "flux_biom_bounds": literal_eval(self.flux_biom_bounds),
            "t_lag": self.t_lag,
            "deg": literal_eval(self.deg),
            "mc": self.mc,
            "iterations": self.iterations,
            "debug_mode": self.debug_mode,
            "path_to_data": path_to_data
        })
        return final_json


if __name__ == "__main__":
    my_app = App()
    my_app.start_app()
