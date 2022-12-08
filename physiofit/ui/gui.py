from ast import literal_eval
import json
import tkinter as tk
from tkinter import filedialog
from pathlib import Path
from copy import copy

import streamlit as st
import pandas as pd

import physiofit
from physiofit.base.io import IoHandler
from physiofit.models.base_model import StandardDevs


class App:
    """
    Physiofit Graphical User Interface
    """

    def __init__(self):

        self.defaults = {
            "iterations" : 100,
            "sd" : StandardDevs()
        }
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
                    f'You can update PhysioFit with: "pip install --upgrade '
                    f'physiofit". Check the documentation for more '
                    f'information.'
                )
        except Exception:
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

        # Check file extension and initialize defaults
        file_extension = self.data_file.name.split(".")[1]
        input_values = self.defaults
        if file_extension == "json":
            # Get parameters from json file
            config = IoHandler.read_json_config(self.data_file)
            input_values.update(config)
        elif file_extension not in ["tsv", "txt"]:
            raise KeyError(
                f"Wrong input file format. Accepted formats are tsv for data "
                f"files or json for configuration files. Detected file: "
                f"{self.data_file.name}")
        else:
            # Load data into io_handler
            data = pd.read_csv(self.data_file, sep="\t")
            self.io_handler.data = data
            self.io_handler.data = self.io_handler.data.sort_values(
                "time", ignore_index=True
            )

            # Initialize default SDs
            self.defaults["sd"].update({
                "X" : 0.5
            })
            self.defaults["sd"].update({
                col: 0.2 for col in self.io_handler.data.columns[2:]
            })
            # Initialize the models
            self.io_handler.get_models()

        # Build menu
        submitted = self._initialize_opt_menu_widgets(
            file_extension
        )

        if submitted:
            session_data = self._get_data_from_session_state()
            # Initialize the fitter object
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
            # Do the work and export results
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

    def _initialize_opt_menu_widgets(self, file_extension):

        # Get model names and give as options to user
        model_options = [
            model.model_name for model in self.io_handler.models
        ]
        model_options.insert(0, "--")
        model_name = st.selectbox(
            label="Model",
            options=model_options,
            key="model_selector",
            help="Select the model to use for flux calculation"
        )

        if model_name != "--":
            # Initialize selected model
            for model in self.io_handler.models:
                if model.model_name == model_name:
                    self.model = model
                    self.model.get_params()
                    break
            expand_basic_settings = st.expander("Basic settings",
                                                expanded=True)
            with expand_basic_settings:
                # Select model parameters
                self.mc = st.checkbox(
                    "Sensitivity analysis (Monte Carlo)",
                    value=True,
                    help="Determine the precision on estimated fluxes by "
                         "Monte Carlo sensitivity analysis."
                )
                enable_mc = False if self.mc else True
                self.iterations = st.number_input(
                    "Number of iterations",
                    value=self.defaults["iterations"],
                    help="Number of iterations for the Monte Carlo analysis.",
                    disabled=enable_mc
                )
                self.debug_mode = st.checkbox(
                    "Verbose logs",
                    help="Useful in case of trouble. Join it to the "
                         "issue on github."
                )

                if file_extension in ["tsv", "txt"]:

                    self._output_directory_selector()

            # Build the form for advanced parameters
            form = st.form("Parameter_form")
            with form:
                expand_run_params = st.expander("Parameters")
                with expand_run_params:
                    st.subheader("Parameters to estimate")
                    col1, col2, col3, col4 = st.columns(4)
                    with col1:
                        st.write("Parameter Name")
                        for key in self.model.parameters_to_estimate:
                            st.text_input(
                                label="label", # Unused
                                label_visibility="collapsed",
                                value=key,
                                key=f"Parameter_name_{key}",
                                disabled=True
                            )
                    with col2:
                        st.write("Parameter Value")
                        for key, value in self.model.parameters_to_estimate.items():
                            st.text_input(
                                label="label", # Unused
                                label_visibility = "collapsed",
                                value=value,
                                key=f"Parameter_value_{key}"
                            )
                    with col3:
                        st.write("Lower Bound")
                        for key, bound in self.model.bounds.items():
                            st.text_input(
                                label="label",  # Unused
                                label_visibility="collapsed",
                                value=bound[0],
                                key=f"Parameter_lower_{key}"
                            )
                    with col4:
                        st.write("Upper Bound")
                        for key, bound in self.model.bounds.items():
                            st.text_input(
                                label="label",  # Unused
                                label_visibility="collapsed",
                                value=bound[1],
                                key=f"Parameter_upper_{key}"
                            )

                    if self.model.fixed_parameters is not None:
                        for param in self.model.fixed_parameters.keys():
                            st.subheader(f"Fixed parameters: {param}")
                            col1, col2 = st.columns(2)
                            with col1:
                                st.write("Parameter Name")
                                for key in self.model.fixed_parameters[param].keys():
                                    st.text_input(
                                        label="label", # Unused
                                        label_visibility="collapsed",
                                        value=key,
                                        key=f"Fixed_{param}_{key}",
                                        disabled=True
                                    )
                            with col2:
                                st.write("Parameter Value")
                                for key, value in self.model.fixed_parameters[param].items():
                                    st.text_input(
                                        label="label", # Unused
                                        label_visibility="collapsed",
                                        value=value,
                                        key=f"Fixed_{param}_value_{key}"
                                    )

                expand_sds = st.expander("Standard Deviations")
                with expand_sds:
                    col1, col2 = st.columns(2)
                    with col1:
                        for key in self.defaults["sd"].keys():
                            st.text_input(
                                label="label",  # Unused
                                label_visibility="collapsed",
                                value=key,
                                key=f"{key}_sd",
                                disabled=True
                            )
                    with col2:
                        for key, value in self.defaults["sd"].items():
                            st.text_input(
                                label="label",  # Unused
                                label_visibility="collapsed",
                                value=value,
                                key=f"Fixed_{key}_sd_value"
                            )

                submitted = st.form_submit_button("Run flux calculation")
            return submitted

    def _get_data_from_session_state(self):
        """
        Get the data from the widgets input
        """

        # Start with estimable parameters
        # Get order of parameter names to build dict
        estimable_parameter_name_order = [key for key in self.model.parameters_to_estimate.keys()]
        for name in estimable_parameter_name_order:
            try:
                # Get values from widgets
                if st.session_state[f"Parameter_value_{name}"] == "0":
                    self.model.parameters_to_estimate[name] = 0
                else:
                    self.model.parameters_to_estimate[name] = literal_eval(
                        st.session_state[f"Parameter_value_{name}"].lstrip("0") # Strip leading zeroes to stop eval errors
                    )
                # Get bounds
                if st.session_state[f"Parameter_lower_{name}"] == "0":
                    st.warning(
                        f"WARNING: {name} has a lower bound at 0. Sometimes this might confuse the optimizer. It is "
                        f"advised to replace 0 with a very small number, 1e-6 for example."
                    )
                    lower_bound = 0
                else:
                    lower_bound = literal_eval(st.session_state[f"Parameter_lower_{name}"].lstrip("0"))
                if st.session_state[f"Parameter_upper_{name}"] == "0":
                    st.warning(
                        f"WARNING: {name} has an upper bound at 0. Sometimes this might confuse the optimizer. It is "
                        f"advised to replace 0 with a very small number, 1e-6 for example."
                    )
                    upper_bound = 0
                else:
                    upper_bound = literal_eval(st.session_state[f"Parameter_upper_{name}"].lstrip("0"))
                self.model.bounds[name] = (
                    lower_bound,
                    upper_bound
                )
            except ValueError:
                st.error(
                    f"ERROR: An error occurred while parsing the input for {name}. Please check that only numbers "
                    f"have been entered."
                )
                raise

        # Do the same for each fixed parameter class
        if self.model.fixed_parameters is not None:
            fixed_parameter_classes = [param for param in self.model.fixed_parameters]
            for param in fixed_parameter_classes:
                for key in self.model.fixed_parameters[param].keys():
                    try:
                        if st.session_state[f"Fixed_{param}_value_{key}"] == "0":
                            self.model.fixed_parameters[param][key] = 0
                        else:
                            self.model.fixed_parameters[param][key] = literal_eval(
                                st.session_state[f"Fixed_{param}_value_{key}"].lstrip("0")
                            )
                    except ValueError:
                        st.error(
                            f"ERROR: An error occurred while parsing the input in the fixed parameter class {param} for "
                            f"{key}. Please check that only numbers have been entered."
                        )
                        raise

        # And finally do the same for Standard Deviations
        self.sd = self.defaults["sd"]
        sd_name_order = [key for key in self.sd.keys()]
        for name in sd_name_order:
            try:
                # Get values from widgets
                self.sd[name] = literal_eval(
                    st.session_state[f"Fixed_{name}_sd_value"].lstrip("0")  # Strip leading zeroes to stop eval errors
                )
            except ValueError:
                st.error(
                    f"ERROR: An error occurred while parsing the input for {name} (SDs). Please check that only numbers "
                    f"have been entered."
                )
                raise

    def _build_fitter_kwargs(self):

        kwargs = {
            "sd": self.sd,
            "model": self.model,
            "mc": self.mc,
            "iterations": self.iterations,
            "debug_mode": self.debug_mode,
        }
        return kwargs

    def _build_internal_json(self, path_to_data):

        final_json = json.dumps({
            "sd": self.sd,
            "model": self.model.model_name,
            "mc": self.mc,
            "iterations": self.iterations,
            "debug_mode": self.debug_mode,
            "path_to_data": path_to_data
        })
        return final_json

    def _output_directory_selector(self):

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
            self.io_handler.res_path = self.io_handler.home_path / (self.io_handler.home_path.name + "_res")
            if not clicked:
                if not self.io_handler.res_path.is_dir():
                    self.io_handler.res_path.mkdir()



if __name__ == "__main__":
    physiofit_app = App()
    physiofit_app.start_app()
