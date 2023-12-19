from ast import literal_eval
import tkinter as tk
from tkinter import filedialog
from pathlib import Path
from copy import copy

import pandas as pd
import streamlit as st

import physiofit
from physiofit.base.io import IoHandler, ConfigParser
from physiofit.models.base_model import StandardDevs


class App:
    """
    Physiofit Graphical User Interface
    """

    def __init__(self):

        self.defaults = {
            "iterations" : 100,
            "sd" : StandardDevs(),
            "mc" : True
        }
        self.select_menu = None
        self.io = None
        self.update_info = None
        self.config_parser = None

    def start_app(self):
        """Launch the application"""

        st.set_page_config(page_title=f"PhysioFit (v{physiofit.__version__})")
        st.title(f"Welcome to PhysioFit (v{physiofit.__version__})")
        st.write("Documentation available at [https://physiofit.readthedocs.io](https://physiofit.readthedocs.io).")
        self.update_info = st.empty()
        self.check_uptodate()
        self.select_menu = st.selectbox(
            "Select a task to execute",
            ["Calculate extracellular fluxes",
             "Calculate degradation constant"]
        )
        if self.select_menu == "Calculate extracellular fluxes":
            self.io = IoHandler()
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
        Initialize all the options. If a yaml file is given as input, options
        are parsed from it. If the input is tsv, the defaults are used instead.
        """

        # Check file extension and initialize defaults
        file_extension = self.data_file.name.split(".")[1]
        if file_extension in ["yaml", "yml"]:
            try:
                # Get parameters from yaml file
                self.config_parser = self.io.read_yaml(self.data_file)
                # Load data into io_handler
                self.io.data = self.io.read_data(self.config_parser.path_to_data)
            except Exception:
                st.error("An error has occurred when reading the yaml configuration file.")
                raise
        elif file_extension in ["tsv", "txt"]:
            try:
                self.io.data = self.io.read_data(self.data_file)
                # Initialize default SDs
                self.defaults["sd"].update({
                    "X": 0.5
                })
                self.defaults["sd"].update({
                    col: 0.2 for col in self.io.data.columns[2:]
                })
            except Exception:
                st.error("An error has occurred when reading the data.")
                raise
        else:
            raise KeyError(
                f"Wrong input file format. Accepted formats are tsv for data "
                f"files or json for configuration files. Detected file: "
                f"{self.data_file.name}")

        if "experiments" not in self.io.data.columns:
            raise ValueError(
                "'experiments' column missing from dataset"
            )
        self.io.data = self.io.data.sort_values(
            ["experiments", "time"], ignore_index=True
        )

        try:
            # Initialize the list of available models
            self.io.get_models()
        except Exception:
            st.error(f"An error has occurred when listing models from the models folder: "
                     f"\n{Path(__file__).parent / 'models'}. Please correct the model or submit an issue at "
                     f"github.com/MetaSys-LISBP/PhysioFit/issues")
            raise

        # Build menu
        submitted = self._initialize_opt_menu_widgets(
            file_extension
        )

        if submitted:
            try:
                self._get_data_from_session_state()
            except Exception:
                st.error("An error has occurred when initializing the model")
                raise
            if not self.io.wkdir:
                raise ValueError("No output directory selected")
            self.config_parser = ConfigParser(
                path_to_data =self.io.wkdir / self.data_file.name,
                selected_model= self.model,
                sds = self.sd,
                mc = self.mc,
                iterations = self.iterations
            )

            full_dataframe = self.io.data.copy()
            results_path = copy(self.io.res_path)
            experiments = list(self.io.data["experiments"].unique())
            self.io.multiple_experiments = []
            for experiment in experiments:
                with st.spinner(f"Running optimization for {experiment}"):
                    # final_table_dict = {}
                    self.model.data = full_dataframe[
                        full_dataframe["experiments"] == experiment
                    ].drop("experiments", axis=1).copy()

                    self.io.res_path = results_path / experiment
                    if not self.io.res_path.is_dir():
                        self.io.res_path.mkdir(parents=True)
                    # Initialize the fitter object
                    self.io.names = self.io.data.columns[1:].to_list()
                    kwargs = self._build_fitter_kwargs()
                    fitter = self.io.initialize_fitter(
                        self.model.data,
                        model=kwargs["model"],
                        mc=kwargs["mc"],
                        iterations=kwargs["iterations"],
                        sd=kwargs["sd"],
                        debug_mode=kwargs["debug_mode"]
                    )
                    # Do the work
                    fitter.optimize()
                    if self.mc:
                        fitter.monte_carlo_analysis()
                    fitter.khi2_test()
                    df = pd.DataFrame.from_dict(
                        fitter.parameter_stats,
                        orient="columns"
                    )
                    df.index = [
                        f"{experiment} {param}" for param in fitter.model.parameters_to_estimate.keys()
                    ]
                    st.write(df)
                    self.io.multiple_experiments.append(df)

                    # Export results
                    self.io.output_report(fitter, self.io.res_path)
                    self.io.plot_data(fitter)
                    self.io.output_plots(fitter, self.io.res_path)
                    with st.expander(f"{experiment} plots"):
                        for fig in self.io.figures:
                            st.pyplot(fig[1])
                    self.io.output_pdf(fitter, self.io.res_path)
                    # Reset figures to free memory
                    self.io.figures = []
                    self.config_parser.export_config(self.io.res_path)
            self.io.data = full_dataframe
            self.io.res_path = results_path
            self.io.output_recap(results_path)

    def silent_sim(self):

        self.model.simulate(
            [param for param in self.model.parameters_to_estimate.values()],
            self.model.experimental_matrix,
            self.model.time_vector,
            self.model.fixed_parameters
        )

    def _initialize_opt_menu_widgets(self, file_extension):

        # Get model names and give as options to user
        model_options = [
            model.model_name for model in self.io.models
        ]
        if self.config_parser:
            if self.config_parser.model:
                try:
                    idx = model_options.index(self.config_parser.model["model_name"])
                except Exception:
                    st.error("Error while reading model name from configuration file")
                    raise
        else:
            idx = None
        model_name = st.selectbox(
            label="Model",
            options=model_options,
            key="model_selector",
            help="Select the flux calculation model",
            index=idx
        )

        #if model_name == "Dynamic system (only substrates)":
        #    st.error("Not yet implemented...")
        if model_name:
            # Initialize selected model
            self.model = self.io.select_model(model_name)
            try:
                self.model.get_params()
                self.silent_sim()
            except Exception as e:
                st.error(f"The following error occurred with the selected model: {e}")
                raise
            else:
                st.success("Model successfully initialized")

            expand_basic_settings = st.expander("Basic settings",
                                                expanded=True)
            with expand_basic_settings:
                # Select model parameters
                self.mc = st.checkbox(
                    "Sensitivity analysis (Monte Carlo)",
                    value=self.defaults["mc"] if self.config_parser is None
                    else self.config_parser.mc,
                    help="Determine the precision on estimated parameters by "
                         "Monte Carlo sensitivity analysis."
                )
                enable_mc = False if self.mc else True
                self.iterations = st.number_input(
                    "Number of iterations",
                    value=self.defaults["iterations"] if self.config_parser is None
                    else self.config_parser.iterations ,
                    help="Number of iterations for Monte Carlo analysis.",
                    disabled=enable_mc
                )
                if self.iterations < 0:
                    st.error("ERROR: Number of Monte-Carlo iterations must be a positive integer")
                self.debug_mode = st.checkbox(
                    "Verbose logs",
                    help="Useful in case of trouble. Join it to the "
                         "issue on github."
                )

                if file_extension in ["tsv", "txt"]:

                    self._output_directory_selector()

                else:
                    self.io.wkdir = Path(self.config_parser.path_to_data).resolve().parent
                    self.io.res_path = self.io.wkdir / (self.io.wkdir.name + "_res")

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
                                value=value if self.config_parser is None
                                else self.config_parser.model["parameters_to_estimate"][key],
                                key=f"Parameter_value_{key}"
                            )
                    with col3:
                        st.write("Lower Bound")
                        for key, bound in self.model.bounds.items():
                            st.text_input(
                                label="label",  # Unused
                                label_visibility="collapsed",
                                value=bound[0] if self.config_parser is None
                                else literal_eval(self.config_parser.model["bounds"][key])[0],
                                key=f"Parameter_lower_{key}"
                            )
                    with col4:
                        st.write("Upper Bound")
                        for key, bound in self.model.bounds.items():
                            st.text_input(
                                label="label",  # Unused
                                label_visibility="collapsed",
                                value=bound[1] if self.config_parser is None
                                else literal_eval(self.config_parser.model["bounds"][key])[1],
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
                                        value=value if self.config_parser is None
                                        else self.config_parser.model["fixed_parameters"][key],
                                        key=f"Fixed_{param}_value_{key}"
                                    )

                expand_sds = st.expander("Standard Deviations")
                # Get origin of sds
                if self.config_parser is None:
                    self.sd = self.defaults["sd"]
                else:
                    self.sd = self.config_parser.sds
                with expand_sds:
                    col1, col2 = st.columns(2)
                    with col1:
                        for key in self.sd.keys():
                            st.text_input(
                                label="label",  # Unused
                                label_visibility="collapsed",
                                value=key,
                                key=f"{key}_sd",
                                disabled=True
                            )
                    with col2:
                        for key, value in self.sd.items():
                            st.text_input(
                                label="label",  # Unused
                                label_visibility="collapsed",
                                value=value if self.config_parser is None
                                else self.config_parser.sds[key],
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
                        f"WARNING: {name} has a lower bound at 0. Sometimes this might confuse the optimizer. We "
                        f"strongly recommend to set the lower bound at a non-zero value, 1e-6 for example."
                    )
                    lower_bound = 0
                else:
                    lower_bound = literal_eval(st.session_state[f"Parameter_lower_{name}"].lstrip("0"))
                if st.session_state[f"Parameter_upper_{name}"] == "0":
                    st.warning(
                        f"WARNING: {name} has an upper bound at 0. Sometimes this might confuse the optimizer. We "
                        f"strongly recommend to set the lower bound at a non-zero value, 1e-6 for example."
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
                            f"ERROR: An error occurred when parsing the input in the fixed parameter class {param} for "
                            f"{key}. Please check that only numbers have been entered."
                        )
                        raise

        # And finally do the same for Standard Deviations
        sd_name_order = [key for key in self.sd.keys()]
        for name in sd_name_order:
            try:
                if st.session_state[f"Fixed_{name}_sd_value"] == "0":
                    self.sd[name] = 0 # will raise Value Error as expected
                else:
                    # Get values from widgets
                    self.sd[name] = literal_eval(
                        st.session_state[f"Fixed_{name}_sd_value"].lstrip("0")  # Strip leading zeroes to stop eval errors
                    )
            except ValueError:
                st.error(
                    f"ERROR: An error occurred when parsing the input for {name} (SDs). Please check that only numbers "
                    f"have been entered and that value is superior to 0."
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
            st.session_state.wkdir = Path(st.text_input(
                "Selected output data directory:",
                filedialog.askdirectory(master=root)
            ))
            if st.session_state.wkdir == Path(".") \
                    or not st.session_state.wkdir.exists():
                raise RuntimeError("Please provide a valid output directory")
            self.io.wkdir = copy(
                st.session_state.wkdir)

        elif hasattr(st.session_state, "wkdir"):

            self.io.wkdir = Path(st.text_input(
                "Selected output data directory:",
                st.session_state.wkdir
            ))
            if self.io.wkdir == Path(".") \
                    or not self.io.wkdir.exists():
                raise RuntimeError("Please provide a valid output directory")

            # Initialize the result export directory
            self.io.res_path = self.io.wkdir / (self.io.wkdir.name + "_res")
            if not clicked:
                if not self.io.res_path.is_dir():
                    self.io.res_path.mkdir()



if __name__ == "__main__":
    physiofit_app = App()
    physiofit_app.start_app()
