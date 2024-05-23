from ast import literal_eval
import tkinter as tk
from tkinter import filedialog
from pathlib import Path
from copy import copy
import logging

import pandas as pd
import streamlit as st

import physiofit
from physiofit.base.io import IoHandler, ConfigParser
from physiofit.models.base_model import StandardDevs

logger = logging.getLogger("physiofit")
logger.setLevel(logging.DEBUG)


class App:
    """
    Physiofit Graphical User Interface
    """

    def __init__(self):

        self.defaults = {
            "iterations": 100,
            "sd": StandardDevs(),
            "mc": True
        }
        self.select_menu = None
        self.io = None
        self.update_info = None
        self.config_parser = None

    def start_app(self):
        """Launch the application"""

        st.set_page_config(page_title=f"PhysioFit (v{physiofit.__version__})")
        st.title(f"Welcome to PhysioFit (v{physiofit.__version__})")
        st.write(
            "Documentation available at [https://physiofit.readthedocs.io]("
            "https://physiofit.readthedocs.io).")
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

    @staticmethod
    def _get_data_path_config_context():

        # Set up tkinter for directory chooser
        root = tk.Tk()
        root.withdraw()

        # Make folder picker dialog appear on top of other windows
        root.wm_attributes('-topmost', 1)

        return str(Path(
            filedialog.askopenfilename(
                master=root,
                title="Select the data file",
                filetypes=[("Data files", "*.tsv *.txt")]
            )
        ))

    @staticmethod
    def clear_input_from_session_state():

        if hasattr(st.session_state, "config_parser_data_path"):
            del st.session_state["config_parser_data_path"]

    def _run(self):
        """
        Initialize all the options. If a yaml file is given as input, options
        are parsed from it. If the input is tsv, the defaults are used instead.
        """

        # Check file extension and initialize defaults
        file_extension = self.data_file.name.split(".")[1]
        if file_extension in ["yaml", "yml"]:
            # We store the path in session state to
            # avoid multiple prompts.
            try:
                # For yaml files, we need to read the configuration file
                # before loading the data.
                self.config_parser = self.io.read_yaml(self.data_file)
                # Check if the data path exists, if not open up prompt to
                # select file
                if not self.config_parser.check_data_path():
                    st.info(
                        "File path in configuration file is incorrect or "
                        "missing. Loading file from prompt."
                    )
                    input_path = st.session_state.get(
                        "config_parser_data_path",
                        self._get_data_path_config_context()
                    )
                    self.config_parser.path_to_data = input_path
                    st.session_state["config_parser_data_path"] = input_path

                st.info(f"Input data: {self.config_parser.path_to_data}")
                self.input_datafile_name = Path(
                    self.config_parser.path_to_data).stem
                # Load data into io_handler
                self.io.data = self.io.read_data(
                    str(self.config_parser.path_to_data)
                )
            except Exception:
                st.error(
                    "An error has occurred when reading the yaml "
                    "configuration file and/or loading the input data."
                )
                raise
        elif file_extension in ["tsv", "txt"]:
            try:
                self.input_datafile_name = self.data_file.name[:-4]
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

        # Reset config_parser path to data if it exists to handle new
        # data file paths
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
            st.error(
                f"An error has occurred when listing models from the "
                f"models folder: \n{Path(__file__).parent / 'models'}. "
                f"Please correct the model or submit an issue at "
                f"github.com/MetaSys-LISBP/PhysioFit/issues")
            raise

        # Build menu
        submitted = self._initialize_opt_menu_widgets()

        if submitted:
            _logger = self._build_logger(self.io.res_path)
            try:
                self._get_data_from_session_state()
            except Exception:
                st.error(
                    "An error has occurred when initializing the model")
                raise
            if not self.io.wkdir:
                raise ValueError("No output directory selected")

            _logger.info("=============================================")
            _logger.info(f"Physiofit version: {physiofit.__version__}")
            _logger.info("Path to results: \n" + str(self.io.res_path))
            _logger.info("=============================================")

            self.config_parser = ConfigParser(
                path_to_data=self.io.wkdir / self.data_file.name,
                selected_model=self.model,
                sds=self.sd,
                mc=self.mc,
                iterations=self.iterations
            )

            full_dataframe = self.io.data.copy()
            results_path = copy(self.io.res_path)
            experiments = list(self.io.data["experiments"].unique())
            self.io.multiple_experiments = []
            for experiment in experiments:
                _logger.info(f"Running optimization for {experiment}")
                with st.spinner(f"Running optimization for {experiment}"):
                    # final_table_dict = {}
                    self.model.data = full_dataframe[
                        full_dataframe["experiments"] == experiment
                        ].drop("experiments", axis=1).copy()

                    self.io.res_path = results_path / str(experiment)
                    if not self.io.res_path.is_dir():
                        self.io.res_path.mkdir(parents=True)
                    # Initialize the fitter object
                    self.io.names = self.io.data.columns[1:].to_list()
                    kwargs = self._build_fitter_kwargs()
                    _logger.info("Run options for the fitter:")
                    for key, value in kwargs.items():
                        _logger.info(f"{key} : {value}")
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
                    printed_df = df.copy()
                    df.index = [
                        f"{experiment} {param}" for param in
                        fitter.model.parameters.keys()
                    ]
                    printed_df.index = [
                        f"{param}" for param in fitter.model.parameters.keys()
                    ]
                    _logger.info(f"Results for {experiment}: \n{df}")
                    self.io.multiple_experiments.append(df)

                    # Export results
                    self.io.output_report(fitter, self.io.res_path)
                    self.io.plot_data(fitter)
                    self.io.output_plots(fitter, self.io.res_path)
                    with st.expander(f"{experiment} plots"):
                        data_col, stat_col = st.columns(2)
                        with data_col:
                            st.write("Parameter statistics:")
                            st.dataframe(printed_df)
                        with stat_col:
                            st.write("Khi² test results:")
                            st.dataframe(fitter.khi2_res)
                        if fitter.khi2_res.at["p_val", "Values"] < 0.95:
                            st.write(
                                f"\n\nAt level of 95% confidence, "
                                f"the model fits the data good enough "
                                f"with respect to the provided "
                                f"measurement SD. Value: "
                                f"{fitter.khi2_res.at['p_val', 'Values']}"
                            )

                        else:
                            st.write(
                                f"\n\nAt level of 95% confidence, "
                                f"the model does not fit the data good "
                                f"enough with respect to the provided "
                                f"measurement SD. Value: "
                                f"{fitter.khi2_res.at['p_val', 'Values']}"
                            )
                        for fig in self.io.figures:
                            st.pyplot(fig[1])
                    self.io.output_pdf(fitter, self.io.res_path)
                    # Reset figures to free memory
                    self.io.figures = []
                    self.config_parser.export_config(self.io.res_path)
            self.io.data = full_dataframe
            _logger.info(f"Resulting dataframe: \n{full_dataframe}")
            self.io.res_path = results_path
            self.io.output_recap(results_path)
            logging.shutdown()

    def _build_flux_menu(self):
        """Build the starting menu with the data upload button"""

        self.data_file = st.file_uploader(
            "Load a data file or a json configuration file",
            key="data_uploader",
            accept_multiple_files=False,
            on_change=self.clear_input_from_session_state
        )

        if self.data_file:
            self._run()

    def silent_sim(self):

        self.model.simulate(
            [param for param in self.model.parameters.values()],
            self.model.experimental_matrix,
            self.model.time_vector,
            self.model.args
        )

    def _build_logger(self, output_path):
        """
        Sets up a logger for the application.

        This method creates a logging handler that writes to a file and a
        stream handler that writes to the console.
        The log level is set to INFO by default, but can be set to DEBUG if
        the debug_mode attribute is True.
        A formatter is also set up to format the log messages.

        :param output_path: The path where the log file will be written.
        :type output_path: Path
        :return: The configured logger.
        :rtype: Logger
        """
        handler = logging.FileHandler(output_path / "log.txt",
                                      "w")
        stream = logging.StreamHandler()
        handler.setLevel(logging.INFO)
        stream.setLevel(logging.INFO)
        if self.debug_mode:
            handler.setLevel(logging.DEBUG)
            stream.setLevel(logging.DEBUG)
        formatter = logging.Formatter('%(asctime)s - %(name)s - %('
                                      'levelname)s - %(message)s')
        handler.setFormatter(formatter)
        if not logger.hasHandlers():
            logger.addHandler(handler)
            logger.addHandler(stream)
        return logger

    def _initialize_opt_menu_widgets(self):

        # Get model names and give as options to user
        model_options = [
            model.name for model in self.io.models
        ]
        idx = None
        if self.config_parser:
            if self.config_parser.model:
                try:
                    idx = model_options.index(
                        self.config_parser.model["model_name"])
                except Exception:
                    st.error(
                        "Error while reading model name from configuration "
                        "file"
                    )
                    raise

        model_name = st.selectbox(
            label="Model",
            options=model_options,
            key="model_selector",
            help="Select the flux calculation model",
            index=idx
        )

        # if model_name == "Dynamic system (only substrates)":
        #    st.error("Not yet implemented...")
        if model_name:
            # Initialize selected model
            self.model = self.io.select_model(model_name)
            try:
                self.model.get_params()
                self.silent_sim()
            except Exception as e:
                st.error(
                    f"The following error occurred with the selected model: "
                    f"{e}")
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
                    value=self.defaults[
                        "iterations"] if self.config_parser is None
                    else self.config_parser.iterations,
                    help="Number of iterations for Monte Carlo analysis.",
                    disabled=enable_mc
                )
                if self.iterations < 0:
                    st.error(
                        "ERROR: Number of Monte-Carlo iterations must be a "
                        "positive integer"
                    )
                self.debug_mode = st.checkbox(
                    "Verbose logs",
                    help="Useful in case of trouble. Join it to the "
                         "issue on github."
                )

                # if file_extension in ["tsv", "txt"]:

                self._output_directory_selector()
                self.io.res_path = self.io.wkdir / (
                        self.input_datafile_name + "_res")
                self.io.res_path.mkdir(parents=True, exist_ok=True)
            # Build the form for advanced parameters
            form = st.form("Parameter_form")
            with form:
                expand_run_params = st.expander("Parameters")
                with expand_run_params:
                    st.subheader("Parameters to estimate")
                    col1, col2, col3, col4 = st.columns(4)
                    with col1:
                        st.write("Parameter Name")
                        for key in self.model.parameters:
                            st.text_input(
                                label="label",  # Unused
                                label_visibility="collapsed",
                                value=key,
                                key=f"Parameter_name_{key}",
                                disabled=True
                            )
                    with col2:
                        st.write("Parameter Value")
                        for key, value in self.model.parameters.items():
                            st.text_input(
                                label="label",  # Unused
                                label_visibility="collapsed",
                                value=value if self.config_parser is None
                                else self.config_parser.model[
                                    "parameters_to_estimate"][key],
                                key=f"Parameter_value_{key}"
                            )
                    with col3:
                        st.write("Lower Bound")
                        for key, bound in self.model.bounds.items():
                            st.text_input(
                                label="label",  # Unused
                                label_visibility="collapsed",
                                value=bound[0] if self.config_parser is None
                                else literal_eval(
                                    self.config_parser.model["bounds"][key])[
                                    0],
                                key=f"Parameter_lower_{key}"
                            )
                    with col4:
                        st.write("Upper Bound")
                        for key, bound in self.model.bounds.items():
                            st.text_input(
                                label="label",  # Unused
                                label_visibility="collapsed",
                                value=bound[1] if self.config_parser is None
                                else literal_eval(
                                    self.config_parser.model["bounds"][key])[
                                    1],
                                key=f"Parameter_upper_{key}"
                            )

                    if self.model.args is not None:
                        for param in self.model.args.keys():
                            st.subheader(f"Fixed parameters: {param}")
                            col1, col2 = st.columns(2)
                            with col1:
                                st.write("Parameter Name")
                                for key in self.model.args[param].keys():
                                    st.text_input(
                                        label="label",  # Unused
                                        label_visibility="collapsed",
                                        value=key,
                                        key=f"Fixed_{param}_{key}",
                                        disabled=True
                                    )
                            with col2:
                                st.write("Parameter Value")
                                for key, value in self.model.args[
                                     param].items():
                                    st.text_input(
                                        label="label",  # Unused
                                        label_visibility="collapsed",
                                        value=value if self.config_parser is
                                        None else
                                        self.config_parser.model["args"][key],
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
        estimable_parameter_name_order = [key for key in
                                          self.model.parameters.keys()]
        for name in estimable_parameter_name_order:
            try:
                # Get values from widgets
                if st.session_state[f"Parameter_value_{name}"] == "0":
                    self.model.parameters[name] = 0
                else:
                    self.model.parameters[name] = literal_eval(
                        st.session_state[f"Parameter_value_{name}"].lstrip("0")
                        # Strip leading zeroes to stop eval errors
                    )
                # Get bounds
                if st.session_state[f"Parameter_lower_{name}"] == "0":
                    st.warning(
                        f"WARNING: {name} has a lower bound at 0. Sometimes "
                        f"this might confuse the optimizer. We strongly "
                        f"recommend to set the lower bound at a non-zero "
                        f"value, 1e-6 for example."
                    )
                    lower_bound = 0
                else:
                    lower_bound = literal_eval(
                        st.session_state[f"Parameter_lower_{name}"].lstrip(
                            "0"))
                if st.session_state[f"Parameter_upper_{name}"] == "0":
                    st.warning(
                        f"WARNING: {name} has an upper bound at 0. Sometimes "
                        f"this might confuse the optimizer. We strongly "
                        f"recommend to set the lower bound at a non-zero "
                        f"value, 1e-6 for example."
                    )
                    upper_bound = 0
                else:
                    upper_bound = literal_eval(
                        st.session_state[f"Parameter_upper_{name}"].lstrip(
                            "0"))
                self.model.bounds[name] = (
                    lower_bound,
                    upper_bound
                )
            except ValueError:
                st.error(
                    f"ERROR: An error occurred while parsing the input for"
                    f" {name}. Please check that only numbers have been "
                    f"entered."
                )
                raise

        # Do the same for each fixed parameter class
        if self.model.args is not None:
            fixed_parameter_classes = [param for param in self.model.args]
            for param in fixed_parameter_classes:
                for key in self.model.args[param].keys():
                    try:
                        if st.session_state[
                                f"Fixed_{param}_value_{key}"] == "0":
                            self.model.args[param][key] = 0
                        else:
                            self.model.args[param][key] = literal_eval(
                                st.session_state[
                                    f"Fixed_{param}_value_{key}"].lstrip("0")
                            )
                    except ValueError:
                        st.error(
                            f"ERROR: An error occurred when parsing the input "
                            f"in the fixed parameter class {param} for {key}. "
                            f"Please check that only numbers have been "
                            f"entered."
                        )
                        raise

        # And finally do the same for Standard Deviations
        sd_name_order = [key for key in self.sd.keys()]
        for name in sd_name_order:
            try:
                if st.session_state[f"Fixed_{name}_sd_value"] == "0":
                    self.sd[name] = 0  # will raise Value Error as expected
                else:
                    # Get values from widgets
                    self.sd[name] = literal_eval(
                        st.session_state[f"Fixed_{name}_sd_value"].lstrip("0")
                        # Strip leading zeroes to stop eval errors
                    )
            except ValueError:
                st.error(
                    f"ERROR: An error occurred when parsing the input for "
                    f"{name} (SDs). Please check that only numbers have been "
                    f"entered and that value is superior to 0."
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

        # Initialize folder picker button and add logic
        clicked = st.button(
            "Select output data directory", key="clicker"
        )
        if clicked:

            # Set up tkinter for directory chooser
            root = tk.Tk()
            root.withdraw()

            # Make folder picker dialog appear on top of other windows
            root.wm_attributes('-topmost', 1)

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
