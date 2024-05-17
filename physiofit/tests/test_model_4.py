import logging

import numpy as np
import pandas as pd
import pytest

import physiofit

logging.getLogger("physiofit").setLevel(logging.ERROR)


@pytest.fixture
def model_4_data():

    return pd.DataFrame(
        {'time': {0: 0.0, 1: 0.2, 2: 0.4, 3: 0.6000000000000001, 4: 0.8,
                  5: 1.0, 6: 1.2000000000000002, 7: 1.4000000000000001, 8: 1.6,
                  9: 1.8, 10: 2.0, 11: 2.2, 12: 2.4000000000000004, 13: 2.6,
                  14: 2.8000000000000003, 15: 3.0, 16: 3.2,
                  17: 3.4000000000000004, 18: 3.6, 19: 3.8000000000000003,
                  20: 4.0, 21: 4.2, 22: 4.4, 23: 4.6000000000000005,
                  24: 4.800000000000001, 25: 5.0, 26: 5.2, 27: 5.4,
                  28: 5.6000000000000005, 29: 5.800000000000001},
         'X': {0: 0.02, 1: 0.023470217419836206, 2: 0.027542555286719145,
               3: 0.03232148804385787, 4: 0.037929617586099036,
               5: 0.044510818569849356, 6: 0.052233929468462365,
               7: 0.06129708406586005, 8: 0.07193279451138565,
               9: 0.08441391633993106, 10: 0.0990606484879023,
               11: 0.1162487478880518, 12: 0.13641916938581505,
               13: 0.16008937828592706, 14: 0.18786662574885568,
               15: 0.2204635276128321, 16: 0.25871634631086166,
               17: 0.3036064448990781, 18: 0.3562854635922441,
               19: 0.4181048647018553, 20: 0.49065060394218707,
               21: 0.5757838175848538, 22: 0.6756885692769915,
               23: 0.7929278814514524, 24: 0.9305094887957849,
               25: 1.0919630006628847, 26: 1.281430451998733,
               27: 1.5037725658404621, 28: 1.7646934535130303,
               29: 2.070886951665622},
         'Glucose': {0: 20.0, 1: 19.96529782580164, 2: 19.92457444713281,
                     3: 19.876785119561422, 4: 19.82070382413901,
                     5: 19.754891814301505, 6: 19.677660705315375,
                     7: 19.587029159341398, 8: 19.480672054886142,
                     9: 19.35586083660069, 10: 19.209393515120976,
                     11: 19.03751252111948, 12: 18.83580830614185,
                     13: 18.59910621714073, 14: 18.32133374251144,
                     15: 17.99536472387168, 16: 17.612836536891383,
                     17: 17.163935551009217, 18: 16.63714536407756,
                     19: 16.018951352981446, 20: 15.293493960578129,
                     21: 14.442161824151462, 22: 13.443114307230086,
                     23: 12.270721185485478, 24: 10.894905112042153,
                     25: 9.280369993371155, 26: 7.385695480012673,
                     27: 5.16227434159538, 28: 2.5530654648697,
                     29: -0.5088695166562189},
         'Acetate': {0: 0.01, 1: 0.02301331532438577, 2: 0.03828458232519679,
                     3: 0.05620558016446701, 4: 0.07723606594787137,
                     5: 0.10191556963693509, 6: 0.13087723550673386,
                     7: 0.1648640652469752, 8: 0.2047479794176962,
                     9: 0.25155218627474146, 10: 0.3064774318296336,
                     11: 0.37093280458019423, 12: 0.4465718851968064,
                     13: 0.5353351685722265, 14: 0.6394998465582088,
                     15: 0.7617382285481203, 16: 0.9051862986657312,
                     17: 1.0735241683715429, 18: 1.2710704884709152,
                     19: 1.5028932426319574, 20: 1.7749397647832015,
                     21: 2.094189315943201, 22: 2.4688321347887174,
                     23: 2.908479555442946, 24: 3.424410582984193,
                     25: 4.0298612524858175, 26: 4.740364194995248,
                     27: 5.574147121901733, 28: 6.552600450673863,
                     29: 7.700826068746083},
         'Glutamate': {0: 0.01, 1: 0.018675543549590515,
                       2: 0.02885638821679786, 3: 0.04080372010964467,
                       4: 0.05482404396524758, 5: 0.07127704642462339,
                       6: 0.0905848236711559, 7: 0.11324271016465011,
                       8: 0.13983198627846413, 9: 0.17103479084982764,
                       10: 0.20765162121975575, 11: 0.25062186972012945,
                       12: 0.30104792346453757, 13: 0.3602234457148176,
                       14: 0.4296665643721392, 15: 0.5111588190320802,
                       16: 0.6067908657771541, 17: 0.7190161122476952,
                       18: 0.8507136589806101, 19: 1.0052621617546382,
                       20: 1.1866265098554676, 21: 1.3994595439621342,
                       22: 1.6492214231924784, 23: 1.9423197036286308,
                       24: 2.2862737219894615, 25: 2.689907501657211,
                       26: 3.1635761299968315, 27: 3.719431414601155,
                       28: 4.371733633782575, 29: 5.1372173791640545}}
    )


def test_model_4_estimation(
        model_4_data: pd.DataFrame,
        sds: physiofit.models.base_model.StandardDevs
):
    io = physiofit.base.io.IoHandler()
    model = io.select_model(
        name="Steady-state batch model",
        data=model_4_data
    )
    model.get_params()
    fitter = io.initialize_fitter(
        data=model.data,
        model=model,
        sd=sds,
        debug_mode=False
    )
    fitter.optimize()

    assert np.allclose(
        a=fitter.parameter_stats["optimal"],
        b=[0.02, 0.8, -8, 20, 3, 0.01, 2, 0.01],
        rtol=1e-3
    )


def test_model_4_simulation(
        placeholder_data,
        parameters,
        model_4_data
):
    io = physiofit.base.io.IoHandler()
    model = io.select_model(
        name="Steady-state batch model",
        data=placeholder_data
    )
    model.get_params()
    model.parameters.update(parameters)
    sim_data = model.simulate(
        list(model.parameters.values()),
        model.data.drop("time", axis=1),
        model.time_vector,
        model.args
    )
    df = pd.DataFrame(
        data=sim_data,
        index=model.time_vector,
        columns=model.name_vector
    )
    df.index.name = "time"

    pd.testing.assert_frame_equal(
        df.reset_index(),
        model_4_data,
        atol=1e-6
    )
